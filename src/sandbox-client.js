// sandbox-client.js - Robust E2B wrapper with dual timeouts and health checks
const { Sandbox } = require('@e2b/code-interpreter');

const TEMPLATE = process.env.E2B_TEMPLATE_ID;
const API_KEY = process.env.E2B_API_KEY;

let sbx = null;
let lastTouch = 0;

async function ensureSandbox({ forceNew = false } = {}) {
  const now = Date.now();
  const expired = !sbx || (now - lastTouch > 55 * 60 * 1000); // Hobby kills at 60m

  if (forceNew || expired || sbx?.closed) {
    if (sbx && !sbx.closed) { 
      try { await sbx.kill(); } catch {} 
    }
    sbx = await Sandbox.create({ 
      templateId: TEMPLATE, 
      apiKey: API_KEY, 
      timeout: 3_600_000 
    });
  }

  // Liveness ping (fast)
  try {
    const { exitCode } = await sbx.runCode('print("__pong__")', { 
      language: 'python', 
      timeout: 5000 
    });
    if (exitCode !== 0) throw new Error('ping failed');
  } catch {
    // recreate once
    sbx = await Sandbox.create({ 
      templateId: TEMPLATE, 
      apiKey: API_KEY, 
      timeout: 3_600_000 
    });
  }

  lastTouch = Date.now();
  return sbx;
}

async function runPythonSafe(source, { 
  runTimeoutMs = 30000, 
  overallDeadlineMs = 35000 
} = {}) {
  const controller = new AbortController();
  const overallTimer = setTimeout(() => controller.abort(), overallDeadlineMs);

  const code = `
import signal, sys
def _alarm_handler(signum,frame): 
    print("__TIMEOUT__", file=sys.stderr); 
    sys.exit(124)
signal.signal(signal.SIGALRM, _alarm_handler)
signal.alarm(${Math.max(1, Math.floor(runTimeoutMs/1000))})
# --- user code starts ---
${source}
# --- user code ends ---
`;

  try {
    const sandbox = await ensureSandbox();
    const res = await sandbox.runCode(code, { 
      language: 'python', 
      timeout: runTimeoutMs, 
      signal: controller.signal 
    });

    lastTouch = Date.now();
    // Normalize outputs (arrays → string)
    const out = (res.logs?.stdout || []).join('');
    const err = (res.logs?.stderr || []).join('');
    return { 
      ok: res.exitCode === 0, 
      exitCode: res.exitCode, 
      stdout: out, 
      stderr: err 
    };
  } catch (e) {
    return { 
      ok: false, 
      exitCode: -1, 
      stdout: '', 
      stderr: String(e?.message || e) 
    };
  } finally {
    clearTimeout(overallTimer);
  }
}

// Cleanup on process exit
process.on('SIGINT', async () => { 
  try { 
    if (sbx && !sbx.closed) await sbx.kill(); 
  } catch {} 
  process.exit(0); 
});

process.on('SIGTERM', async () => { 
  try { 
    if (sbx && !sbx.closed) await sbx.kill(); 
  } catch {} 
  process.exit(0); 
});

module.exports = {
  ensureSandbox,
  runPythonSafe,
  getStatus: () => ({
    healthy: !!sbx && !sbx.closed,
    lastTouch,
    templateId: TEMPLATE
  })
};