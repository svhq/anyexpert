// npf-smoke.mjs
import { Sandbox } from '@e2b/code-interpreter';

const sbx = await Sandbox.create({
  template: process.env.E2B_TEMPLATE_ID || 'prod-all',
  apiKey: process.env.E2B_API_KEY,
  timeout: 60_000
});

console.log(`Created sandbox ${sbx.sandboxId} with template ${process.env.E2B_TEMPLATE_ID || 'prod-all'}`);

const code = `
import sys, pkgutil, os
import numpy as np
print("numpy:", np.__version__)
print("has numpy_financial:", pkgutil.find_loader("numpy_financial") is not None)
print("BUILD_TAG:", os.environ.get("BUILD_TAG", "missing"))
try:
    import numpy_financial as npf
    print("pv test:", round(npf.pv(0.05/12, 36, 300), 2))
except ImportError as e:
    print("IMPORT ERROR:", e)
`;

const r = await sbx.runCode(code, { language: 'python', timeout: 15_000 });
console.log(r.logs.stdout.join(''));
console.error(r.logs.stderr.join(''));
await sbx.kill();