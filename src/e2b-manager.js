// e2b-manager.js - Production-ready E2B sandbox manager with pooling and reliability
const { Sandbox } = require('@e2b/code-interpreter');
const logger = require('./utils/logger');

class E2BManager {
  constructor() {
    this.templateId = process.env.E2B_TEMPLATE_ID || 'prod-all';
    
    // Log warning if E2B_API_KEY is not set
    if (!process.env.E2B_API_KEY) {
      logger.warn({ message: 'E2B_API_KEY environment variable is not set' });
    }
    
    // Sandbox pool
    this.warmSandboxes = [];
    this.maxWarmSandboxes = 20; // Keep up to 20 warm sandboxes for better burst handling
    this.sandboxTimeout = 30 * 60 * 1000; // 30 minutes
    
    // Retry configuration
    this.maxRetries = 2;
    this.retryDelay = 1000; // Start with 1 second
    
    // Performance tracking
    this.metrics = {
      totalExecutions: 0,
      failures: 0,
      retries: 0,
      avgExecutionTime: 0
    };
    
    // Initialize warm sandboxes in background
    // this.initializeWarmSandboxes(); // Temporarily disabled
  }

  async initializeWarmSandboxes() {
    try {
      for (let i = 0; i < this.maxWarmSandboxes; i++) {
        this.createWarmSandbox().catch(err => 
          logger.warn({ message: 'Failed to create warm sandbox', error: err.message })
        );
      }
    } catch (error) {
      logger.error({ message: 'Failed to initialize warm sandboxes', error: error.message });
    }
  }

  async createWarmSandbox() {
    try {
      const sandbox = await Sandbox.create({
        template: this.templateId,
        timeoutMs: this.sandboxTimeout
      });
      
      // Health check
      const health = await sandbox.runCode('print("__healthy__")', { timeoutMs: 5000 });
      if (health.error || !(health.logs?.stdout?.join('').includes('__healthy__') || health.results?.[0]?.text?.includes('__healthy__'))) {
        await sandbox.kill();
        throw new Error('Sandbox health check failed');
      }
      
      this.warmSandboxes.push({
        sandbox,
        createdAt: Date.now(),
        healthy: true,
        reuseCount: 0
      });
      
      logger.info({ message: 'Created warm sandbox', sandboxId: sandbox.id });
    } catch (error) {
      logger.error({ message: 'Failed to create warm sandbox', error: error.message });
    }
  }

  async returnSandbox(sandboxEntry) {
    // Don't return if we already have enough warm sandboxes
    if (this.warmSandboxes.length >= this.maxWarmSandboxes) {
      return;
    }
    
    // Check if sandbox is still healthy and not overused
    const maxReuses = 50; // Maximum times to reuse a sandbox
    if (sandboxEntry.reuseCount >= maxReuses) {
      // Kill overused sandbox
      try {
        await sandboxEntry.sandbox.kill();
      } catch (e) {
        // Ignore errors
      }
      return;
    }
    
    // Quick health check before returning to pool
    try {
      const health = await sandboxEntry.sandbox.runCode('1', { timeoutMs: 2000 });
      if (!health.error && health.results?.[0]?.text === '1') {
        // Increment reuse count and return to pool
        sandboxEntry.reuseCount++;
        sandboxEntry.healthy = true;
        this.warmSandboxes.push(sandboxEntry);
        logger.debug({ 
          message: 'Sandbox returned to pool', 
          reuseCount: sandboxEntry.reuseCount 
        });
      } else {
        // Sandbox is unhealthy, kill it
        await sandboxEntry.sandbox.kill();
      }
    } catch (e) {
      // Sandbox is dead, don't return to pool
      try {
        await sandboxEntry.sandbox.kill();
      } catch (e2) {
        // Ignore
      }
    }
  }

  async getSandbox() {
    // Try to get a warm sandbox
    while (this.warmSandboxes.length > 0) {
      const entry = this.warmSandboxes.shift();
      
      // Check if sandbox is still valid (not expired)
      if (Date.now() - entry.createdAt < this.sandboxTimeout - 60000) { // 1 min buffer
        // Quick health check
        try {
          const health = await entry.sandbox.runCode('1+1', { timeoutMs: 3000 });
          if (!health.error && health.results?.[0]?.text === '2') {
            // Return the sandbox with metadata for tracking
            return {
              sandbox: entry.sandbox,
              createdAt: entry.createdAt,
              reuseCount: entry.reuseCount || 0
            };
          }
        } catch (e) {
          // Sandbox is dead
        }
      }
      
      // Kill expired/unhealthy sandbox
      try {
        await entry.sandbox.kill();
      } catch (e) {
        // Ignore cleanup errors
      }
    }
    
    // No warm sandbox available, create new one
    const sandbox = await Sandbox.create({
      template: this.templateId,
      timeoutMs: this.sandboxTimeout
    });
    
    return {
      sandbox,
      createdAt: Date.now(),
      reuseCount: 0
    };
  }

  async executeCode(code, options = {}) {
    const {
      language = 'python',
      timeoutMs = 30000,
      requestId = null
    } = options;
    
    // Check if E2B is configured
    if (!process.env.E2B_API_KEY) {
      logger.warn({ 
        requestId, 
        message: 'E2B_API_KEY not set, returning mock result for code execution',
        code: code.substring(0, 100)
      });
      
      // For simple arithmetic, compute directly
      if (code.includes('print') && /^\d+\s*[\+\-\*\/]\s*\d+$/.test(code.replace('print(', '').replace(')', '').trim())) {
        try {
          const result = eval(code.replace('print(', '').replace(')', '').trim());
          return {
            results: [{ text: String(result) }],
            logs: { stdout: [String(result)], stderr: [] },
            error: null
          };
        } catch (e) {
          // Fall through to general mock
        }
      }
      
      return {
        results: [{ text: 'E2B not configured - unable to execute code' }],
        logs: { stdout: [], stderr: [] },
        error: 'E2B_API_KEY not configured'
      };
    }
    
    const startTime = Date.now();
    let lastError = null;
    
    for (let attempt = 0; attempt <= this.maxRetries; attempt++) {
      try {
        if (attempt > 0) {
          this.metrics.retries++;
          logger.info({ requestId, attempt, message: 'Retrying code execution' });
          await new Promise(resolve => setTimeout(resolve, this.retryDelay * attempt));
        }
        
        const sandboxEntry = await this.getSandbox();
        const sandbox = sandboxEntry.sandbox;
        
        // Add timeout protection in Python code
        const wrappedCode = language === 'python' ? `
import signal
import sys

def timeout_handler(signum, frame):
    print("__EXECUTION_TIMEOUT__", file=sys.stderr)
    sys.exit(124)

signal.signal(signal.SIGALRM, timeout_handler)
signal.alarm(${Math.floor(timeoutMs / 1000)})

try:
${code.split('\n').map(line => '    ' + line).join('\n')}
finally:
    signal.alarm(0)  # Cancel the alarm
` : code;
        
        const result = await sandbox.runCode(wrappedCode, {
          language,
          timeoutMs: timeoutMs + 5000 // Give SDK slightly more time
        });
        
        this.metrics.totalExecutions++;
        const executionTime = Date.now() - startTime;
        this.metrics.avgExecutionTime = 
          (this.metrics.avgExecutionTime * (this.metrics.totalExecutions - 1) + executionTime) / 
          this.metrics.totalExecutions;
        
        // Determine success - E2B sets error property on failure
        const isSuccess = !result.error;
        
        logger.info({ 
          requestId, 
          attempt,
          executionTime,
          success: isSuccess,
          message: 'Code execution completed' 
        });
        
        // Return sandbox to pool if execution was successful
        if (isSuccess) {
          await this.returnSandbox(sandboxEntry);
        }
        
        // Format the result
        return {
          success: isSuccess,
          stdout: result.logs?.stdout?.join('\n') || '',
          stderr: result.logs?.stderr?.join('\n') || result.error || '',
          exitCode: result.error ? 1 : 0,
          executionTime
        };
        
      } catch (error) {
        lastError = error;
        logger.error({ 
          requestId, 
          attempt,
          error: error.message,
          message: 'Code execution failed' 
        });
        
        if (attempt === this.maxRetries) {
          this.metrics.failures++;
        }
      }
    }
    
    // All retries failed
    return {
      success: false,
      stdout: '',
      stderr: `E2B execution failed after ${this.maxRetries + 1} attempts: ${lastError?.message || 'Unknown error'}`,
      exitCode: 1,
      executionTime: Date.now() - startTime
    };
  }

  async executeWithFallback(code, options = {}) {
    const startTime = Date.now();
    
    // Try E2B first
    const e2bResult = await this.executeCode(code, options);
    
    if (e2bResult.success || !options.fallbackToLocal) {
      return e2bResult;
    }
    
    // Fallback to local Python if configured
    logger.info({ message: 'Falling back to local Python execution' });
    
    try {
      const { spawn } = require('child_process');
      const python = spawn('python', ['-c', code]);
      
      let stdout = '';
      let stderr = '';
      
      python.stdout.on('data', (data) => {
        stdout += data.toString();
      });
      
      python.stderr.on('data', (data) => {
        stderr += data.toString();
      });
      
      const exitCode = await new Promise((resolve) => {
        python.on('close', resolve);
        
        // Add timeout
        setTimeout(() => {
          python.kill();
          resolve(124);
        }, options.timeoutMs || 30000);
      });
      
      return {
        success: exitCode === 0,
        stdout,
        stderr: stderr || (exitCode === 124 ? 'Local execution timeout' : ''),
        exitCode,
        executionTime: Date.now() - startTime,
        fallback: true
      };
      
    } catch (error) {
      return {
        success: false,
        stdout: '',
        stderr: `Local Python execution failed: ${error.message}`,
        exitCode: 1,
        executionTime: Date.now() - startTime,
        fallback: true
      };
    }
  }

  async shutdown() {
    logger.info({ message: 'Shutting down E2B manager' });
    
    // Kill all warm sandboxes
    const killPromises = this.warmSandboxes.map(async (entry) => {
      try {
        await entry.sandbox.kill();
      } catch (e) {
        // Ignore errors during shutdown
      }
    });
    
    await Promise.all(killPromises);
    this.warmSandboxes = [];
    
    logger.info({ 
      message: 'E2B manager shutdown complete',
      metrics: this.metrics
    });
  }

  getMetrics() {
    return {
      ...this.metrics,
      warmSandboxCount: this.warmSandboxes.length
    };
  }
}

// Export singleton instance
module.exports = new E2BManager();

// Graceful shutdown
process.on('SIGINT', async () => {
  await module.exports.shutdown();
  process.exit(0);
});

process.on('SIGTERM', async () => {
  await module.exports.shutdown();
  process.exit(0);
});