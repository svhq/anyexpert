// e2b-manager-v2.js - Production-ready E2B sandbox manager with proper lifecycle management
const { Sandbox } = require('@e2b/code-interpreter');
const logger = require('./utils/logger');

class E2BManagerV2 {
  constructor() {
    this.templateId = process.env.E2B_TEMPLATE_ID || 'prod-all';
    
    // Critical: Production configuration for concurrent users
    this.config = {
      maxTotalSandboxes: 20,          // Hard limit to prevent hitting E2B's 100 limit
      maxWarmSandboxes: 3,            // Keep 3 warm for burst traffic
      maxSandboxAge: 5 * 60 * 1000,   // Kill sandboxes after 5 minutes (not 30)
      maxReuseCount: 10,              // Reuse sandbox max 10 times before refresh
      healthCheckInterval: 60 * 1000,  // Check health every minute
      cleanupInterval: 30 * 1000      // Cleanup dead sandboxes every 30s
    };
    
    // Sandbox tracking
    this.activeSandboxes = new Map(); // All sandboxes we've created
    this.warmPool = [];               // Available for immediate use
    this.inUse = new Set();          // Currently executing code
    
    // Metrics
    this.metrics = {
      created: 0,
      reused: 0,
      killed: 0,
      failed: 0,
      concurrent: 0,
      peakConcurrent: 0
    };
    
    // Start background tasks
    this.startBackgroundTasks();
  }

  startBackgroundTasks() {
    // Periodic cleanup of old/dead sandboxes
    this.cleanupTimer = setInterval(() => {
      this.cleanupSandboxes().catch(err => 
        logger.error({ message: 'Cleanup failed', error: err.message })
      );
    }, this.config.cleanupInterval);
    
    // Periodic health checks
    this.healthTimer = setInterval(() => {
      this.healthCheckPool().catch(err => 
        logger.error({ message: 'Health check failed', error: err.message })
      );
    }, this.config.healthCheckInterval);
  }

  async cleanupSandboxes() {
    const now = Date.now();
    const toKill = [];
    
    // Check all active sandboxes
    for (const [id, entry] of this.activeSandboxes) {
      const age = now - entry.createdAt;
      
      // Kill if: too old, too many reuses, or unhealthy
      if (age > this.config.maxSandboxAge || 
          entry.reuseCount > this.config.maxReuseCount ||
          !entry.healthy) {
        toKill.push(id);
      }
    }
    
    // Kill identified sandboxes
    for (const id of toKill) {
      await this.killSandbox(id);
    }
    
    // Ensure we don't exceed max total sandboxes
    if (this.activeSandboxes.size > this.config.maxTotalSandboxes) {
      const excess = this.activeSandboxes.size - this.config.maxTotalSandboxes;
      const sorted = Array.from(this.activeSandboxes.entries())
        .sort((a, b) => a[1].createdAt - b[1].createdAt); // Oldest first
      
      for (let i = 0; i < excess; i++) {
        await this.killSandbox(sorted[i][0]);
      }
    }
    
    logger.info({ 
      message: 'Cleanup complete',
      killed: toKill.length,
      active: this.activeSandboxes.size,
      pool: this.warmPool.length
    });
  }

  async healthCheckPool() {
    const unhealthy = [];
    
    for (const entry of this.warmPool) {
      try {
        const result = await entry.sandbox.runCode('1', { timeoutMs: 2000 });
        if (result.error || result.results?.[0]?.text !== '1') {
          unhealthy.push(entry);
        }
      } catch (e) {
        unhealthy.push(entry);
      }
    }
    
    // Remove unhealthy sandboxes from pool
    for (const entry of unhealthy) {
      const idx = this.warmPool.indexOf(entry);
      if (idx !== -1) {
        this.warmPool.splice(idx, 1);
        await this.killSandbox(entry.sandbox.id);
      }
    }
  }

  async killSandbox(sandboxId) {
    const entry = this.activeSandboxes.get(sandboxId);
    if (!entry) return;
    
    try {
      await entry.sandbox.kill();
    } catch (e) {
      // Ignore kill errors
    }
    
    this.activeSandboxes.delete(sandboxId);
    this.inUse.delete(sandboxId);
    this.metrics.killed++;
    
    // Remove from warm pool if present
    const poolIdx = this.warmPool.findIndex(e => e.sandbox.id === sandboxId);
    if (poolIdx !== -1) {
      this.warmPool.splice(poolIdx, 1);
    }
  }

  async getSandbox() {
    // Update concurrent usage metrics
    this.metrics.concurrent = this.inUse.size;
    if (this.metrics.concurrent > this.metrics.peakConcurrent) {
      this.metrics.peakConcurrent = this.metrics.concurrent;
    }
    
    // Try to get from warm pool first
    while (this.warmPool.length > 0) {
      const entry = this.warmPool.shift();
      
      // Quick health check
      try {
        const result = await entry.sandbox.runCode('1', { timeoutMs: 1000 });
        if (!result.error && result.results?.[0]?.text === '1') {
          this.inUse.add(entry.sandbox.id);
          this.metrics.reused++;
          logger.debug({ message: 'Reusing sandbox', id: entry.sandbox.id });
          return entry;
        }
      } catch (e) {
        // Sandbox is dead
      }
      
      // Remove dead sandbox
      await this.killSandbox(entry.sandbox.id);
    }
    
    // Check if we can create a new sandbox
    if (this.activeSandboxes.size >= this.config.maxTotalSandboxes) {
      // Try to kill oldest unused sandbox
      const unused = Array.from(this.activeSandboxes.entries())
        .filter(([id]) => !this.inUse.has(id))
        .sort((a, b) => a[1].createdAt - b[1].createdAt);
      
      if (unused.length > 0) {
        await this.killSandbox(unused[0][0]);
      } else {
        throw new Error(`Sandbox limit reached (${this.config.maxTotalSandboxes}). Too many concurrent executions.`);
      }
    }
    
    // Create new sandbox
    try {
      const sandbox = await Sandbox.create({
        template: this.templateId,
        timeoutMs: this.config.maxSandboxAge
      });
      
      const entry = {
        sandbox,
        createdAt: Date.now(),
        reuseCount: 0,
        healthy: true
      };
      
      this.activeSandboxes.set(sandbox.id, entry);
      this.inUse.add(sandbox.id);
      this.metrics.created++;
      
      logger.info({ 
        message: 'Created new sandbox',
        id: sandbox.id,
        active: this.activeSandboxes.size
      });
      
      return entry;
    } catch (error) {
      this.metrics.failed++;
      throw error;
    }
  }

  async returnSandbox(entry) {
    const sandboxId = entry.sandbox.id;
    this.inUse.delete(sandboxId);
    
    // Don't return if we have enough in warm pool
    if (this.warmPool.length >= this.config.maxWarmSandboxes) {
      await this.killSandbox(sandboxId);
      return;
    }
    
    // Don't return if sandbox is old or overused
    const age = Date.now() - entry.createdAt;
    if (age > this.config.maxSandboxAge || entry.reuseCount >= this.config.maxReuseCount) {
      await this.killSandbox(sandboxId);
      return;
    }
    
    // Health check before returning to pool
    try {
      const result = await entry.sandbox.runCode('1', { timeoutMs: 1000 });
      if (!result.error && result.results?.[0]?.text === '1') {
        entry.reuseCount++;
        this.warmPool.push(entry);
        logger.debug({ 
          message: 'Returned to pool',
          id: sandboxId,
          reuseCount: entry.reuseCount
        });
      } else {
        await this.killSandbox(sandboxId);
      }
    } catch (e) {
      await this.killSandbox(sandboxId);
    }
  }

  async executeCode(code, options = {}) {
    const { language = 'python', timeoutMs = 30000, requestId = null } = options;
    const startTime = Date.now();
    
    let sandboxEntry;
    try {
      sandboxEntry = await this.getSandbox();
      
      // Wrap Python code with timeout protection
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
      
      const result = await sandboxEntry.sandbox.runCode(wrappedCode, {
        language,
        timeoutMs: timeoutMs + 5000
      });
      
      const isSuccess = !result.error;
      
      logger.info({ 
        requestId,
        executionTime: Date.now() - startTime,
        success: isSuccess,
        message: 'Code execution completed'
      });
      
      // Return sandbox to pool if successful
      if (isSuccess) {
        await this.returnSandbox(sandboxEntry);
      } else {
        await this.killSandbox(sandboxEntry.sandbox.id);
      }
      
      return {
        success: isSuccess,
        stdout: result.logs?.stdout?.join('\n') || '',
        stderr: result.logs?.stderr?.join('\n') || result.error || '',
        exitCode: result.error ? 1 : 0,
        executionTime: Date.now() - startTime
      };
      
    } catch (error) {
      if (sandboxEntry) {
        await this.killSandbox(sandboxEntry.sandbox.id);
      }
      
      return {
        success: false,
        stdout: '',
        stderr: `E2B execution failed: ${error.message}`,
        exitCode: 1,
        executionTime: Date.now() - startTime
      };
    }
  }

  async getStatus() {
    return {
      active: this.activeSandboxes.size,
      inUse: this.inUse.size,
      warmPool: this.warmPool.length,
      metrics: this.metrics,
      canCreate: this.activeSandboxes.size < this.config.maxTotalSandboxes
    };
  }

  async shutdown() {
    // Clear timers
    if (this.cleanupTimer) clearInterval(this.cleanupTimer);
    if (this.healthTimer) clearInterval(this.healthTimer);
    
    // Kill all sandboxes
    const killPromises = [];
    for (const [id] of this.activeSandboxes) {
      killPromises.push(this.killSandbox(id));
    }
    
    await Promise.all(killPromises);
    
    logger.info({ 
      message: 'E2B manager shutdown complete',
      metrics: this.metrics
    });
  }
}

// Export singleton instance
module.exports = new E2BManagerV2();

// Graceful shutdown
process.on('SIGINT', async () => {
  await module.exports.shutdown();
  process.exit(0);
});

process.on('SIGTERM', async () => {
  await module.exports.shutdown();
  process.exit(0);
});