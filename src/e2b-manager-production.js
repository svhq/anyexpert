// e2b-manager-production.js - Production-ready with safe cleanup
const { Sandbox } = require('@e2b/code-interpreter');
const logger = require('./utils/logger');

class E2BManagerProduction {
  constructor() {
    this.templateId = process.env.E2B_TEMPLATE_ID || 'prod-all';
    
    // Sandbox tracking - THREE states
    this.sandboxes = {
      inUse: new Map(),      // Currently executing code - NEVER kill these
      warmPool: [],          // Available for reuse
      idle: new Map()        // Returned but not in warm pool
    };
    
    // More lenient configuration for production
    this.config = {
      maxWarmPool: 5,              // Keep 5 ready for quick reuse
      idleTimeout: 2 * 60 * 1000,  // Kill idle sandboxes after 2 minutes
      maxAge: 10 * 60 * 1000,      // Kill any sandbox after 10 minutes
      cleanupInterval: 30 * 1000,   // Run cleanup every 30 seconds
      maxReuses: 25                 // Reuse sandbox up to 25 times
    };
    
    // Metrics
    this.metrics = {
      created: 0,
      reused: 0,
      killed: 0,
      peakConcurrent: 0,
      activeConcurrent: 0
    };
    
    // Start cleanup timer
    this.startCleanup();
  }

  startCleanup() {
    this.cleanupTimer = setInterval(() => {
      this.cleanup().catch(err => 
        logger.error({ message: 'Cleanup error', error: err.message })
      );
    }, this.config.cleanupInterval);
  }

  async cleanup() {
    const now = Date.now();
    const toKill = [];
    
    logger.debug({ 
      message: 'Cleanup running',
      inUse: this.sandboxes.inUse.size,
      warmPool: this.sandboxes.warmPool.length,
      idle: this.sandboxes.idle.size
    });
    
    // 1. NEVER touch sandboxes that are in use
    // These are actively running user code
    
    // 2. Clean up idle sandboxes (not in use, not in warm pool)
    for (const [id, entry] of this.sandboxes.idle) {
      const idleTime = now - entry.lastUsed;
      const age = now - entry.createdAt;
      
      // Kill if idle too long OR too old OR unhealthy
      if (idleTime > this.config.idleTimeout || 
          age > this.config.maxAge ||
          entry.reuseCount > this.config.maxReuses) {
        toKill.push({ id, reason: 'idle-timeout' });
      }
    }
    
    // 3. Check warm pool health (but don't kill if healthy)
    const unhealthyPool = [];
    for (const entry of this.sandboxes.warmPool) {
      const age = now - entry.createdAt;
      if (age > this.config.maxAge || entry.reuseCount > this.config.maxReuses) {
        unhealthyPool.push(entry);
      }
    }
    
    // Remove unhealthy from warm pool
    for (const entry of unhealthyPool) {
      const idx = this.sandboxes.warmPool.indexOf(entry);
      if (idx !== -1) {
        this.sandboxes.warmPool.splice(idx, 1);
        toKill.push({ id: entry.sandbox.id, reason: 'pool-unhealthy' });
      }
    }
    
    // Kill identified sandboxes
    let killed = 0;
    for (const { id, reason } of toKill) {
      try {
        const entry = this.sandboxes.idle.get(id) || 
                     this.sandboxes.warmPool.find(e => e.sandbox.id === id);
        if (entry) {
          await entry.sandbox.kill();
          this.sandboxes.idle.delete(id);
          killed++;
          logger.info({ 
            message: 'Killed sandbox',
            id,
            reason,
            reuseCount: entry.reuseCount,
            age: Math.floor((now - entry.createdAt) / 1000) + 's'
          });
        }
      } catch (e) {
        // Ignore kill errors
      }
    }
    
    this.metrics.killed += killed;
    
    if (killed > 0) {
      logger.info({ 
        message: 'Cleanup complete',
        killed,
        remaining: {
          inUse: this.sandboxes.inUse.size,
          warmPool: this.sandboxes.warmPool.length,
          idle: this.sandboxes.idle.size
        }
      });
    }
  }

  async getSandbox() {
    // Update concurrent metrics
    this.metrics.activeConcurrent = this.sandboxes.inUse.size;
    if (this.metrics.activeConcurrent > this.metrics.peakConcurrent) {
      this.metrics.peakConcurrent = this.metrics.activeConcurrent;
    }
    
    // Try warm pool first
    while (this.sandboxes.warmPool.length > 0) {
      const entry = this.sandboxes.warmPool.shift();
      
      try {
        // Quick health check
        const result = await entry.sandbox.runCode('1', { timeoutMs: 1000 });
        if (!result.error && result.results?.[0]?.text === '1') {
          // Move to inUse
          this.sandboxes.inUse.set(entry.sandbox.id, entry);
          this.metrics.reused++;
          logger.debug({ 
            message: 'Reusing sandbox from pool',
            id: entry.sandbox.id,
            reuseCount: entry.reuseCount
          });
          return entry;
        }
      } catch (e) {
        // Dead sandbox
      }
      
      // Kill dead sandbox
      try {
        await entry.sandbox.kill();
      } catch (e) {}
    }
    
    // Create new sandbox
    try {
      const sandbox = await Sandbox.create({
        template: this.templateId,
        timeoutMs: this.config.maxAge
      });
      
      const entry = {
        sandbox,
        createdAt: Date.now(),
        lastUsed: Date.now(),
        reuseCount: 0
      };
      
      this.sandboxes.inUse.set(sandbox.id, entry);
      this.metrics.created++;
      
      logger.info({ 
        message: 'Created new sandbox',
        id: sandbox.id,
        totalActive: this.getTotalSandboxes()
      });
      
      return entry;
    } catch (error) {
      logger.error({ 
        message: 'Failed to create sandbox',
        error: error.message,
        totalActive: this.getTotalSandboxes()
      });
      throw error;
    }
  }

  async returnSandbox(entry, wasSuccessful) {
    const id = entry.sandbox.id;
    
    // Remove from inUse
    this.sandboxes.inUse.delete(id);
    entry.lastUsed = Date.now();
    
    // Don't return if execution failed
    if (!wasSuccessful) {
      try {
        await entry.sandbox.kill();
        this.metrics.killed++;
      } catch (e) {}
      return;
    }
    
    // Don't return if too old or overused
    const age = Date.now() - entry.createdAt;
    if (age > this.config.maxAge || entry.reuseCount >= this.config.maxReuses) {
      try {
        await entry.sandbox.kill();
        this.metrics.killed++;
      } catch (e) {}
      return;
    }
    
    // If warm pool has space, add there
    if (this.sandboxes.warmPool.length < this.config.maxWarmPool) {
      entry.reuseCount++;
      this.sandboxes.warmPool.push(entry);
      logger.debug({ 
        message: 'Returned to warm pool',
        id,
        reuseCount: entry.reuseCount,
        poolSize: this.sandboxes.warmPool.length
      });
    } else {
      // Otherwise, mark as idle
      this.sandboxes.idle.set(id, entry);
      logger.debug({ 
        message: 'Marked as idle',
        id,
        idleCount: this.sandboxes.idle.size
      });
    }
  }

  async executeCode(code, options = {}) {
    const { language = 'python', timeoutMs = 30000, requestId = null } = options;
    const startTime = Date.now();
    
    let sandboxEntry;
    try {
      sandboxEntry = await this.getSandbox();
      
      // Add timeout protection
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
    signal.alarm(0)
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
      
      // Return sandbox for reuse
      await this.returnSandbox(sandboxEntry, isSuccess);
      
      return {
        success: isSuccess,
        stdout: result.logs?.stdout?.join('\n') || '',
        stderr: result.logs?.stderr?.join('\n') || result.error || '',
        exitCode: result.error ? 1 : 0,
        executionTime: Date.now() - startTime
      };
      
    } catch (error) {
      // Make sure sandbox is removed from inUse if error occurs
      if (sandboxEntry) {
        this.sandboxes.inUse.delete(sandboxEntry.sandbox.id);
        try {
          await sandboxEntry.sandbox.kill();
        } catch (e) {}
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

  getTotalSandboxes() {
    return this.sandboxes.inUse.size + 
           this.sandboxes.warmPool.length + 
           this.sandboxes.idle.size;
  }

  getStatus() {
    return {
      total: this.getTotalSandboxes(),
      inUse: this.sandboxes.inUse.size,
      warmPool: this.sandboxes.warmPool.length,
      idle: this.sandboxes.idle.size,
      metrics: this.metrics
    };
  }

  async shutdown() {
    if (this.cleanupTimer) {
      clearInterval(this.cleanupTimer);
    }
    
    // Kill all sandboxes
    const allSandboxes = [
      ...Array.from(this.sandboxes.inUse.values()),
      ...this.sandboxes.warmPool,
      ...Array.from(this.sandboxes.idle.values())
    ];
    
    const killPromises = allSandboxes.map(async entry => {
      try {
        await entry.sandbox.kill();
      } catch (e) {}
    });
    
    await Promise.all(killPromises);
    
    logger.info({ 
      message: 'E2B manager shutdown complete',
      metrics: this.metrics
    });
  }
}

module.exports = E2BManagerProduction;