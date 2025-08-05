// e2b-manager-improved-v2.js - Production-ready E2B sandbox manager with serialization and pre-warming
const { Sandbox } = require('@e2b/code-interpreter');
const PQueue = require('p-queue').default;
const logger = require('./utils/logger');

class E2BManagerImproved {
  constructor() {
    this.templateId = process.env.E2B_TEMPLATE_ID || 'prod-all';
    
    // Log warning if E2B_API_KEY is not set
    if (!process.env.E2B_API_KEY) {
      logger.warn({ message: 'E2B_API_KEY environment variable is not set' });
    }
    
    // Sandbox pool configuration
    this.warmSandboxes = [];
    this.targetWarmSandboxes = 3; // Keep 2-3 warm sandboxes
    this.maxWarmSandboxes = 10; // Allow growth under load
    this.sandboxTimeout = 30 * 60 * 1000; // 30 minutes
    
    // Serialization: Only 1 sandbox creation at a time
    this.createQueue = new PQueue({ concurrency: 1 });
    
    // Circuit breaker configuration
    this.circuitBreaker = {
      failures: 0,
      lastFailureTime: 0,
      threshold: 5, // Open after 5 failures
      timeout: 30000, // 30 second cooldown
      state: 'closed' // closed, open, half-open
    };
    
    // Retry configuration with exponential backoff
    this.maxRetries = 3;
    this.baseRetryDelay = 250; // Start with 250ms
    this.maxRetryDelay = 2000; // Cap at 2 seconds
    
    // Performance tracking
    this.metrics = {
      totalExecutions: 0,
      successfulExecutions: 0,
      failures: 0,
      retries: 0,
      creationAttempts: 0,
      creationSuccesses: 0,
      creationFailures: 0,
      avgExecutionTime: 0
    };
    
    // Initialize warm sandboxes in background with proper spacing
    this.initializeWarmSandboxes();
    
    // Start monitoring
    this.startMonitoring();
  }

  async initializeWarmSandboxes() {
    logger.info({ 
      message: 'Starting sandbox pre-warming',
      target: this.targetWarmSandboxes 
    });
    
    // Create sandboxes one at a time with delays
    for (let i = 0; i < this.targetWarmSandboxes; i++) {
      // Add jitter to prevent thundering herd
      const delay = 500 + Math.random() * 500; // 500-1000ms
      await new Promise(resolve => setTimeout(resolve, delay));
      
      this.createWarmSandbox().catch(err => {
        logger.warn({ 
          message: 'Failed to create warm sandbox during initialization',
          error: err.message,
          attempt: i + 1
        });
      });
    }
  }

  async createWarmSandbox() {
    // Check circuit breaker
    if (this.isCircuitOpen()) {
      throw new Error('Circuit breaker is open - sandbox creation temporarily disabled');
    }
    
    // Queue the creation to ensure serialization
    return this.createQueue.add(async () => {
      const startTime = Date.now();
      this.metrics.creationAttempts++;
      
      try {
        // Exponential backoff for retries
        let lastError;
        for (let attempt = 1; attempt <= this.maxRetries; attempt++) {
          try {
            if (attempt > 1) {
              const backoffMs = Math.min(
                this.baseRetryDelay * Math.pow(2, attempt - 1),
                this.maxRetryDelay
              );
              const jitter = Math.random() * 100;
              await new Promise(resolve => setTimeout(resolve, backoffMs + jitter));
              logger.info({ 
                message: 'Retrying sandbox creation',
                attempt,
                backoffMs: backoffMs + jitter
              });
            }
            
            const sandbox = await Sandbox.create({
              template: this.templateId,
              timeoutMs: this.sandboxTimeout
            });
            
            // Health check
            const health = await sandbox.runCode('print("__healthy__")', { timeoutMs: 5000 });
            if (health.error || !health.logs?.stdout?.join('').includes('__healthy__')) {
              await sandbox.kill();
              throw new Error('Sandbox health check failed');
            }
            
            const sandboxEntry = {
              sandbox,
              createdAt: Date.now(),
              healthy: true,
              reuseCount: 0,
              inUse: false
            };
            
            this.warmSandboxes.push(sandboxEntry);
            this.metrics.creationSuccesses++;
            
            // Reset circuit breaker on success
            this.circuitBreaker.failures = 0;
            
            logger.info({ 
              message: 'Created warm sandbox',
              sandboxId: sandbox.id,
              duration: Date.now() - startTime,
              warmCount: this.warmSandboxes.length
            });
            
            return sandboxEntry;
            
          } catch (error) {
            lastError = error;
            if (attempt === this.maxRetries) {
              throw error;
            }
          }
        }
        throw lastError;
        
      } catch (error) {
        this.metrics.creationFailures++;
        this.handleCreationFailure(error);
        throw error;
      }
    });
  }

  handleCreationFailure(error) {
    this.circuitBreaker.failures++;
    this.circuitBreaker.lastFailureTime = Date.now();
    
    if (this.circuitBreaker.failures >= this.circuitBreaker.threshold) {
      this.circuitBreaker.state = 'open';
      logger.error({ 
        message: 'Circuit breaker opened due to repeated failures',
        failures: this.circuitBreaker.failures,
        error: error.message
      });
    }
  }

  isCircuitOpen() {
    if (this.circuitBreaker.state === 'closed') {
      return false;
    }
    
    // Check if cooldown period has passed
    const timeSinceLastFailure = Date.now() - this.circuitBreaker.lastFailureTime;
    if (timeSinceLastFailure > this.circuitBreaker.timeout) {
      this.circuitBreaker.state = 'half-open';
      this.circuitBreaker.failures = 0;
      return false;
    }
    
    return true;
  }

  async returnSandbox(sandboxEntry) {
    sandboxEntry.inUse = false;
    
    // Don't keep too many warm sandboxes
    if (this.warmSandboxes.length > this.maxWarmSandboxes) {
      const index = this.warmSandboxes.indexOf(sandboxEntry);
      if (index > -1) {
        this.warmSandboxes.splice(index, 1);
        try {
          await sandboxEntry.sandbox.kill();
        } catch (e) {
          // Ignore errors
        }
      }
      return;
    }
    
    // Check if sandbox is still healthy and not overused
    const maxReuses = 50;
    if (sandboxEntry.reuseCount >= maxReuses) {
      const index = this.warmSandboxes.indexOf(sandboxEntry);
      if (index > -1) {
        this.warmSandboxes.splice(index, 1);
      }
      try {
        await sandboxEntry.sandbox.kill();
      } catch (e) {
        // Ignore errors
      }
      return;
    }
    
    // Quick health check before keeping in pool
    try {
      const health = await sandboxEntry.sandbox.runCode('1', { timeoutMs: 2000 });
      if (!health.error && health.results?.[0]?.text === '1') {
        sandboxEntry.reuseCount++;
        sandboxEntry.healthy = true;
        logger.debug({ 
          message: 'Sandbox returned to pool',
          reuseCount: sandboxEntry.reuseCount,
          warmCount: this.warmSandboxes.length
        });
      } else {
        // Remove unhealthy sandbox
        const index = this.warmSandboxes.indexOf(sandboxEntry);
        if (index > -1) {
          this.warmSandboxes.splice(index, 1);
        }
        await sandboxEntry.sandbox.kill();
      }
    } catch (e) {
      // Remove dead sandbox
      const index = this.warmSandboxes.indexOf(sandboxEntry);
      if (index > -1) {
        this.warmSandboxes.splice(index, 1);
      }
      try {
        await sandboxEntry.sandbox.kill();
      } catch (e2) {
        // Ignore
      }
    }
    
    // Maintain target warm sandbox count
    this.maintainWarmPool();
  }

  async maintainWarmPool() {
    const availableCount = this.warmSandboxes.filter(s => !s.inUse).length;
    
    if (availableCount < this.targetWarmSandboxes) {
      // Create one more, but don't block
      this.createWarmSandbox().catch(err => {
        logger.debug({ 
          message: 'Failed to maintain warm pool',
          error: err.message
        });
      });
    }
  }

  async getSandbox() {
    // Try to get an available warm sandbox
    const available = this.warmSandboxes.find(s => !s.inUse && s.healthy);
    
    if (available) {
      // Check if sandbox is still valid (not expired)
      if (Date.now() - available.createdAt < this.sandboxTimeout - 60000) { // 1 min buffer
        available.inUse = true;
        return available;
      } else {
        // Remove expired sandbox
        const index = this.warmSandboxes.indexOf(available);
        if (index > -1) {
          this.warmSandboxes.splice(index, 1);
        }
        try {
          await available.sandbox.kill();
        } catch (e) {
          // Ignore
        }
      }
    }
    
    // No warm sandbox available, create new one
    logger.info({ message: 'No warm sandbox available, creating new one' });
    const newSandbox = await this.createWarmSandbox();
    newSandbox.inUse = true;
    return newSandbox;
  }

  async executeCode(code, options = {}) {
    const {
      language = 'python',
      timeoutMs = 30000,
      requestId = null
    } = options;
    
    const startTime = Date.now();
    let sandboxEntry = null;
    
    try {
      this.metrics.totalExecutions++;
      
      // Get a sandbox
      sandboxEntry = await this.getSandbox();
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
      
      const executionTime = Date.now() - startTime;
      this.metrics.successfulExecutions++;
      this.metrics.avgExecutionTime = 
        (this.metrics.avgExecutionTime * (this.metrics.successfulExecutions - 1) + executionTime) / 
        this.metrics.successfulExecutions;
      
      // Determine success
      const isSuccess = !result.error;
      
      logger.info({ 
        requestId,
        executionTime,
        success: isSuccess,
        message: 'Code execution completed',
        warmPoolSize: this.warmSandboxes.length
      });
      
      // Format the result
      return {
        success: isSuccess,
        stdout: result.logs?.stdout?.join('\n') || '',
        stderr: result.logs?.stderr?.join('\n') || result.error || '',
        exitCode: result.error ? 1 : 0,
        executionTime
      };
      
    } catch (error) {
      this.metrics.failures++;
      logger.error({ 
        requestId,
        error: error.message,
        message: 'Code execution failed'
      });
      
      return {
        success: false,
        stdout: '',
        stderr: `E2B execution failed: ${error.message}`,
        exitCode: 1,
        executionTime: Date.now() - startTime
      };
      
    } finally {
      // Always return sandbox to pool
      if (sandboxEntry) {
        await this.returnSandbox(sandboxEntry);
      }
    }
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

  startMonitoring() {
    // Log metrics every 30 seconds
    setInterval(() => {
      const warmCount = this.warmSandboxes.length;
      const availableCount = this.warmSandboxes.filter(s => !s.inUse).length;
      
      logger.info({
        message: 'E2B Manager Status',
        warmSandboxes: warmCount,
        availableSandboxes: availableCount,
        metrics: this.getMetrics(),
        circuitBreaker: this.circuitBreaker.state,
        queueSize: this.createQueue.size,
        queuePending: this.createQueue.pending
      });
    }, 30000);
  }

  async shutdown() {
    logger.info({ message: 'Shutting down E2B manager' });
    
    // Clear the queue
    this.createQueue.clear();
    
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
    const warmCount = this.warmSandboxes.length;
    const availableCount = this.warmSandboxes.filter(s => !s.inUse).length;
    
    return {
      ...this.metrics,
      warmSandboxCount: warmCount,
      availableSandboxCount: availableCount,
      successRate: this.metrics.totalExecutions > 0 
        ? ((this.metrics.successfulExecutions / this.metrics.totalExecutions) * 100).toFixed(1) + '%'
        : 'N/A',
      creationSuccessRate: this.metrics.creationAttempts > 0
        ? ((this.metrics.creationSuccesses / this.metrics.creationAttempts) * 100).toFixed(1) + '%'
        : 'N/A'
    };
  }
}

// Export singleton instance
module.exports = new E2BManagerImproved();

// Graceful shutdown
process.on('SIGINT', async () => {
  await module.exports.shutdown();
  process.exit(0);
});

process.on('SIGTERM', async () => {
  await module.exports.shutdown();
  process.exit(0);
});