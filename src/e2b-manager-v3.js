// e2b-manager-v3.js - E2B Manager using the new orchestrator
const orchestrator = require('./e2b-orchestrator');
const logger = require('./utils/logger');

class E2BManagerV3 {
  constructor() {
    this.orchestrator = orchestrator;
    
    // For backward compatibility
    this.maxRetries = 2;
    this.sandboxTimeout = 30 * 60 * 1000;
  }
  
  async executeCode(code, options = {}) {
    const {
      language = 'python',
      timeoutMs = 30000,
      requestId = null
    } = options;
    
    const userId = options.userId || 'default';
    
    try {
      const result = await this.orchestrator.execute(userId, code, {
        language,
        timeoutMs,
        requestId
      });
      
      // Ensure backward compatibility with expected format
      return {
        success: result.success,
        stdout: result.stdout || '',
        stderr: result.stderr || '',
        exitCode: result.exitCode || 0,
        executionTime: result.executionTime || result.totalTime || 0,
        cached: result.cached || false,
        executor: result.executor || 'unknown'
      };
      
    } catch (error) {
      logger.error({
        requestId,
        error: error.message,
        message: 'E2B execution failed'
      });
      
      return {
        success: false,
        stdout: '',
        stderr: error.message,
        exitCode: 1,
        executionTime: 0
      };
    }
  }
  
  async executeWithFallback(code, options = {}) {
    // The orchestrator already handles fallback through circuit breaker
    return this.executeCode(code, options);
  }
  
  getMetrics() {
    const metrics = this.orchestrator.getMetrics();
    
    // Flatten metrics for backward compatibility
    return {
      ...metrics.pool,
      ...metrics.orchestrator,
      rateLimiter: metrics.rateLimiter,
      queue: metrics.queue,
      circuitBreaker: metrics.circuitBreaker,
      cache: metrics.cache
    };
  }
  
  async shutdown() {
    await this.orchestrator.shutdown();
  }
}

// Export singleton instance
module.exports = new E2BManagerV3();

// Handle process termination
process.on('SIGINT', async () => {
  await module.exports.shutdown();
  process.exit(0);
});

process.on('SIGTERM', async () => {
  await module.exports.shutdown();
  process.exit(0);
});