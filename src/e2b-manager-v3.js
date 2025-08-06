// e2b-manager-v3.js - E2B Manager with fail-fast for missing API key
const e2bManager = require('./e2b-manager');
const logger = require('./utils/logger');

class E2BManagerV3 {
  constructor() {
    this.e2bManager = e2bManager;
    this.templateId = process.env.E2B_TEMPLATE_ID || 'prod-all';
    
    // For backward compatibility
    this.maxRetries = process.env.NODE_ENV === 'production' ? 0 : 2;
    this.sandboxTimeout = 30 * 60 * 1000;
  }
  
  async executeCode(code, options = {}) {
    const {
      language = 'python',
      timeoutMs = 30000,
      requestId = null
    } = options;
    
    // Fail fast if E2B is not configured
    if (!process.env.E2B_API_KEY) {
      logger.warn({ requestId, message: 'E2B_API_KEY not configured, skipping code execution' });
      return {
        success: false,
        stdout: '',
        stderr: 'E2B not configured. Using reasoning instead.',
        exitCode: 1,
        executionTime: 0,
        skipped: true
      };
    }
    
    const userId = options.userId || 'default';
    
    try {
      // Delegate to the actual e2b-manager
      const result = await this.e2bManager.executeCode(code, {
        language,
        timeoutMs,
        requestId,
        userId
      });
      
      // Ensure backward compatibility with expected format
      return {
        success: result.success,
        stdout: result.stdout || '',
        stderr: result.stderr || '',
        exitCode: result.exitCode || 0,
        executionTime: result.executionTime || 0,
        cached: result.cached || false
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
    return this.e2bManager.getMetrics();
  }
  
  async shutdown() {
    await this.e2bManager.shutdown();
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