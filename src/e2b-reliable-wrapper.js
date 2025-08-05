const fetch = require('node-fetch');
const logger = require('./utils/logger');

/**
 * Reliable wrapper for E2B service that handles common failure modes
 */
class E2BReliableWrapper {
  constructor(serviceUrl = 'http://localhost:3001') {
    this.serviceUrl = serviceUrl;
    this.maxRetries = 2;
    this.healthCheckInterval = 60000; // Check health every minute
    this.lastHealthCheck = null;
    this.isHealthy = false;
  }

  /**
   * Check if E2B service is healthy
   */
  async checkHealth() {
    try {
      const response = await fetch(`${this.serviceUrl}/health`, {
        method: 'GET',
        signal: AbortSignal.timeout(3000)
      });
      
      if (!response.ok) {
        this.isHealthy = false;
        return false;
      }
      
      const data = await response.json();
      this.isHealthy = data.sandboxReady === true;
      this.lastHealthCheck = Date.now();
      return this.isHealthy;
    } catch (error) {
      this.isHealthy = false;
      return false;
    }
  }

  /**
   * Execute code with automatic retry and fallback
   */
  async runCode(language, source, timeout = 30000, requestId = null) {
    // Check health if needed
    if (!this.lastHealthCheck || Date.now() - this.lastHealthCheck > this.healthCheckInterval) {
      await this.checkHealth();
    }

    // If not healthy, return error immediately
    if (!this.isHealthy) {
      logger.warn({ requestId, message: 'E2B service not healthy, skipping code execution' });
      return {
        success: false,
        error: 'E2B service not available',
        stdout: '',
        stderr: 'E2B service is not healthy. Code execution skipped.',
        exitCode: 1
      };
    }

    // Try to execute code with retries
    for (let attempt = 0; attempt <= this.maxRetries; attempt++) {
      try {
        const response = await fetch(`${this.serviceUrl}/run_code`, {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ language, source, timeout }),
          signal: AbortSignal.timeout(timeout + 5000) // Add 5s buffer
        });

        if (!response.ok) {
          throw new Error(`HTTP ${response.status}: ${response.statusText}`);
        }

        const result = await response.json();
        
        // Check for sandbox death
        if (result.stderr && result.stderr.includes('Sandbox is probably not running')) {
          this.isHealthy = false;
          logger.warn({ requestId, attempt, message: 'E2B sandbox died' });
          
          // Try to restart sandbox
          if (attempt < this.maxRetries) {
            await this.restartSandbox();
            continue;
          }
        }

        return {
          success: result.exitCode === 0,
          ...result
        };

      } catch (error) {
        logger.error({ 
          requestId, 
          attempt, 
          error: error.message,
          isTimeout: error.name === 'AbortError'
        });

        if (attempt === this.maxRetries) {
          return {
            success: false,
            error: error.message,
            stdout: '',
            stderr: `E2B execution failed: ${error.message}`,
            exitCode: 1
          };
        }

        // Wait before retry
        await new Promise(resolve => setTimeout(resolve, 1000 * (attempt + 1)));
      }
    }
  }

  /**
   * Attempt to restart the sandbox
   */
  async restartSandbox() {
    try {
      logger.info({ message: 'Attempting to restart E2B sandbox' });
      
      const response = await fetch(`${this.serviceUrl}/restart_sandbox`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        signal: AbortSignal.timeout(10000)
      });

      if (response.ok) {
        // Wait for sandbox to initialize
        await new Promise(resolve => setTimeout(resolve, 5000));
        await this.checkHealth();
      }
    } catch (error) {
      logger.error({ message: 'Failed to restart sandbox', error: error.message });
    }
  }

  /**
   * Get service status
   */
  getStatus() {
    return {
      healthy: this.isHealthy,
      lastHealthCheck: this.lastHealthCheck,
      serviceUrl: this.serviceUrl
    };
  }
}

// Export singleton instance
module.exports = new E2BReliableWrapper();