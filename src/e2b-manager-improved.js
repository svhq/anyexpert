// e2b-manager-improved.js - Better error handling and retry logic
const { Sandbox } = require('@e2b/code-interpreter');
const logger = require('./utils/logger');

class E2BManagerImproved {
  constructor() {
    this.templateId = process.env.E2B_TEMPLATE_ID || 'prod-all';
    
    // Improved retry configuration
    this.maxRetries = 3;
    this.initialRetryDelay = 1000;
    this.maxRetryDelay = 8000;
    this.backoffMultiplier = 2;
    
    // Better tracking
    this.metrics = {
      totalExecutions: 0,
      successes: 0,
      failures: 0,
      fallbacks: 0,
      retries: 0,
      avgExecutionTime: 0
    };
  }

  async executeCode(code, options = {}) {
    const { language = 'python', timeoutMs = 30000, requestId = null } = options;
    const startTime = Date.now();
    
    let lastError = null;
    let attempt = 0;
    
    while (attempt <= this.maxRetries) {
      try {
        if (attempt > 0) {
          const delay = Math.min(
            this.initialRetryDelay * Math.pow(this.backoffMultiplier, attempt - 1),
            this.maxRetryDelay
          );
          
          logger.info({ 
            requestId, 
            attempt, 
            delay,
            message: 'Retrying E2B execution with exponential backoff' 
          });
          
          await new Promise(resolve => setTimeout(resolve, delay));
          this.metrics.retries++;
        }
        
        // Create sandbox with better error handling
        let sandbox;
        try {
          sandbox = await Sandbox.create({
            template: this.templateId,
            timeoutMs: 30000 // 30s timeout for creation
          });
        } catch (createError) {
          logger.error({
            requestId,
            attempt,
            error: createError.message,
            message: 'Failed to create E2B sandbox'
          });
          
          // Check for specific error types
          if (createError.message.includes('rate limit')) {
            lastError = new Error('E2B rate limit exceeded. Please try again later.');
          } else if (createError.message.includes('timeout')) {
            lastError = new Error('E2B sandbox creation timed out.');
          } else {
            lastError = createError;
          }
          
          throw lastError;
        }
        
        // Execute code with proper timeout handling
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
        
        const result = await sandbox.runCode(wrappedCode, {
          language,
          timeoutMs: timeoutMs + 5000
        });
        
        // Clean up sandbox
        try {
          await sandbox.kill();
        } catch (e) {
          // Ignore cleanup errors
        }
        
        const executionTime = Date.now() - startTime;
        this.metrics.totalExecutions++;
        
        if (!result.error) {
          this.metrics.successes++;
          this.metrics.avgExecutionTime = 
            (this.metrics.avgExecutionTime * (this.metrics.successes - 1) + executionTime) / 
            this.metrics.successes;
          
          logger.info({ 
            requestId,
            attempt,
            executionTime,
            message: 'E2B execution successful'
          });
          
          return {
            success: true,
            stdout: result.logs?.stdout?.join('\n') || '',
            stderr: result.logs?.stderr?.join('\n') || '',
            exitCode: 0,
            executionTime,
            e2b: true
          };
        } else {
          // Code execution error (not E2B failure)
          return {
            success: false,
            stdout: result.logs?.stdout?.join('\n') || '',
            stderr: result.logs?.stderr?.join('\n') || result.error || '',
            exitCode: 1,
            executionTime,
            e2b: true
          };
        }
        
      } catch (error) {
        lastError = error;
        attempt++;
        
        logger.error({ 
          requestId,
          attempt,
          error: error.message,
          message: 'E2B execution attempt failed'
        });
      }
    }
    
    // All retries exhausted
    this.metrics.failures++;
    
    const errorMessage = lastError?.message || 'Unknown E2B error';
    logger.error({
      requestId,
      totalAttempts: attempt,
      error: errorMessage,
      message: 'E2B execution failed after all retries'
    });
    
    return {
      success: false,
      stdout: '',
      stderr: `E2B Error: ${errorMessage}\n\nThe code execution service is temporarily unavailable. This may be due to high demand or network issues.`,
      exitCode: 1,
      executionTime: Date.now() - startTime,
      e2b: false,
      e2bError: errorMessage
    };
  }

  async executeWithFallback(code, options = {}) {
    // First try E2B
    const e2bResult = await this.executeCode(code, options);
    
    if (e2bResult.success || !options.fallbackToLocal) {
      return e2bResult;
    }
    
    // Log fallback
    this.metrics.fallbacks++;
    logger.warn({
      requestId: options.requestId,
      e2bError: e2bResult.e2bError,
      message: 'E2B failed, falling back to local Python'
    });
    
    // Try local Python
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
        fallback: true,
        e2bError: e2bResult.e2bError
      };
      
    } catch (localError) {
      return {
        success: false,
        stdout: '',
        stderr: `Both E2B and local execution failed.\nE2B: ${e2bResult.e2bError}\nLocal: ${localError.message}`,
        exitCode: 1,
        executionTime: Date.now() - startTime,
        fallback: true,
        bothFailed: true
      };
    }
  }

  getMetrics() {
    return {
      ...this.metrics,
      successRate: this.metrics.totalExecutions > 0 
        ? (this.metrics.successes / this.metrics.totalExecutions * 100).toFixed(1) + '%'
        : 'N/A',
      fallbackRate: this.metrics.totalExecutions > 0
        ? (this.metrics.fallbacks / this.metrics.totalExecutions * 100).toFixed(1) + '%'
        : 'N/A'
    };
  }
}

module.exports = E2BManagerImproved;