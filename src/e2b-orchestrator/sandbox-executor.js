// sandbox-executor.js - Per-sandbox execution queue with fairness
const PQueue = require('p-queue').default || require('p-queue');

class SandboxExecutor {
  constructor(sandbox, sandboxId) {
    this.sandbox = sandbox;
    this.sandboxId = sandboxId;
    this.queue = new PQueue({ concurrency: 1 }); // Enforce serial execution
    this.lastUsed = Date.now();
    this.createdAt = Date.now();
    this.executionCount = 0;
    this.totalExecutionTime = 0;
    this.errors = 0;
    
    // Track current state
    this.state = 'idle'; // idle, busy, unhealthy
  }
  
  async execute(code, options = {}) {
    const { 
      priority = 0,
      timeoutMs = 30000,
      language = 'python'
    } = options;
    
    if (this.state === 'unhealthy') {
      throw new Error(`Sandbox ${this.sandboxId} is unhealthy`);
    }
    
    const task = async () => {
      this.state = 'busy';
      const startTime = Date.now();
      
      try {
        // Add timeout protection for Python
        const wrappedCode = language === 'python' ? this.wrapPythonCode(code, timeoutMs) : code;
        
        const result = await this.sandbox.runCode(wrappedCode, {
          language,
          timeoutMs: timeoutMs + 5000 // Give SDK extra time
        });
        
        const executionTime = Date.now() - startTime;
        this.executionCount++;
        this.totalExecutionTime += executionTime;
        this.lastUsed = Date.now();
        
        // Check if result indicates success
        const success = !result.error;
        
        if (!success) {
          this.errors++;
          // Mark unhealthy if too many errors
          if (this.errors > 5) {
            this.state = 'unhealthy';
          }
        }
        
        return {
          success,
          stdout: result.logs?.stdout?.join('\n') || '',
          stderr: result.logs?.stderr?.join('\n') || result.error || '',
          exitCode: result.error ? 1 : 0,
          executionTime,
          sandboxId: this.sandboxId,
          queueDepth: this.queue.size
        };
        
      } catch (error) {
        this.errors++;
        if (this.errors > 5) {
          this.state = 'unhealthy';
        }
        
        throw error;
        
      } finally {
        if (this.queue.size === 0) {
          this.state = 'idle';
        }
      }
    };
    
    return this.queue.add(task, { priority });
  }
  
  wrapPythonCode(code, timeoutMs) {
    // Auto-install wrapper disabled - numpy_financial is pre-installed in template
    // Previously: auto-installed numpy_financial if ImportError
    // Now: expecting it to be available via template BUILD_TAG npf-20250805-003130

    return `
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
`;
  }
  
  async healthCheck() {
    try {
      // Check both health and BUILD_TAG to ensure we have the latest template
      const healthCode = `
import os
print("__healthy__")
print(f"BUILD_TAG:{os.environ.get('BUILD_TAG', 'missing')}")
# Check for numpy_financial
try:
    import numpy_financial
    print("NPF:installed")
except ImportError:
    print("NPF:missing")
`;
      
      const result = await this.sandbox.runCode(healthCode, { 
        timeoutMs: 5000,
        language: 'python'
      });
      
      const output = result.logs?.stdout?.join('') || '';
      const hasHealthy = output.includes('__healthy__');
      const hasBuildTag = output.includes('BUILD_TAG:') && !output.includes('BUILD_TAG:missing');
      const hasNumpyFinancial = output.includes('NPF:installed');
      
      // Sandbox is healthy if basic health check passes
      // TODO: Re-enable numpy_financial check once template is updated
      const healthy = !result.error && hasHealthy;
      
      if (healthy) {
        this.state = 'idle';
        this.errors = 0;
      } else {
        this.state = 'unhealthy';
        logger.warn({
          message: 'Sandbox failed health check',
          sandboxId: this.sandboxId,
          hasHealthy,
          hasBuildTag,
          hasNumpyFinancial,
          output
        });
      }
      
      return healthy;
    } catch (error) {
      this.state = 'unhealthy';
      return false;
    }
  }
  
  async getBuildTag() {
    // Return cached BUILD_TAG if available
    if (this.buildTag) {
      return this.buildTag;
    }
    
    // Otherwise fetch it
    try {
      const code = `import os\nprint(os.environ.get('BUILD_TAG', 'missing'))`;
      const result = await this.sandbox.runCode(code, { 
        timeoutMs: 5000,
        language: 'python'
      });
      
      if (!result.error && result.logs?.stdout) {
        this.buildTag = result.logs.stdout.join('').trim();
        return this.buildTag;
      }
    } catch (error) {
      logger.error({ 
        message: 'Error getting BUILD_TAG',
        sandboxId: this.sandboxId,
        error: error.message
      });
    }
    
    return 'missing';
  }
  
  async kill() {
    try {
      await this.sandbox.kill();
      this.state = 'dead';
    } catch (error) {
      // Ignore errors during cleanup
    }
  }
  
  // Metrics and scoring
  get load() {
    return this.queue.size + (this.state === 'busy' ? 1 : 0);
  }
  
  get idleTime() {
    return Date.now() - this.lastUsed;
  }
  
  get age() {
    return Date.now() - this.createdAt;
  }
  
  get avgExecutionTime() {
    return this.executionCount > 0 
      ? this.totalExecutionTime / this.executionCount 
      : 0;
  }
  
  get errorRate() {
    return this.executionCount > 0 
      ? this.errors / this.executionCount 
      : 0;
  }
  
  // Score for load balancing (higher = better candidate)
  get score() {
    if (this.state === 'unhealthy') return -1;
    if (this.state === 'busy' && this.queue.size > 2) return 0;
    
    // Factors:
    // - Prefer idle sandboxes (weight: 0.4)
    // - Prefer low queue depth (weight: 0.3)
    // - Prefer low error rate (weight: 0.2)
    // - Slight preference for newer sandboxes (weight: 0.1)
    
    const idleScore = Math.min(1, this.idleTime / 60000); // Max at 1 minute idle
    const queueScore = 1 / (this.load + 1);
    const errorScore = 1 - this.errorRate;
    const ageScore = Math.max(0, 1 - (this.age / 3600000)); // Decay over 1 hour
    
    return (
      idleScore * 0.4 +
      queueScore * 0.3 +
      errorScore * 0.2 +
      ageScore * 0.1
    );
  }
  
  getMetrics() {
    return {
      sandboxId: this.sandboxId,
      state: this.state,
      load: this.load,
      queueSize: this.queue.size,
      executionCount: this.executionCount,
      avgExecutionTime: this.avgExecutionTime,
      errorRate: this.errorRate,
      idleTime: this.idleTime,
      age: this.age,
      score: this.score
    };
  }
}

module.exports = SandboxExecutor;