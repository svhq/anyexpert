// pool-manager.js - Dynamic pool management with load balancing
const { Sandbox } = require('@e2b/code-interpreter');
const SandboxExecutor = require('./sandbox-executor');
const logger = require('../utils/logger');
const templateConfig = require('../../config/template-config');

class DynamicPoolManager {
  constructor(rateLimiter, config = {}) {
    this.rateLimiter = rateLimiter;
    this.pool = new Map(); // sandboxId -> SandboxExecutor
    
    // Configuration
    this.config = {
      templateId: templateConfig.TEMPLATE_ID,
      expectedBuildTag: templateConfig.EXPECTED_BUILD_TAG,
      minSize: config.minSize || 2,
      targetSize: config.targetSize || 3,
      maxSize: config.maxSize || 10,
      sandboxTimeout: config.sandboxTimeout || 30 * 60 * 1000, // 30 minutes
      maxReuseCount: config.maxReuseCount || 50,
      maxAge: config.maxAge || 60 * 60 * 1000, // 1 hour
      scaleUpThreshold: config.scaleUpThreshold || 0.7,
      scaleDownThreshold: config.scaleDownThreshold || 0.3,
      maintenanceInterval: config.maintenanceInterval !== undefined ? config.maintenanceInterval : 30000
    };
    
    // Metrics
    this.metrics = {
      created: 0,
      destroyed: 0,
      reused: 0,
      errors: 0
    };
    
    // Start maintenance
    if (this.config.maintenanceInterval > 0) {
      this.maintenanceInterval = setInterval(() => this.maintain(), this.config.maintenanceInterval);
    }
  }
  
  async initialize() {
    // Log template configuration
    templateConfig.logTemplateInfo();
    
    logger.info({ 
      message: 'Initializing dynamic pool',
      targetSize: this.config.targetSize,
      expectedBuildTag: this.config.expectedBuildTag
    });
    
    // First, try to discover and reuse existing sandboxes
    try {
      const existingSandboxes = await Sandbox.list();
      logger.info({ 
        message: 'Found existing sandboxes',
        count: existingSandboxes.length 
      });
      
      // Connect to existing sandboxes up to our target size
      let connected = 0;
      for (const sandboxInfo of existingSandboxes) {
        if (connected >= this.config.targetSize) break;
        
        try {
          const sandbox = await Sandbox.connect(sandboxInfo.sandboxId);
          const executor = new SandboxExecutor(sandbox, sandboxInfo.sandboxId);
          
          // Health check including BUILD_TAG verification
          const healthy = await executor.healthCheck();
          const buildTag = await executor.getBuildTag();
          
          if (healthy && buildTag === this.config.expectedBuildTag) {
            this.pool.set(sandboxInfo.sandboxId, executor);
            connected++;
            logger.info({ 
              message: 'Connected to existing sandbox',
              sandboxId: sandboxInfo.sandboxId,
              buildTag: buildTag,
              poolSize: this.pool.size
            });
          } else {
            await sandbox.kill();
            logger.warn({ 
              message: 'Existing sandbox rejected',
              sandboxId: sandboxInfo.sandboxId,
              reason: !healthy ? 'unhealthy' : 'wrong BUILD_TAG',
              actualBuildTag: buildTag,
              expectedBuildTag: this.config.expectedBuildTag
            });
          }
        } catch (error) {
          logger.warn({ 
            message: 'Failed to connect to existing sandbox',
            sandboxId: sandboxInfo.sandboxId,
            error: error.message
          });
        }
      }
      
      // If we connected to enough sandboxes, we're done
      if (this.pool.size >= this.config.minSize) {
        logger.info({ 
          message: 'Successfully reused existing sandboxes',
          poolSize: this.pool.size
        });
        return;
      }
    } catch (error) {
      logger.warn({ 
        message: 'Failed to list existing sandboxes',
        error: error.message
      });
    }
    
    // Create new sandboxes only if needed
    const sandboxesNeeded = Math.max(0, this.config.minSize - this.pool.size);
    for (let i = 0; i < sandboxesNeeded; i++) {
      if (i > 0 || this.pool.size > 0) {
        await new Promise(resolve => setTimeout(resolve, 500 + Math.random() * 500));
      }
      
      try {
        await this.createSandbox();
      } catch (error) {
        logger.error({ 
          message: 'Failed to create initial sandbox',
          error: error.message,
          attempt: i + 1
        });
      }
    }
  }
  
  async createSandbox() {
    // Check rate limit - wait up to 10 seconds for creation token
    await this.rateLimiter.acquireToken('creation', 1, 10000);
    
    this.metrics.created++;
    
    const sandbox = await Sandbox.create(this.config.templateId, {
      timeoutMs: this.config.sandboxTimeout
    });
    
    const executor = new SandboxExecutor(sandbox, sandbox.id);
    
    // Health check
    const healthy = await executor.healthCheck();
    if (!healthy) {
      await executor.kill();
      throw new Error('Sandbox failed health check');
    }
    
    this.pool.set(sandbox.id, executor);
    
    logger.info({ 
      message: 'Created sandbox',
      sandboxId: sandbox.id,
      poolSize: this.pool.size
    });
    
    return executor;
  }
  
  async getExecutor(userId) {
    // Get all executors sorted by score
    const executors = Array.from(this.pool.values())
      .filter(e => e.state !== 'unhealthy')
      .sort((a, b) => b.score - a.score);
    
    // Find best available executor
    for (const executor of executors) {
      if (executor.load < 2) {
        this.metrics.reused++;
        return executor;
      }
    }
    
    // Check if we should scale up
    const avgLoad = this.getAverageLoad();
    if (this.pool.size < this.config.maxSize && avgLoad > this.config.scaleUpThreshold) {
      try {
        const newExecutor = await this.createSandbox();
        return newExecutor;
      } catch (error) {
        logger.warn({ 
          message: 'Failed to scale up',
          error: error.message,
          currentSize: this.pool.size
        });
      }
    }
    
    // Return least loaded executor (might be busy)
    if (executors.length > 0) {
      this.metrics.reused++;
      return executors[0];
    }
    
    // Emergency: create new sandbox
    if (this.pool.size === 0) {
      return await this.createSandbox();
    }
    
    throw new Error('No available executors');
  }
  
  async destroySandbox(sandboxId) {
    const executor = this.pool.get(sandboxId);
    if (!executor) return;
    
    this.pool.delete(sandboxId);
    this.metrics.destroyed++;
    
    try {
      await executor.kill();
    } catch (error) {
      logger.error({ 
        message: 'Error killing sandbox',
        sandboxId,
        error: error.message
      });
    }
  }
  
  async maintain() {
    const now = Date.now();
    const avgLoad = this.getAverageLoad();
    
    // Clean up unhealthy or old sandboxes
    for (const [id, executor] of this.pool.entries()) {
      const shouldRemove = 
        executor.state === 'unhealthy' ||
        executor.age > this.config.maxAge ||
        executor.executionCount > this.config.maxReuseCount;
      
      if (shouldRemove) {
        logger.info({ 
          message: 'Removing sandbox',
          sandboxId: id,
          reason: executor.state === 'unhealthy' ? 'unhealthy' : 'age/reuse'
        });
        await this.destroySandbox(id);
      }
    }
    
    // Scale down if load is low
    if (avgLoad < this.config.scaleDownThreshold && this.pool.size > this.config.minSize) {
      // Find most idle sandbox
      const mostIdle = Array.from(this.pool.values())
        .sort((a, b) => b.idleTime - a.idleTime)[0];
      
      if (mostIdle && mostIdle.idleTime > 5 * 60 * 1000) { // 5 minutes idle
        logger.info({ 
          message: 'Scaling down - removing idle sandbox',
          sandboxId: mostIdle.sandboxId,
          idleTime: mostIdle.idleTime
        });
        await this.destroySandbox(mostIdle.sandboxId);
      }
    }
    
    // Scale up if needed
    if (avgLoad > this.config.scaleUpThreshold && this.pool.size < this.config.targetSize) {
      try {
        await this.createSandbox();
      } catch (error) {
        logger.warn({ 
          message: 'Failed to scale up during maintenance',
          error: error.message
        });
      }
    }
    
    // Log pool status
    logger.info({
      message: 'Pool maintenance complete',
      poolSize: this.pool.size,
      avgLoad,
      metrics: this.getMetrics()
    });
  }
  
  getAverageLoad() {
    if (this.pool.size === 0) return 0;
    
    const totalLoad = Array.from(this.pool.values())
      .reduce((sum, executor) => sum + executor.load, 0);
    
    return totalLoad / this.pool.size;
  }
  
  getMetrics() {
    const executors = Array.from(this.pool.values());
    
    return {
      poolSize: this.pool.size,
      avgLoad: this.getAverageLoad(),
      healthyCount: executors.filter(e => e.state !== 'unhealthy').length,
      busyCount: executors.filter(e => e.state === 'busy').length,
      totalExecutions: executors.reduce((sum, e) => sum + e.executionCount, 0),
      avgQueueDepth: executors.reduce((sum, e) => sum + e.queue.size, 0) / Math.max(1, this.pool.size),
      ...this.metrics
    };
  }
  
  async clearAllSandboxes() {
    logger.info({ message: 'Clearing all sandboxes from E2B' });
    
    try {
      // First kill all sandboxes in our pool
      const poolPromises = Array.from(this.pool.keys()).map(id => this.destroySandbox(id));
      await Promise.all(poolPromises);
      
      // Then kill any remaining sandboxes in E2B
      const allSandboxes = await Sandbox.list();
      logger.info({ 
        message: 'Found sandboxes to clear',
        count: allSandboxes.length 
      });
      
      for (const sandboxInfo of allSandboxes) {
        try {
          const sandbox = await Sandbox.connect(sandboxInfo.sandboxId);
          await sandbox.kill();
          logger.info({ 
            message: 'Killed sandbox',
            sandboxId: sandboxInfo.sandboxId
          });
        } catch (error) {
          logger.warn({ 
            message: 'Failed to kill sandbox',
            sandboxId: sandboxInfo.sandboxId,
            error: error.message
          });
        }
      }
      
      logger.info({ message: 'All sandboxes cleared' });
    } catch (error) {
      logger.error({ 
        message: 'Error clearing sandboxes',
        error: error.message
      });
    }
  }
  
  async shutdown() {
    if (this.maintenanceInterval) {
      clearInterval(this.maintenanceInterval);
    }
    
    // Kill all sandboxes
    const promises = Array.from(this.pool.keys()).map(id => this.destroySandbox(id));
    await Promise.all(promises);
    
    logger.info({ 
      message: 'Pool manager shutdown complete',
      metrics: this.metrics
    });
  }
}

module.exports = DynamicPoolManager;