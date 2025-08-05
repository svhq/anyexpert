// intelligent-queue.js - Fair request queuing with priority management
const PQueue = require('p-queue').default;

class IntelligentQueue {
  constructor(config = {}) {
    this.config = {
      globalConcurrency: config.globalConcurrency || 6,
      perUserConcurrency: config.perUserConcurrency || 2,
      maxQueueSize: config.maxQueueSize || 50,
      maxWaitTime: config.maxWaitTime || 30000, // 30 seconds
      ...config
    };
    
    // Global queue for all requests
    this.globalQueue = new PQueue({ 
      concurrency: this.config.globalConcurrency,
      timeout: this.config.maxWaitTime,
      throwOnTimeout: true
    });
    
    // Per-user queues for fairness
    this.userQueues = new Map();
    
    // Track user activity for priority calculation
    this.userActivity = new Map();
    
    // Metrics
    this.metrics = {
      enqueued: 0,
      completed: 0,
      failed: 0,
      timedOut: 0,
      rejected: 0
    };
    
    // Clean up old user data periodically
    const cleanupInterval = config.cleanupInterval !== undefined ? config.cleanupInterval : 60000;
    if (cleanupInterval > 0) {
      this.cleanupInterval = setInterval(() => this.cleanup(), cleanupInterval);
    }
  }
  
  async enqueue(userId, task, options = {}) {
    const {
      complexity = 0.5,
      priority: userPriority,
      metadata = {}
    } = options;
    
    this.metrics.enqueued++;
    
    // Check global queue size
    if (this.globalQueue.size >= this.config.maxQueueSize) {
      this.metrics.rejected++;
      throw new Error(`Queue full. Please try again in a few moments. Current queue size: ${this.globalQueue.size}`);
    }
    
    // Get or create user queue
    if (!this.userQueues.has(userId)) {
      this.userQueues.set(userId, new PQueue({ 
        concurrency: this.config.perUserConcurrency 
      }));
      this.userActivity.set(userId, {
        requests: [],
        totalComplexity: 0
      });
    }
    
    const userQueue = this.userQueues.get(userId);
    const userActivityData = this.userActivity.get(userId);
    
    // Track request
    const request = {
      timestamp: Date.now(),
      complexity,
      metadata
    };
    
    userActivityData.requests.push(request);
    userActivityData.totalComplexity += complexity;
    
    // Calculate dynamic priority
    const priority = userPriority !== undefined 
      ? userPriority 
      : this.calculatePriority(userId, complexity, userActivityData);
    
    // Create wrapped task that goes through both queues
    const wrappedTask = async () => {
      const startTime = Date.now();
      
      try {
        // Add to global queue with calculated priority
        const result = await this.globalQueue.add(
          async () => {
            // Check if request has waited too long
            const waitTime = Date.now() - startTime;
            if (waitTime > this.config.maxWaitTime) {
              throw new Error(`Request timed out after ${waitTime}ms wait`);
            }
            
            return await task();
          },
          { priority }
        );
        
        this.metrics.completed++;
        
        // Update user activity
        const executionTime = Date.now() - startTime;
        request.executionTime = executionTime;
        request.completed = true;
        
        return result;
        
      } catch (error) {
        if (error.message.includes('timed out')) {
          this.metrics.timedOut++;
        } else {
          this.metrics.failed++;
        }
        
        request.failed = true;
        request.error = error.message;
        
        throw error;
      }
    };
    
    // Add to user queue (ensures per-user concurrency)
    return userQueue.add(wrappedTask);
  }
  
  calculatePriority(userId, complexity, userActivity) {
    const now = Date.now();
    const recentWindow = 60000; // 1 minute
    
    // Get recent requests
    const recentRequests = userActivity.requests.filter(
      r => (now - r.timestamp) < recentWindow
    );
    
    // Calculate factors
    const userLoad = recentRequests.length / 10; // Normalize to 0-1 (10 requests/min = 1.0)
    const avgComplexity = recentRequests.length > 0
      ? recentRequests.reduce((sum, r) => sum + r.complexity, 0) / recentRequests.length
      : 0.5;
    
    // Calculate wait bonus (increases priority for waiting requests)
    const oldestPending = recentRequests.find(r => !r.completed && !r.failed);
    const waitBonus = oldestPending 
      ? Math.min((now - oldestPending.timestamp) / 30000, 1) // Max bonus at 30s
      : 0;
    
    // Priority formula (higher = better)
    // - Prefer simple queries (40%)
    // - Prefer users with low recent activity (30%)
    // - Boost long-waiting requests (30%)
    const priority = (
      (1 - complexity) * 0.4 +
      (1 - Math.min(userLoad, 1)) * 0.3 +
      waitBonus * 0.3
    );
    
    return priority;
  }
  
  getUserLoad(userId) {
    const userActivity = this.userActivity.get(userId);
    if (!userActivity) return 0;
    
    const now = Date.now();
    const recentWindow = 60000; // 1 minute
    
    const recentRequests = userActivity.requests.filter(
      r => (now - r.timestamp) < recentWindow
    );
    
    return recentRequests.length;
  }
  
  cleanup() {
    const now = Date.now();
    const staleThreshold = 5 * 60000; // 5 minutes
    
    // Clean up old user data
    for (const [userId, activity] of this.userActivity.entries()) {
      // Remove old requests
      activity.requests = activity.requests.filter(
        r => (now - r.timestamp) < staleThreshold
      );
      
      // Remove user if no recent activity
      if (activity.requests.length === 0) {
        this.userActivity.delete(userId);
        const userQueue = this.userQueues.get(userId);
        if (userQueue && userQueue.size === 0) {
          this.userQueues.delete(userId);
        }
      }
    }
  }
  
  getMetrics() {
    return {
      ...this.metrics,
      globalQueueSize: this.globalQueue.size,
      globalQueuePending: this.globalQueue.pending,
      activeUsers: this.userQueues.size,
      successRate: this.metrics.enqueued > 0
        ? ((this.metrics.completed / this.metrics.enqueued) * 100).toFixed(1) + '%'
        : 'N/A',
      timeoutRate: this.metrics.enqueued > 0
        ? ((this.metrics.timedOut / this.metrics.enqueued) * 100).toFixed(1) + '%'
        : 'N/A',
      rejectionRate: this.metrics.enqueued > 0
        ? ((this.metrics.rejected / (this.metrics.enqueued + this.metrics.rejected)) * 100).toFixed(1) + '%'
        : 'N/A',
      userLoads: Array.from(this.userActivity.keys()).map(userId => ({
        userId,
        load: this.getUserLoad(userId)
      })).sort((a, b) => b.load - a.load)
    };
  }
  
  async shutdown() {
    if (this.cleanupInterval) {
      clearInterval(this.cleanupInterval);
    }
    this.globalQueue.clear();
    
    for (const queue of this.userQueues.values()) {
      queue.clear();
    }
  }
}

module.exports = IntelligentQueue;