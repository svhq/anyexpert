// token-bucket.js - Adaptive rate limiter with token bucket algorithm
const logger = require('../utils/logger');

class TokenBucket {
  constructor(capacity, refillRate, refillInterval = 1000) {
    this.capacity = capacity;
    this.tokens = capacity;
    this.refillRate = refillRate;
    this.refillInterval = refillInterval;
    this.lastRefill = Date.now();
    
    // Start refill timer
    this.refillTimer = setInterval(() => this.refill(), this.refillInterval);
  }
  
  refill() {
    const now = Date.now();
    const timePassed = now - this.lastRefill;
    const tokensToAdd = (timePassed / this.refillInterval) * this.refillRate;
    
    this.tokens = Math.min(this.capacity, this.tokens + tokensToAdd);
    this.lastRefill = now;
  }
  
  async tryAcquire(count = 1) {
    this.refill();
    
    if (this.tokens >= count) {
      this.tokens -= count;
      return true;
    }
    
    return false;
  }
  
  async acquire(count = 1) {
    // Try immediately
    if (await this.tryAcquire(count)) {
      return true;
    }
    
    // Wait for tokens
    const waitTime = ((count - this.tokens) / this.refillRate) * this.refillInterval;
    await new Promise(resolve => setTimeout(resolve, waitTime));
    
    return this.tryAcquire(count);
  }
  
  getAvailableTokens() {
    this.refill();
    return Math.floor(this.tokens);
  }
  
  destroy() {
    clearInterval(this.refillTimer);
  }
}

class AdaptiveRateLimiter {
  constructor(config = {}) {
    const { 
      hobby = false,
      sharedIP = true 
    } = config;
    
    // E2B documented limits
    this.buckets = {
      // Sandbox creation: 1/sec (Hobby) or 5/sec (Pro)
      creation: new TokenBucket(
        hobby ? 1 : 5,
        hobby ? 1 : 5,
        1000
      ),
      
      // Operations: 40k/60s = 667/sec per IP
      operations: new TokenBucket(
        40000,
        667,
        1000
      ),
      
      // Lifecycle API: 20k/30s = 667/sec
      lifecycle: new TokenBucket(
        20000,
        667,
        1000
      )
    };
    
    // Per-sandbox serialization tracking
    this.sandboxQueues = new Map();
    
    // Adaptive health score
    this.healthScore = 1.0;
    this.failureWindow = [];
    this.windowSize = 60000; // 1 minute window
  }
  
  async acquireToken(type, count = 1, maxWaitMs = 5000) {
    const bucket = this.buckets[type];
    if (!bucket) {
      throw new Error(`Unknown rate limit type: ${type}`);
    }
    
    // Apply adaptive rate based on health
    const adaptedCount = count / this.healthScore;
    
    // First try to acquire immediately
    const immediateAcquired = await bucket.tryAcquire(adaptedCount);
    if (immediateAcquired) {
      this.recordSuccess();
      return true;
    }
    
    // Calculate wait time needed
    const waitTime = ((adaptedCount - bucket.getAvailableTokens()) / bucket.refillRate) * bucket.refillInterval;
    
    if (waitTime > maxWaitMs) {
      this.recordFailure();
      throw new Error(`Rate limit exceeded for ${type}. Would need to wait ${Math.ceil(waitTime)}ms. Available: ${bucket.getAvailableTokens()}`);
    }
    
    // Wait for tokens to become available
    logger.debug({ 
      message: 'Waiting for rate limit tokens',
      type,
      waitTime: Math.ceil(waitTime),
      available: bucket.getAvailableTokens()
    });
    
    const acquired = await bucket.acquire(adaptedCount);
    
    if (acquired) {
      this.recordSuccess();
      return true;
    } else {
      this.recordFailure();
      throw new Error(`Rate limit exceeded for ${type} after waiting. Available: ${bucket.getAvailableTokens()}`);
    }
  }
  
  async acquireForSandbox(sandboxId, operation = 'runCode') {
    // Ensure per-sandbox serialization
    if (!this.sandboxQueues.has(sandboxId)) {
      this.sandboxQueues.set(sandboxId, {
        busy: false,
        queueDepth: 0
      });
    }
    
    const sandboxState = this.sandboxQueues.get(sandboxId);
    
    if (sandboxState.busy) {
      throw new Error(`Sandbox ${sandboxId} is busy. Cannot execute concurrent operations.`);
    }
    
    // Acquire operation token
    await this.acquireToken('operations');
    
    // Mark sandbox as busy
    sandboxState.busy = true;
    
    return {
      release: () => {
        sandboxState.busy = false;
      }
    };
  }
  
  recordSuccess() {
    this.cleanWindow();
    
    // Gradually improve health score
    this.healthScore = Math.min(1.0, this.healthScore * 1.02);
  }
  
  recordFailure() {
    this.cleanWindow();
    this.failureWindow.push(Date.now());
    
    // Calculate failure rate
    const failureRate = this.failureWindow.length / (this.windowSize / 1000);
    
    // Adjust health score based on failure rate
    if (failureRate > 10) {
      this.healthScore = Math.max(0.1, this.healthScore * 0.5);
    } else if (failureRate > 5) {
      this.healthScore = Math.max(0.3, this.healthScore * 0.8);
    } else {
      this.healthScore = Math.max(0.5, this.healthScore * 0.95);
    }
  }
  
  cleanWindow() {
    const cutoff = Date.now() - this.windowSize;
    this.failureWindow = this.failureWindow.filter(time => time > cutoff);
  }
  
  getStatus() {
    this.cleanWindow();
    
    return {
      healthScore: this.healthScore,
      failureRate: this.failureWindow.length / (this.windowSize / 1000),
      availableTokens: {
        creation: this.buckets.creation.getAvailableTokens(),
        operations: this.buckets.operations.getAvailableTokens(),
        lifecycle: this.buckets.lifecycle.getAvailableTokens()
      },
      busySandboxes: Array.from(this.sandboxQueues.entries())
        .filter(([_, state]) => state.busy)
        .map(([id]) => id)
    };
  }
  
  destroy() {
    Object.values(this.buckets).forEach(bucket => bucket.destroy());
  }
}

module.exports = { TokenBucket, AdaptiveRateLimiter };