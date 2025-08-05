// circuit-breaker.js - Gradual circuit breaker with adaptive recovery
class GradualCircuitBreaker {
  constructor(config = {}) {
    this.config = {
      failureThreshold: config.failureThreshold || 5,
      successThreshold: config.successThreshold || 3,
      timeout: config.timeout || 30000, // 30 seconds
      volumeThreshold: config.volumeThreshold || 10, // Min requests before opening
      halfOpenMaxRequests: config.halfOpenMaxRequests || 3,
      ...config
    };
    
    // States
    this.states = {
      CLOSED: 'CLOSED',
      OPEN: 'OPEN',
      HALF_OPEN: 'HALF_OPEN'
    };
    
    this.state = this.states.CLOSED;
    this.failures = 0;
    this.successes = 0;
    this.lastFailureTime = null;
    this.nextAttempt = null;
    this.halfOpenRequests = 0;
    
    // Metrics
    this.metrics = {
      requests: 0,
      successful: 0,
      failed: 0,
      rejected: 0,
      stateChanges: []
    };
    
    // Success rate tracking (sliding window)
    this.window = [];
    this.windowSize = 100;
    this.windowDuration = 60000; // 1 minute
  }
  
  async execute(fn, fallbackFn = null) {
    this.metrics.requests++;
    
    // Check if we should attempt
    if (!this.canAttempt()) {
      this.metrics.rejected++;
      
      if (fallbackFn) {
        return await fallbackFn();
      }
      
      throw new Error(`Circuit breaker is ${this.state}. Service temporarily unavailable.`);
    }
    
    try {
      // Execute the function
      const result = await fn();
      
      this.recordSuccess();
      return result;
      
    } catch (error) {
      this.recordFailure();
      
      // Try fallback if available
      if (fallbackFn && this.state === this.states.OPEN) {
        return await fallbackFn();
      }
      
      throw error;
    }
  }
  
  canAttempt() {
    switch (this.state) {
      case this.states.CLOSED:
        return true;
        
      case this.states.OPEN:
        // Check if timeout has passed
        if (Date.now() >= this.nextAttempt) {
          this.transitionTo(this.states.HALF_OPEN);
          return true;
        }
        return false;
        
      case this.states.HALF_OPEN:
        // Allow limited requests in half-open state
        return this.halfOpenRequests < this.config.halfOpenMaxRequests;
        
      default:
        return false;
    }
  }
  
  recordSuccess() {
    this.metrics.successful++;
    this.updateWindow(true);
    
    switch (this.state) {
      case this.states.CLOSED:
        this.failures = 0;
        break;
        
      case this.states.HALF_OPEN:
        this.successes++;
        this.halfOpenRequests++;
        
        if (this.successes >= this.config.successThreshold) {
          this.transitionTo(this.states.CLOSED);
        }
        break;
        
      case this.states.OPEN:
        // Shouldn't happen, but handle gracefully
        this.transitionTo(this.states.HALF_OPEN);
        break;
    }
  }
  
  recordFailure() {
    this.metrics.failed++;
    this.updateWindow(false);
    this.lastFailureTime = Date.now();
    
    switch (this.state) {
      case this.states.CLOSED:
        this.failures++;
        
        // Check if we should open
        const recentRequests = this.getRecentRequests();
        if (recentRequests.length >= this.config.volumeThreshold) {
          const failureRate = this.calculateFailureRate();
          
          if (this.failures >= this.config.failureThreshold || failureRate > 0.5) {
            this.transitionTo(this.states.OPEN);
          }
        }
        break;
        
      case this.states.HALF_OPEN:
        // Failure in half-open state immediately opens circuit
        this.transitionTo(this.states.OPEN);
        break;
        
      case this.states.OPEN:
        // Reset timeout on new failures
        this.nextAttempt = Date.now() + this.config.timeout;
        break;
    }
  }
  
  transitionTo(newState) {
    const oldState = this.state;
    this.state = newState;
    
    this.metrics.stateChanges.push({
      from: oldState,
      to: newState,
      timestamp: Date.now(),
      failureRate: this.calculateFailureRate()
    });
    
    switch (newState) {
      case this.states.CLOSED:
        this.failures = 0;
        this.successes = 0;
        this.nextAttempt = null;
        break;
        
      case this.states.OPEN:
        this.nextAttempt = Date.now() + this.config.timeout;
        this.successes = 0;
        break;
        
      case this.states.HALF_OPEN:
        this.failures = 0;
        this.successes = 0;
        this.halfOpenRequests = 0;
        break;
    }
  }
  
  updateWindow(success) {
    const now = Date.now();
    
    // Add new entry
    this.window.push({
      timestamp: now,
      success
    });
    
    // Remove old entries
    const cutoff = now - this.windowDuration;
    this.window = this.window.filter(entry => entry.timestamp > cutoff);
    
    // Keep window size limited
    if (this.window.length > this.windowSize) {
      this.window = this.window.slice(-this.windowSize);
    }
  }
  
  getRecentRequests() {
    const now = Date.now();
    const cutoff = now - this.windowDuration;
    return this.window.filter(entry => entry.timestamp > cutoff);
  }
  
  calculateFailureRate() {
    const recent = this.getRecentRequests();
    if (recent.length === 0) return 0;
    
    const failures = recent.filter(entry => !entry.success).length;
    return failures / recent.length;
  }
  
  getMetrics() {
    return {
      state: this.state,
      ...this.metrics,
      failureRate: (this.calculateFailureRate() * 100).toFixed(1) + '%',
      successRate: ((1 - this.calculateFailureRate()) * 100).toFixed(1) + '%',
      recentRequests: this.getRecentRequests().length,
      nextAttempt: this.nextAttempt ? new Date(this.nextAttempt).toISOString() : null,
      timeUntilNextAttempt: this.nextAttempt ? Math.max(0, this.nextAttempt - Date.now()) : null
    };
  }
  
  reset() {
    this.transitionTo(this.states.CLOSED);
    this.window = [];
    this.metrics = {
      requests: 0,
      successful: 0,
      failed: 0,
      rejected: 0,
      stateChanges: []
    };
  }
}

module.exports = GradualCircuitBreaker;