// index.js - E2B Orchestrator that brings all components together
const { AdaptiveRateLimiter } = require('./token-bucket');
const { SmartRouter } = require('./smart-router');
const DynamicPoolManager = require('./pool-manager');
const IntelligentQueue = require('./intelligent-queue');
const GradualCircuitBreaker = require('./circuit-breaker');
const SmartCache = require('./smart-cache');
const logger = require('../utils/logger');
const orchestratorConfig = require('../../config/orchestrator-config');

class E2BOrchestrator {
  constructor(config = {}) {
    this.config = {
      hobby: config.hobby !== false, // Default to hobby tier
      enableCache: config.enableCache !== false,
      enableLocalExecution: false, // DISABLED - only use E2B
      ...config
    };
    
    // Initialize components
    this.rateLimiter = new AdaptiveRateLimiter({
      hobby: this.config.hobby
    });
    
    this.poolManager = new DynamicPoolManager(this.rateLimiter, orchestratorConfig.pool);
    
    this.queue = new IntelligentQueue(orchestratorConfig.queue);
    
    this.circuitBreaker = new GradualCircuitBreaker({
      failureThreshold: 5,
      timeout: 30000
    });
    
    this.cache = new SmartCache({
      maxSize: 1000,
      ttl: 5 * 60 * 1000
    });
    
    // Router needs a reference to this orchestrator
    this.router = new SmartRouter(this);
    
    // Metrics
    this.metrics = {
      totalRequests: 0,
      cachedRequests: 0,
      localRequests: 0,
      e2bRequests: 0,
      errors: 0
    };
    
    // Initialize pool
    this.initialized = false;
    this.initPromise = null;
  }
  
  async initialize() {
    if (this.initialized) return;
    if (this.initPromise) return this.initPromise;
    
    this.initPromise = this._initialize();
    await this.initPromise;
    this.initialized = true;
  }
  
  async _initialize() {
    logger.info({ message: 'Initializing E2B Orchestrator' });
    
    try {
      await this.poolManager.initialize();
      logger.info({ message: 'E2B Orchestrator initialized successfully' });
    } catch (error) {
      logger.error({ 
        message: 'Failed to initialize E2B Orchestrator',
        error: error.message
      });
      throw error;
    }
  }
  
  async execute(userId, code, options = {}) {
    await this.initialize();
    
    this.metrics.totalRequests++;
    const startTime = Date.now();
    
    try {
      // Check cache if enabled
      if (this.config.enableCache) {
        const cached = await this.cache.getOrCompute(
          code,
          async () => this._executeUncached(userId, code, options),
          options
        );
        
        if (cached.cached) {
          this.metrics.cachedRequests++;
          logger.debug({ 
            message: 'Cache hit',
            userId,
            cacheAge: cached.cacheAge
          });
        }
        
        return {
          ...cached,
          totalTime: Date.now() - startTime
        };
      }
      
      // Execute without cache
      return await this._executeUncached(userId, code, options);
      
    } catch (error) {
      this.metrics.errors++;
      logger.error({
        message: 'Orchestrator execution failed',
        userId,
        error: error.message
      });
      throw error;
    }
  }
  
  async _executeUncached(userId, code, options) {
    // Analyze complexity for routing
    const analysis = this.router.analyzer.analyze(code);
    
    // Route based on complexity and configuration
    if (this.config.enableLocalExecution && analysis.recommendation !== 'e2b') {
      const localResult = await this.router.route(code, options);
      
      if (localResult.routing === 'local-math' || localResult.routing === 'local-safe') {
        this.metrics.localRequests++;
        return localResult;
      }
    }
    
    // E2B execution path
    this.metrics.e2bRequests++;
    
    // Queue the request
    return await this.queue.enqueue(
      userId,
      async () => {
        // Check circuit breaker
        return await this.circuitBreaker.execute(
          async () => {
            // Get rate limit token for operations
            await this.rateLimiter.acquireToken('operations');
            
            // Get executor from pool
            const executor = await this.poolManager.getExecutor(userId);
            
            // Acquire per-sandbox lock
            const lock = await this.rateLimiter.acquireForSandbox(
              executor.sandboxId,
              'runCode'
            );
            
            try {
              // Execute code
              const result = await executor.execute(code, options);
              
              return {
                ...result,
                executor: 'e2b',
                complexity: analysis
              };
              
            } finally {
              // Release sandbox lock
              lock.release();
            }
          },
          // Fallback function if circuit is open
          this.config.enableLocalExecution 
            ? async () => {
                logger.warn({ 
                  message: 'Circuit breaker open, falling back to local execution',
                  userId
                });
                return this.router.route(code, { ...options, forceLocal: true });
              }
            : null
        );
      },
      { complexity: analysis.score }
    );
  }
  
  getMetrics() {
    const routerStats = this.router.getStats();
    const rateLimiterStatus = this.rateLimiter.getStatus();
    const poolMetrics = this.poolManager.getMetrics();
    const queueMetrics = this.queue.getMetrics();
    const circuitMetrics = this.circuitBreaker.getMetrics();
    const cacheMetrics = this.cache.getMetrics();
    
    return {
      orchestrator: {
        ...this.metrics,
        cacheHitRate: this.metrics.totalRequests > 0
          ? ((this.metrics.cachedRequests / this.metrics.totalRequests) * 100).toFixed(1) + '%'
          : 'N/A',
        localOffloadRate: this.metrics.totalRequests > 0
          ? ((this.metrics.localRequests / this.metrics.totalRequests) * 100).toFixed(1) + '%'
          : 'N/A',
        e2bRate: this.metrics.totalRequests > 0
          ? ((this.metrics.e2bRequests / this.metrics.totalRequests) * 100).toFixed(1) + '%'
          : 'N/A'
      },
      router: routerStats,
      rateLimiter: rateLimiterStatus,
      pool: poolMetrics,
      queue: queueMetrics,
      circuitBreaker: circuitMetrics,
      cache: cacheMetrics
    };
  }
  
  // Convenience method for E2B manager compatibility
  async executeCode(code, options = {}) {
    const userId = options.userId || 'default';
    return this.execute(userId, code, options);
  }
  
  async shutdown() {
    logger.info({ message: 'Shutting down E2B Orchestrator' });
    
    await this.poolManager.shutdown();
    await this.queue.shutdown();
    this.rateLimiter.destroy();
    
    logger.info({ 
      message: 'E2B Orchestrator shutdown complete',
      metrics: this.getMetrics()
    });
  }
}

// Export singleton instance
module.exports = new E2BOrchestrator();

// Graceful shutdown
process.on('SIGINT', async () => {
  await module.exports.shutdown();
  process.exit(0);
});

process.on('SIGTERM', async () => {
  await module.exports.shutdown();
  process.exit(0);
});