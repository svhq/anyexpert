/**
 * E2B Orchestrator Configuration
 * Centralized configuration for all orchestrator components
 */

const isTest = process.env.NODE_ENV === 'test' || process.argv.includes('test');
const isProduction = process.env.NODE_ENV === 'production';

module.exports = {
  // Pool Manager Configuration
  pool: {
    minSize: parseInt(process.env.E2B_POOL_MIN_SIZE) || 2,
    targetSize: parseInt(process.env.E2B_POOL_TARGET_SIZE) || 3,
    maxSize: parseInt(process.env.E2B_POOL_MAX_SIZE) || 10,
    sandboxTimeout: parseInt(process.env.E2B_SANDBOX_TIMEOUT) || 30 * 60 * 1000, // 30 minutes
    maxReuseCount: parseInt(process.env.E2B_MAX_REUSE_COUNT) || 50,
    maxAge: parseInt(process.env.E2B_MAX_AGE) || 60 * 60 * 1000, // 1 hour
    scaleUpThreshold: parseFloat(process.env.E2B_SCALE_UP_THRESHOLD) || 0.7,
    scaleDownThreshold: parseFloat(process.env.E2B_SCALE_DOWN_THRESHOLD) || 0.3,
    maintenanceInterval: isTest ? 0 : 30000 // Disable in tests
  },
  
  // Queue Configuration
  queue: {
    globalConcurrency: parseInt(process.env.E2B_GLOBAL_CONCURRENCY) || 6,
    perUserConcurrency: parseInt(process.env.E2B_PER_USER_CONCURRENCY) || 2,
    cleanupInterval: isTest ? 0 : 60000 // Disable in tests
  },
  
  // Circuit Breaker Configuration
  circuitBreaker: {
    failureThreshold: parseInt(process.env.E2B_FAILURE_THRESHOLD) || 5,
    timeout: parseInt(process.env.E2B_CIRCUIT_TIMEOUT) || 30000
  },
  
  // Cache Configuration
  cache: {
    enabled: process.env.E2B_CACHE_ENABLED !== 'false',
    maxSize: parseInt(process.env.E2B_CACHE_MAX_SIZE) || 1000,
    ttl: parseInt(process.env.E2B_CACHE_TTL) || 5 * 60 * 1000 // 5 minutes
  },
  
  // Rate Limiter Configuration
  rateLimiter: {
    hobby: process.env.E2B_TIER !== 'pro',
    refillInterval: isTest ? 0 : 1000 // Disable in tests
  },
  
  // General Configuration
  enableLocalExecution: false, // Always use E2B
  isTest,
  isProduction
};