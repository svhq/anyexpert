// smart-cache.js - Intelligent caching with determinism assessment
const { LRUCache } = require('lru-cache');
const crypto = require('crypto');

class SmartCache {
  constructor(config = {}) {
    this.config = {
      maxSize: config.maxSize || 1000,
      ttl: config.ttl || 5 * 60 * 1000, // 5 minutes
      updateAgeOnGet: config.updateAgeOnGet || true,
      ...config
    };
    
    // LRU cache with TTL
    this.cache = new LRUCache({
      max: this.config.maxSize,
      ttl: this.config.ttl,
      updateAgeOnGet: this.config.updateAgeOnGet
    });
    
    // Track cache performance
    this.metrics = {
      hits: 0,
      misses: 0,
      evictions: 0,
      stores: 0
    };
    
    // Track determinism scores for adaptive caching
    this.determinismScores = new Map();
  }
  
  async getOrCompute(code, computeFn, options = {}) {
    const key = this.generateKey(code, options);
    
    // Check cache
    const cached = this.cache.get(key);
    
    if (cached) {
      this.metrics.hits++;
      
      // Check if result still has high confidence
      if (cached.confidence > 0.7) {
        return {
          ...cached.result,
          cached: true,
          cacheAge: Date.now() - cached.timestamp
        };
      }
    }
    
    this.metrics.misses++;
    
    // Compute result
    const startTime = Date.now();
    const result = await computeFn();
    const computeTime = Date.now() - startTime;
    
    // Assess determinism
    const determinism = this.assessDeterminism(code, result);
    
    // Update rolling determinism score
    this.updateDeterminismScore(key, determinism);
    
    // Cache if deterministic enough
    const avgDeterminism = this.determinismScores.get(key) || determinism;
    
    if (avgDeterminism > 0.5 && result.success) {
      this.metrics.stores++;
      
      this.cache.set(key, {
        result,
        confidence: avgDeterminism,
        timestamp: Date.now(),
        computeTime
      });
    }
    
    return {
      ...result,
      cached: false,
      computeTime
    };
  }
  
  generateKey(code, options = {}) {
    // Create a unique key based on code and relevant options
    const normalized = this.normalizeCode(code);
    const optionsStr = JSON.stringify({
      language: options.language || 'python',
      // Don't include transient options like timeoutMs
    });
    
    const hash = crypto
      .createHash('sha256')
      .update(normalized + optionsStr)
      .digest('hex');
    
    return hash.substring(0, 16); // Use first 16 chars
  }
  
  normalizeCode(code) {
    // Normalize code to improve cache hits
    return code
      .replace(/\s+/g, ' ') // Normalize whitespace
      .replace(/^\s+|\s+$/gm, '') // Trim lines
      .replace(/\n{2,}/g, '\n') // Remove multiple newlines
      .toLowerCase(); // Case insensitive
  }
  
  assessDeterminism(code, result) {
    let score = 1.0;
    
    // Check for non-deterministic patterns
    const patterns = [
      { regex: /random|rand|shuffle|sample/i, penalty: 0.5 },
      { regex: /datetime|time\.|now\(/i, penalty: 0.3 },
      { regex: /uuid|guid/i, penalty: 0.4 },
      { regex: /requests|urllib|fetch/i, penalty: 0.5 },
      { regex: /input\s*\(/i, penalty: 0.6 },
      { regex: /open\s*\(.*['"](r|rb|w|wb|a)/i, penalty: 0.4 }
    ];
    
    for (const { regex, penalty } of patterns) {
      if (regex.test(code)) {
        score *= (1 - penalty);
      }
    }
    
    // Check output variability (if we've seen this before)
    const key = this.generateKey(code);
    const cached = this.cache.get(key);
    
    if (cached && result.stdout) {
      // Compare outputs
      const similarity = this.calculateSimilarity(
        cached.result.stdout,
        result.stdout
      );
      
      score *= similarity;
    }
    
    return Math.max(0, Math.min(1, score));
  }
  
  calculateSimilarity(str1, str2) {
    if (str1 === str2) return 1.0;
    if (!str1 || !str2) return 0.0;
    
    // Simple similarity based on common characters
    const set1 = new Set(str1.split(''));
    const set2 = new Set(str2.split(''));
    
    const intersection = new Set([...set1].filter(x => set2.has(x)));
    const union = new Set([...set1, ...set2]);
    
    return intersection.size / union.size;
  }
  
  updateDeterminismScore(key, newScore) {
    const current = this.determinismScores.get(key);
    
    if (current === undefined) {
      this.determinismScores.set(key, newScore);
    } else {
      // Rolling average with more weight on recent scores
      const updated = current * 0.7 + newScore * 0.3;
      this.determinismScores.set(key, updated);
    }
    
    // Limit map size
    if (this.determinismScores.size > this.config.maxSize * 2) {
      // Remove oldest entries
      const toRemove = this.determinismScores.size - this.config.maxSize;
      const keys = Array.from(this.determinismScores.keys());
      
      for (let i = 0; i < toRemove; i++) {
        this.determinismScores.delete(keys[i]);
      }
    }
  }
  
  invalidate(pattern) {
    let invalidated = 0;
    
    if (pattern instanceof RegExp) {
      // Invalidate entries matching pattern
      for (const [key, value] of this.cache.entries()) {
        if (pattern.test(value.result.stdout || '')) {
          this.cache.delete(key);
          invalidated++;
        }
      }
    } else if (typeof pattern === 'string') {
      // Invalidate specific key
      if (this.cache.delete(pattern)) {
        invalidated++;
      }
    } else {
      // Clear all
      const size = this.cache.size;
      this.cache.clear();
      invalidated = size;
    }
    
    return invalidated;
  }
  
  getMetrics() {
    const total = this.metrics.hits + this.metrics.misses;
    
    return {
      ...this.metrics,
      hitRate: total > 0 
        ? ((this.metrics.hits / total) * 100).toFixed(1) + '%'
        : 'N/A',
      size: this.cache.size,
      maxSize: this.config.maxSize,
      avgDeterminismScore: this.determinismScores.size > 0
        ? (Array.from(this.determinismScores.values())
            .reduce((a, b) => a + b, 0) / this.determinismScores.size).toFixed(3)
        : 'N/A'
    };
  }
  
  clear() {
    this.cache.clear();
    this.determinismScores.clear();
    this.metrics = {
      hits: 0,
      misses: 0,
      evictions: 0,
      stores: 0
    };
  }
}

module.exports = SmartCache;