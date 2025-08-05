// smart-router.js - Intelligent request routing based on complexity analysis
const { spawn } = require('child_process');
const math = require('mathjs');

class ComplexityAnalyzer {
  analyze(code) {
    // Score components (0-1 scale)
    const scores = {
      imports: this.analyzeImports(code),
      loops: this.analyzeLoops(code),
      dataSize: this.analyzeDataSize(code),
      operations: this.analyzeOperations(code),
      io: this.analyzeIO(code)
    };
    
    // Weighted average
    const weights = {
      imports: 0.3,
      loops: 0.2,
      dataSize: 0.2,
      operations: 0.2,
      io: 0.1
    };
    
    let totalScore = 0;
    for (const [key, weight] of Object.entries(weights)) {
      totalScore += scores[key] * weight;
    }
    
    return {
      score: totalScore,
      details: scores,
      recommendation: this.getRecommendation(totalScore, scores)
    };
  }
  
  analyzeImports(code) {
    const heavyLibraries = [
      'numpy', 'pandas', 'scipy', 'sklearn', 'scikit-learn',
      'tensorflow', 'torch', 'keras', 'matplotlib', 'seaborn',
      'nltk', 'opencv', 'pillow', 'rdkit', 'biopython'
    ];
    
    const imports = code.match(/(?:import|from)\s+(\w+)/g) || [];
    const heavyImports = imports.filter(imp => 
      heavyLibraries.some(lib => imp.includes(lib))
    );
    
    return Math.min(1, heavyImports.length / 2);
  }
  
  analyzeLoops(code) {
    const loops = code.match(/(?:for|while)\s+/g) || [];
    const nestedLoops = code.match(/(?:for|while)[\s\S]*?(?:for|while)/g) || [];
    const largeRanges = code.match(/range\s*\(\s*(\d+)/g) || [];
    
    let score = loops.length * 0.1;
    score += nestedLoops.length * 0.3;
    
    largeRanges.forEach(match => {
      const num = parseInt(match.match(/\d+/)[0]);
      if (num > 10000) score += 0.3;
      else if (num > 1000) score += 0.2;
      else if (num > 100) score += 0.1;
    });
    
    return Math.min(1, score);
  }
  
  analyzeDataSize(code) {
    // Check for large data structures
    const arrays = code.match(/\[[\s\S]*?\]/g) || [];
    let score = 0;
    
    arrays.forEach(arr => {
      const elements = arr.split(',').length;
      if (elements > 100) score += 0.3;
      else if (elements > 50) score += 0.2;
      else if (elements > 20) score += 0.1;
    });
    
    // Check for data generation
    if (/np\.random|np\.arange|np\.linspace/i.test(code)) score += 0.2;
    if (/DataFrame|Series/i.test(code)) score += 0.2;
    
    return Math.min(1, score);
  }
  
  analyzeOperations(code) {
    let score = 0;
    
    // Mathematical operations
    if (/\*\*|pow|exp|log|sqrt/i.test(code)) score += 0.1;
    if (/fft|fourier|transform/i.test(code)) score += 0.3;
    if (/solve|integrate|differentiate/i.test(code)) score += 0.3;
    
    // ML operations
    if (/fit|predict|train|model/i.test(code)) score += 0.4;
    if (/neural|network|layer/i.test(code)) score += 0.5;
    
    // Statistical operations
    if (/mean|std|var|correlation|anova/i.test(code)) score += 0.2;
    
    return Math.min(1, score);
  }
  
  analyzeIO(code) {
    let score = 0;
    
    // File operations
    if (/open\s*\(|read|write|save/i.test(code)) score += 0.3;
    
    // Network operations
    if (/requests|urllib|http|api/i.test(code)) score += 0.5;
    
    // Plotting
    if (/plot|scatter|hist|show|savefig/i.test(code)) score += 0.2;
    
    return Math.min(1, score);
  }
  
  getRecommendation(score, details) {
    if (score < 0.2) return 'local-math';
    if (score < 0.4 && details.imports < 0.3) return 'local-safe';
    return 'e2b';
  }
}

class LocalMathExecutor {
  async execute(code) {
    try {
      // Parse and evaluate simple math expressions
      const sanitized = this.sanitizeCode(code);
      const result = math.evaluate(sanitized);
      
      return {
        success: true,
        stdout: String(result),
        stderr: '',
        exitCode: 0,
        executionTime: 0,
        executor: 'local-math'
      };
    } catch (error) {
      return {
        success: false,
        stdout: '',
        stderr: error.message,
        exitCode: 1,
        executionTime: 0,
        executor: 'local-math'
      };
    }
  }
  
  sanitizeCode(code) {
    // Remove print statements and extract expressions
    const lines = code.split('\n')
      .filter(line => !line.trim().startsWith('#'))
      .map(line => line.replace(/print\s*\((.*)\)/, '$1'))
      .filter(line => line.trim());
    
    return lines.join(';');
  }
}

class LocalSafeExecutor {
  async execute(code) {
    // Execute simple Python code locally with restrictions
    const safeCode = `
import sys
import math

# Restricted environment
__builtins__ = {
    'print': print,
    'len': len,
    'range': range,
    'int': int,
    'float': float,
    'str': str,
    'round': round,
    'abs': abs,
    'min': min,
    'max': max,
    'sum': sum,
}

${code}
`;
    
    return new Promise((resolve) => {
      const python = spawn('python', ['-c', safeCode]);
      let stdout = '';
      let stderr = '';
      const startTime = Date.now();
      
      python.stdout.on('data', (data) => {
        stdout += data.toString();
      });
      
      python.stderr.on('data', (data) => {
        stderr += data.toString();
      });
      
      python.on('close', (exitCode) => {
        resolve({
          success: exitCode === 0,
          stdout,
          stderr,
          exitCode,
          executionTime: Date.now() - startTime,
          executor: 'local-safe'
        });
      });
      
      // Kill if takes too long
      setTimeout(() => {
        python.kill();
      }, 5000);
    });
  }
}

class SmartRouter {
  constructor(e2bManager) {
    this.e2bManager = e2bManager;
    this.localMath = new LocalMathExecutor();
    this.localSafe = new LocalSafeExecutor();
    this.analyzer = new ComplexityAnalyzer();
    
    // Metrics
    this.routingStats = {
      'local-math': 0,
      'local-safe': 0,
      'e2b': 0
    };
  }
  
  async route(code, options = {}) {
    const analysis = this.analyzer.analyze(code);
    const { forceE2B = false } = options;
    
    // Get current E2B load
    const e2bLoad = this.e2bManager.getMetrics().load || 0;
    
    // Dynamic routing decision
    let executor = analysis.recommendation;
    
    if (forceE2B) {
      executor = 'e2b';
    } else if (analysis.score < 0.3 && e2bLoad < 0.5) {
      // Use local for simple stuff unless E2B is idle
      executor = 'local-math';
    } else if (analysis.score < 0.5 && e2bLoad > 0.8) {
      // Offload medium complexity when E2B is busy
      executor = 'local-safe';
    }
    
    // Track routing decision
    this.routingStats[executor]++;
    
    // Execute based on routing
    let result;
    switch (executor) {
      case 'local-math':
        result = await this.localMath.execute(code);
        break;
      case 'local-safe':
        result = await this.localSafe.execute(code);
        break;
      default:
        result = await this.e2bManager.execute(code, options);
    }
    
    return {
      ...result,
      complexity: analysis,
      routing: executor
    };
  }
  
  getStats() {
    const total = Object.values(this.routingStats).reduce((a, b) => a + b, 0);
    return {
      total,
      routing: this.routingStats,
      e2bOffloadRate: total > 0 
        ? ((this.routingStats['local-math'] + this.routingStats['local-safe']) / total * 100).toFixed(1) + '%'
        : '0%'
    };
  }
}

module.exports = { SmartRouter, ComplexityAnalyzer, LocalMathExecutor, LocalSafeExecutor };