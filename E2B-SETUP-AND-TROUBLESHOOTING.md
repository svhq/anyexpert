# E2B Setup and Troubleshooting Guide

## Overview
This document captures all E2B-related changes, issues encountered, and solutions implemented during the integration of E2B sandbox into the Ask Any Expert system.

## Timeline of Changes

### Phase 1: Initial Microservice Architecture (Failed)
- **What we tried**: Created a separate microservice (`microservices/run-code`) that wrapped E2B
- **Issues**: 
  - Added unnecessary complexity
  - E2B service hanging (taking 161-257 seconds)
  - Token generation issues (140K tokens)
  - Service not starting reliably
- **Lesson**: Don't wrap E2B in another service layer - use SDK directly

### Phase 2: Direct SDK Integration (Current - Working)
- **What we implemented**: Direct E2B SDK usage with manager pattern
- **Key components**:
  - `src/e2b-manager.js` - Handles sandbox lifecycle, pooling, retries
  - Removed microservice layer entirely
  - Sandbox pooling for performance
  - Health checks and retry logic

## Common Issues and Solutions

### 1. Environment Variable Loading Issue
**Problem**: E2B API key not found, "401: Invalid API key" errors

**Root Cause**: When running scripts from parent directory, dotenv loads wrong .env file
```
Working directory: C:\Users\cc
Script directory: C:\Users\cc\askanyexpert
dotenv loads: C:\Users\cc\.env (wrong!)
Should load: C:\Users\cc\askanyexpert\.env
```

**Solution**: Explicitly specify .env path
```javascript
const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
```

### 2. Missing E2B Package
**Problem**: Module not found errors

**Solution**: Install E2B SDK
```bash
npm install @e2b/code-interpreter
```

### 3. API Key Parameter Issue
**Problem**: Passing `apiKey` to Sandbox.create() causes errors

**Solution**: E2B SDK reads from environment variable automatically
```javascript
// Wrong
const sandbox = await Sandbox.create({
  template: this.templateId,
  apiKey: this.apiKey,  // Don't do this!
  timeoutMs: this.sandboxTimeout
});

// Correct
const sandbox = await Sandbox.create({
  template: this.templateId,
  timeoutMs: this.sandboxTimeout
});
```

### 4. Sandbox Cleanup
**Problem**: Sandboxes not closing properly

**Solution**: Use proper cleanup method
```javascript
// Wrong
await sandbox.close();  // This method doesn't exist

// Correct
await sandbox.kill();
```

## Environment Setup

### Required Environment Variables
```env
# E2B Configuration
E2B_API_KEY=e2b_f80dbb11e8fd8f540b4bc6fce4a4b7aab22425c2
E2B_TEMPLATE_ID=prod-all  # Custom template with scientific libraries
```

### Custom Template
Our `prod-all` template includes:
- Python 3.x
- NumPy, SciPy, Pandas
- Scikit-learn
- Matplotlib
- SymPy
- RDKit (for chemistry)
- And more scientific computing libraries

## Architecture Overview

### Current Implementation
```
User Query
    ↓
Unified Agent (src/unified-agent.js)
    ↓
E2B Manager (src/e2b-manager.js)
    ↓
E2B SDK → E2B Cloud Sandbox
```

### Key Features
1. **Sandbox Pooling**: Pre-warmed sandboxes for instant execution
2. **Health Checks**: Verify sandbox is responsive before use
3. **Retry Logic**: Exponential backoff for transient failures
4. **Local Fallback**: Falls back to local Python if E2B fails
5. **Telemetry**: Track execution times and success rates

## Testing E2B Integration

### Test File Template
```javascript
const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const unifiedAgent = require('./src/unified-agent');
const e2bManager = require('./src/e2b-manager');

// Your test query here
const query = `...`;

async function test() {
  try {
    const result = await unifiedAgent.process(query, {}, {
      requestId: 'test-' + Date.now()
    });
    console.log('Result:', result);
  } catch (error) {
    console.error('Error:', error);
  }
  await e2bManager.shutdown();
}

test();
```

## Performance Metrics

From our testing:
- **Sandbox creation**: ~1-2 seconds (with warm pool: instant)
- **Code execution**: 1-3 seconds typically
- **Total overhead**: Minimal with pooling
- **Success rate**: 100% with proper setup

## Debugging Tips

1. **Check environment variables**:
   ```javascript
   console.log('E2B_API_KEY exists:', !!process.env.E2B_API_KEY);
   ```

2. **Enable verbose logging**:
   - E2B manager logs all operations
   - Check for warm sandbox creation
   - Monitor retry attempts

3. **Test direct SDK usage**:
   ```javascript
   const { Sandbox } = require('@e2b/code-interpreter');
   const sandbox = await Sandbox.create({ template: 'prod-all' });
   ```

## Lessons Learned

1. **Keep it simple**: Direct SDK usage > wrapping in microservices
2. **Environment matters**: Always verify .env is loading correctly
3. **Pool resources**: Pre-warm sandboxes eliminate startup delays
4. **Handle failures gracefully**: Implement retries and fallbacks
5. **Monitor performance**: Track metrics to identify issues early

## Future Improvements

1. Consider implementing sandbox reuse for multiple executions
2. Add more sophisticated error categorization
3. Implement sandbox lifecycle optimization
4. Add support for other languages beyond Python

## Successful E2B Query Tests

### Query 1: Mathematical Function Composition
- **Task**: Calculate f(g(x,y), g(f(x,y))) for complex functions
- **E2B Usage**: 2 executions using SymPy
- **Result**: Correctly calculated 2x⁴ + 2y⁴
- **Performance**: ~1.9 seconds average

### Query 2: Chemistry Calculation (RDKit)
- **Task**: Calculate molecular properties of Ibuprofen
- **E2B Usage**: 1 execution
- **Result**: MW = 206.285 g/mol
- **Performance**: ~1 second

### Query 3: Machine Learning (Scikit-learn)
- **Task**: Linear regression for house prices
- **E2B Usage**: 2 executions
- **Result**: Slope=100, Intercept=50000, R²=1.0
- **Performance**: ~2.9 seconds average

### Query 4: Data Visualization (Matplotlib)
- **Task**: Create scatter plot with trend line
- **E2B Usage**: 2 executions
- **Result**: Correlation=0.977, Optimal hours=6.25
- **Performance**: ~2.8 seconds average

### Query 5: Binary Search Tree Implementation
- **Task**: Implement BST with insert, delete, min/max, height operations
- **E2B Usage**: 0 executions (model provided implementation without running)
- **Result**: Complete BST implementation with all operations
- **Note**: Model sometimes provides code without execution

### Query 6: Code Optimization (Prime Numbers)
- **Task**: Analyze and optimize prime finding algorithm
- **E2B Usage**: 1 execution
- **Result**: Implemented Sieve of Eratosthenes, analyzed O(n²) vs O(n log log n)
- **Performance**: ~0.8 seconds

### Query 7: Sorting Algorithms Benchmark
- **Task**: Implement and benchmark 4 sorting algorithms
- **E2B Usage**: 1 execution (comprehensive benchmark)
- **Result**: Complete performance comparison table for different data sizes/types
- **Performance**: ~20.9 seconds (due to extensive benchmarking)
- **Key Finding**: Python's Timsort outperforms others, Quick Sort best for large random data

## Key Success Factors

1. **Proper Environment Setup**: Ensure .env is loaded from correct directory
2. **Sandbox Pooling**: Pre-warmed sandboxes eliminate cold starts
3. **Error Handling**: Graceful fallbacks prevent total failures
4. **Library Availability**: Custom template includes all scientific libraries
5. **Performance Varies**: Simple calculations ~1s, complex benchmarks can take 20+s