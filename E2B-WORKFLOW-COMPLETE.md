# E2B Workflow - Complete Documentation

## Overview

This document provides a complete explanation of how the E2B code execution workflow operates in the Ask Any Expert system, including all activated files and the step-by-step process flow.

## System Components & File Structure

### 1. Entry Points
- **`src/unified-agent.js`** - Main agent that processes user queries
- **`src/workflow-engine.js`** - Orchestrates the overall workflow
- **`index.js`** - Express server entry point

### 2. E2B Integration Layer
- **`src/e2b-manager-v3.js`** - Simplified interface to E2B orchestrator
- **`src/e2b-orchestrator/index.js`** - Main orchestration controller
- **`config/orchestrator-config.js`** - Configuration settings
- **`config/template-config.js`** - E2B template configuration

### 3. Orchestrator Components
- **`src/e2b-orchestrator/token-bucket.js`** - Rate limiting
- **`src/e2b-orchestrator/smart-router.js`** - Request routing (currently disabled)
- **`src/e2b-orchestrator/pool-manager.js`** - Sandbox pool management
- **`src/e2b-orchestrator/intelligent-queue.js`** - Fair request queuing
- **`src/e2b-orchestrator/circuit-breaker.js`** - Fault tolerance
- **`src/e2b-orchestrator/smart-cache.js`** - Result caching
- **`src/e2b-orchestrator/sandbox-executor.js`** - Code execution

### 4. E2B Template
- **`my-e2b-template/e2b.Dockerfile`** - Custom sandbox image
- **`my-e2b-template/e2b.toml`** - Template configuration

## Complete Workflow - Step by Step

### Step 1: User Query Reception
```
User → Express Server (index.js) → Workflow Engine
```

The query enters through the API endpoint and is passed to the workflow engine.

### Step 2: Agent Processing
```javascript
// src/unified-agent.js
async process(userQuery, chatHistory, options) {
  // 1. Initialize request
  const requestId = options.requestId || this.generateRequestId();
  
  // 2. Multi-step reasoning loop
  for (let stepNum = 0; stepNum < this.maxSteps; stepNum++) {
    // Plan next action
    const action = await this.planNextAction(userQuery, steps);
    
    // Execute action (search_web or run_code)
    const result = await this.executeAction(action, ...);
    
    // Check confidence
    if (confidence >= this.confidenceThreshold) break;
  }
}
```

### Step 3: Code Execution Request
When the agent decides to run code:

```javascript
// Tool call from agent
{
  name: 'run_code',
  arguments: {
    code: `
import numpy_financial as npf
# Calculate monthly payment
payment = npf.pmt(0.075/12, 12*15, 200000)
print(f"Monthly payment: $" + f"{-payment:.2f}")
    `,
    timeout: 30000
  }
}
```

### Step 4: E2B Manager Processing
```javascript
// src/e2b-manager-v3.js
async executeCode(code, options = {}) {
  const result = await this.orchestrator.execute({
    code,
    language: options.language || 'python',
    userId: options.userId || 'default',
    metadata: options.metadata || {}
  });
  
  return this.formatResult(result);
}
```

### Step 5: Orchestrator Flow

#### 5.1 Rate Limiting
```javascript
// src/e2b-orchestrator/token-bucket.js
await this.rateLimiter.acquireToken('operations', 1);
```

#### 5.2 Cache Check
```javascript
// src/e2b-orchestrator/smart-cache.js
const cacheKey = this.generateKey(code);
const cached = this.cache.get(cacheKey);
if (cached) return cached;
```

#### 5.3 Queue Management
```javascript
// src/e2b-orchestrator/intelligent-queue.js
await this.queue.enqueue(userId, async () => {
  // Execute in sandbox
});
```

#### 5.4 Circuit Breaker
```javascript
// src/e2b-orchestrator/circuit-breaker.js
return await this.circuitBreaker.execute(async () => {
  // Protected execution
});
```

#### 5.5 Pool Assignment
```javascript
// src/e2b-orchestrator/pool-manager.js
const executor = await this.poolManager.getExecutor();
// Assigns from pool or creates new sandbox
```

#### 5.6 Sandbox Execution
```javascript
// src/e2b-orchestrator/sandbox-executor.js
async execute(code, options) {
  // Wrap code with timeout protection
  const wrappedCode = this.wrapPythonCode(code, options.timeout);
  
  // Execute in E2B sandbox
  const result = await this.sandbox.runCode(wrappedCode, {
    language: 'python'
  });
  
  return this.formatOutput(result);
}
```

### Step 6: Result Processing
The result flows back through the stack:
1. Sandbox returns output
2. Cache stores result
3. Queue completes request
4. Orchestrator formats response
5. Manager returns to agent
6. Agent incorporates into response

## E2B Sandbox Details

### Template: `prod-all`
- **ID**: `izdx3u0apwnbdatk6pmh`
- **BUILD_TAG**: `npf-20250805-003130`

### Pre-installed Libraries
```python
# Scientific Computing
numpy, numpy_financial, scipy, sympy, mpmath

# Data Analysis
pandas, statistics, statsmodels

# Visualization
matplotlib, seaborn, plotly

# Machine Learning
scikit-learn

# Bioinformatics
biopython, splicekit, deeptools

# Web & Utilities
beautifulsoup4, requests, openpyxl, xlrd
```

### Sandbox Lifecycle
1. **Creation**: New sandboxes created with `prod-all` template
2. **Warm Pool**: 2-3 sandboxes kept ready
3. **Reuse**: Each sandbox handles up to 50 executions
4. **Health Checks**: BUILD_TAG verified on each use
5. **Retirement**: After 1 hour or 50 uses

## Performance Optimizations

### 1. Connection Pooling
- Maintains 2-3 warm sandboxes
- ~10x faster than creating new sandboxes
- Instant code execution

### 2. Smart Caching
- LRU cache for deterministic computations
- 5-minute TTL
- Instant results for repeated calculations

### 3. Request Batching
- Parallel tool execution when possible
- Fair queuing prevents user monopolization

### 4. Rate Limiting
- Respects E2B's limits
- Adaptive health scoring
- Graceful degradation

## Monitoring & Metrics

### Real-time Metrics
```javascript
{
  orchestrator: {
    totalRequests: 156,
    cacheHitRate: "23.5%",
    e2bRate: "100.0%"
  },
  pool: {
    poolSize: 3,
    created: 12,
    reused: 144,
    errors: 0
  },
  queue: {
    successRate: "98.7%",
    avgWaitTime: 125
  }
}
```

### Logging
All major events logged with structured data:
- Request start/complete
- Sandbox creation/assignment
- Execution success/failure
- Performance metrics

## Error Handling

### Retry Strategy
1. Transient failures: Automatic retry with backoff
2. Import errors: Adaptive behavior in prompts
3. Timeout protection: 30-second default
4. Circuit breaker: Prevents cascade failures

### Graceful Degradation
- Falls back to new sandbox if pool exhausted
- Queues requests when at capacity
- Returns helpful error messages

## Testing

### Test Files
- **`test-utils.js`** - Test lifecycle management
- **`test-loan-clean.js`** - Financial calculation tests
- **`test-workflow-demonstration.js`** - Full workflow demo

### Running Tests
```bash
# Single test with cleanup
node test-loan-clean.js

# Multiple finance calculations
node test-quick-finance.js

# Full workflow demonstration
node test-workflow-demonstration.js
```

## Configuration

### Environment Variables
```env
# E2B API Configuration
E2B_API_KEY=e2b_f80dbb11e8fd8f540b4bc6fce4a4b7aab22425c2
E2B_TEMPLATE_ID=prod-all

# OpenRouter Configuration
OPENROUTER_API_KEY=sk-or-v1-...
OPENROUTER_MODEL=google/gemini-2.5-flash-lite

# Serper API (for web search)
SERPER_API_KEY=7dd5ae9d21e0e339c6eaf29609d2b9abc3b658a4
```

### Orchestrator Settings
```javascript
// config/orchestrator-config.js
{
  pool: {
    minSize: 2,
    targetSize: 3,
    maxSize: 10
  },
  queue: {
    globalConcurrency: 6,
    perUserConcurrency: 2
  },
  cache: {
    enabled: true,
    maxSize: 1000,
    ttl: 300000  // 5 minutes
  }
}
```

## Troubleshooting

### Common Issues & Solutions

1. **numpy_financial Import Errors**
   - Solution: Rebuild template with BUILD_TAG
   - Verify: Check sandbox BUILD_TAG matches expected

2. **Timeout Errors in Tests**
   - Solution: Use test-utils.js for proper cleanup
   - Set NODE_ENV=test to disable maintenance intervals

3. **Template Not Found**
   - Verify E2B_TEMPLATE_ID=prod-all in .env
   - Check template exists: `npx e2b template list`

4. **Sandbox Creation Failures**
   - Check rate limits in metrics
   - Verify API key is valid
   - Monitor circuit breaker state

## Summary

The E2B workflow provides a robust, scalable solution for secure code execution with:
- ✅ Intelligent request routing and queuing
- ✅ Dynamic sandbox pool management
- ✅ Comprehensive error handling
- ✅ Performance optimizations
- ✅ Real-time monitoring
- ✅ Secure isolated execution

The system efficiently handles financial calculations, data analysis, and complex computations while maintaining security and performance.