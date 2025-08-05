# E2B Orchestrator Architecture Documentation

## Overview

The E2B Orchestrator is a sophisticated system for managing code execution in secure sandboxed environments. It provides intelligent request routing, dynamic pool management, rate limiting, caching, and fault tolerance for running Python code through E2B's cloud infrastructure.

## System Architecture

### High-Level Flow

```
User Query → Unified Agent → E2B Manager v3 → E2B Orchestrator → E2B Sandbox
                    ↓                              ↓
              Web Search Tool                 Local Cache
                                                   ↓
                                            Circuit Breaker
                                                   ↓
                                             Rate Limiter
                                                   ↓
                                            Request Queue
                                                   ↓
                                             Pool Manager
                                                   ↓
                                           Sandbox Executor
```

## Core Components

### 1. **Unified Agent** (`src/unified-agent.js`)
- Entry point for all user queries
- Implements multi-step reasoning with tool usage
- Manages conversation flow and confidence scoring
- Tools available:
  - `search_web`: Web search via Serper API
  - `run_code`: Code execution via E2B

### 2. **E2B Manager v3** (`src/e2b-manager-v3.js`)
- Simplified interface to the orchestrator
- Handles code execution requests
- Returns formatted results

### 3. **E2B Orchestrator** (`src/e2b-orchestrator/index.js`)
The main orchestration layer with these sub-components:

#### a. **Token Bucket Rate Limiter** (`token-bucket.js`)
- Implements E2B's rate limits:
  - Sandbox creation: 1/sec (Hobby) or 5/sec (Pro)
  - Operations: 40k/60s per IP
  - Lifecycle API: 20k/30s
- Adaptive health scoring based on failure rates

#### b. **Smart Router** (`smart-router.js`)
- Routes requests based on complexity analysis
- Currently disabled (all requests go to E2B)
- Prepared for future local execution options

#### c. **Dynamic Pool Manager** (`pool-manager.js`)
- Maintains a pool of warm sandboxes for fast execution
- Features:
  - Sandbox reuse (up to 50 executions per sandbox)
  - Health checks with BUILD_TAG verification
  - Auto-scaling based on load (2-10 sandboxes)
  - Graceful sandbox retirement after 1 hour

#### d. **Intelligent Queue** (`intelligent-queue.js`)
- Fair queuing with per-user concurrency limits
- Global concurrency: 6 requests
- Per-user concurrency: 2 requests
- Priority-based scheduling

#### e. **Circuit Breaker** (`circuit-breaker.js`)
- Gradual failure handling
- States: CLOSED → OPEN → HALF_OPEN
- Prevents cascade failures

#### f. **Smart Cache** (`smart-cache.js`)
- LRU cache for deterministic computations
- TTL: 5 minutes
- Max size: 1000 entries
- Determinism scoring for cache eligibility

#### g. **Sandbox Executor** (`sandbox-executor.js`)
- Manages individual sandbox execution
- Features:
  - Code wrapping with timeout protection
  - Error handling and output formatting
  - BUILD_TAG verification
  - Resource cleanup

## E2B Custom Template

### Template Configuration (`my-e2b-template/`)
- **Template Name**: `prod-all`
- **Template ID**: `izdx3u0apwnbdatk6pmh`
- **BUILD_TAG**: `npf-20250805-003130`

### Dockerfile (`e2b.Dockerfile`)
```dockerfile
FROM e2bdev/code-interpreter:latest

# Python science stack
RUN pip install --no-cache-dir \
    "numpy<2" numpy-financial scipy sympy mpmath pandas \
    matplotlib seaborn scikit-learn pillow networkx plotly \
    statsmodels yfinance biopython splicekit deeptools \
    beautifulsoup4 requests openpyxl xlrd python-dateutil

# Verify installations
RUN python -c "import numpy_financial; print('numpy_financial version:', numpy_financial.__version__)"
```

## Request Lifecycle

### 1. **Query Processing**
```javascript
// User query enters unified agent
const result = await unifiedAgent.process(query);
```

### 2. **Tool Selection**
- Agent analyzes query complexity
- Selects appropriate tools (search, code, or both)
- May execute tools in parallel for efficiency

### 3. **Code Execution Flow**
```javascript
// Agent calls run_code tool
{
  name: 'run_code',
  arguments: {
    code: 'import numpy_financial as npf\n...',
    timeout: 30000
  }
}
```

### 4. **Orchestrator Processing**
1. **Rate Limit Check**: Acquire tokens from appropriate bucket
2. **Cache Lookup**: Check if result exists for deterministic code
3. **Queue Management**: Enqueue request with fair scheduling
4. **Circuit Breaker**: Verify system health
5. **Pool Assignment**: Get available sandbox from pool
6. **Execution**: Run code in assigned sandbox
7. **Result Processing**: Format output and handle errors

### 5. **Sandbox Management**
- Sandboxes are pre-warmed for instant execution
- Each sandbox tracked with:
  - Execution count
  - Last used timestamp
  - Health status
  - BUILD_TAG version
- Automatic cleanup after max reuse or age

## Configuration

### Environment Variables (`.env`)
```env
E2B_API_KEY=your_api_key
E2B_TEMPLATE_ID=prod-all
```

### Orchestrator Config (`config/orchestrator-config.js`)
```javascript
{
  pool: {
    minSize: 2,
    targetSize: 3,
    maxSize: 10,
    maintenanceInterval: 30000  // Disabled in tests
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

## Monitoring & Metrics

The orchestrator provides comprehensive metrics:

```javascript
const metrics = orchestrator.getMetrics();
// Returns:
{
  orchestrator: { totalRequests, cacheHitRate, e2bRate },
  pool: { poolSize, created, reused, errors },
  queue: { enqueued, completed, successRate },
  circuitBreaker: { state, failureRate },
  cache: { hits, misses, hitRate },
  rateLimiter: { healthScore, availableTokens }
}
```

## Error Handling

### Retry Strategy
- Automatic retries for transient failures
- Exponential backoff with jitter
- Circuit breaker prevents repeated failures

### Graceful Degradation
- Falls back to new sandbox creation if pool exhausted
- Queues requests when at capacity
- Returns helpful error messages

## Testing

### Test Utilities (`test-utils.js`)
- Proper lifecycle management
- Automatic cleanup after tests
- Process handler isolation

### Example Test
```javascript
const { runTest } = require('./test-utils');

async function myTest() {
  const result = await unifiedAgent.process(query);
  // Assertions...
}

runTest('Test Name', myTest);
```

## Performance Optimizations

1. **Sandbox Pooling**: ~10x faster than creating new sandboxes
2. **Request Batching**: Parallel tool execution when possible
3. **Smart Caching**: Instant results for repeated calculations
4. **Connection Reuse**: Persistent sandbox connections

## Security

- All code runs in isolated E2B sandboxes
- No access to host filesystem
- Network isolation
- Resource limits enforced
- Automatic timeout protection

## Troubleshooting

### Common Issues

1. **Template Not Found**
   - Ensure `E2B_TEMPLATE_ID=prod-all` in `.env`
   - Verify template built successfully

2. **Import Errors**
   - Check BUILD_TAG matches expected version
   - Rebuild template if needed

3. **Timeout Errors**
   - Increase timeout in request
   - Check sandbox health

4. **Rate Limits**
   - Monitor token bucket status
   - Implement request throttling

## Future Enhancements

1. **Local Execution Router**: Enable safe local code execution
2. **Multi-Region Support**: Reduce latency with regional pools
3. **Advanced Caching**: Semantic similarity matching
4. **WebSocket Support**: Real-time execution updates