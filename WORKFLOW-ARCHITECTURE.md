# Modular System Workflow Architecture

## Overview

This document provides a comprehensive guide to understanding how queries flow through the Ask Any Expert modular system, detailing every file activation and decision point.

## System Architecture Diagram

```
┌─────────────────┐
│   Client/User   │
│ (demo-*.js)     │
└────────┬────────┘
         │ query
         ▼
┌─────────────────────────┐
│ workflow-engine-modular │ ← config.js
│ - Route to agent        │ ← tool-config.js
│ - Generate requestId    │
└────────┬────────────────┘
         │
         ▼
┌─────────────────────────┐
│ unified-agent-modular   │ ← system-prompt-modular.js
│ - Planning loop         │
│ - Tool execution        │
│ - Progress assessment   │
└──┬──────────────────┬──┘
   │                  │
   ▼                  ▼
┌─────────────┐  ┌──────────────┐
│   Tools     │  │  Synthesis   │
├─────────────┤  ├──────────────┤
│ • Search    │  │ • Info Synth │
│ • Scrape    │  │ • Final Ans  │
│ • Code      │  └──────────────┘
│ • Reason    │
└─────────────┘
```

## Detailed Query Flow

### 1. Client Entry Point

**File**: `demo-askanyexpert.js` (or any client application)

```javascript
const workflowEngine = require('./src/workflow-engine');
const response = await workflowEngine.answer(query, chatHistory);
```

**Purpose**: Initiates the query processing

### 2. Workflow Engine Router

**File**: `src/workflow-engine-modular.js`

**Key Functions**:
- `constructor()` - Determines modular vs original system
- `answer(userQuery, chatHistory)` - Main entry point

**Process**:
1. Checks `USE_MODULAR_SYSTEM` environment variable
2. If true, loads `unified-agent-modular`
3. Generates unique `requestId` for tracking
4. Calls `agent.process(userQuery, chatHistory, { requestId })`
5. Enhances response with backward compatibility fields
6. Returns to client

**Configuration Files Loaded**:
- `config.js` - Environment variables
- `config/tool-config.js` - API mode and tool definitions

### 3. Unified Agent Processing

**File**: `src/unified-agent-modular.js`

**Key Functions**:
- `process(userQuery, chatHistory, options)` - Main processing loop
- `planNextAction()` - Determines which tool to use
- `executeAction()` - Routes to specific tool handler
- `assessProgress()` - Evaluates confidence level

**Main Processing Loop**:
```javascript
while (stepNum < this.maxRounds && confidence < this.confidenceThreshold) {
    1. Plan next action
    2. Check tool availability
    3. Execute action(s)
    4. Store results
    5. Assess progress
}
```

### 4. Tool Selection and Execution

#### 4.1 Planning Phase

**Function**: `planNextAction(userQuery, steps, chatHistory)`

**Process**:
1. Builds planning prompt with available actions
2. Calls LLM to determine next action
3. Returns action object with type and rationale

**Decision Logic**:
- First query → Usually search or code based on content
- Has search results → May scrape for details
- Need calculations → Code execution
- Sufficient info → Synthesize final answer

#### 4.2 Tool Execution Paths

##### A. Search Path
```
planNextAction → "search"
↓
executeSearch()
↓
search-planner.js → generate queries
↓
web-search.js → runBatch()
↓
Format results with sources
```

**Files**:
- `src/search-planner.js` - Generates 3-5 search queries
- `src/web-search.js` - Executes via Serper API

##### B. Scrape Path
```
planNextAction → "scrape"
↓
executeScrape()
↓
web-search.js → scrape(url)
↓
Format content (truncate if >10K chars)
```

**Files**:
- `src/web-search.js` - scrape() method

##### C. Code Path
```
planNextAction → "code"
↓
executeCode()
↓
Choose: executeCodeViaTools or executeCodeViaExtraction
↓
e2b-manager.js → executeCode()
↓
e2b-orchestrator → pool management
↓
Return execution results
```

**Files**:
- `src/e2b-manager.js` - Singleton manager
- `src/e2b-orchestrator/index.js` - Advanced orchestration
- `src/e2b-orchestrator/pool-manager.js` - Sandbox pooling
- `src/e2b-orchestrator/sandbox-executor.js` - Actual execution

##### D. Reason Path
```
planNextAction → "reason"
↓
executeReason()
↓
Build reasoning prompt with context
↓
Call LLM with expert persona
```

**Files**:
- `src/system-prompt-modular.js` - Expert persona system

### 5. Progress Assessment

**Function**: `assessProgress(userQuery, steps)`

**Logic**:
- 0 steps → 0% confidence
- Has search results → 70-80% confidence
- Has code execution → 80-90% confidence
- Multiple steps → Higher confidence
- Synthesized answer → 95%+ confidence

### 6. Response Synthesis

**When**: confidence >= threshold OR max rounds reached

**Process**:
```
synthesizeFinalAnswer()
↓
Build synthesis prompt with all steps
↓
Call LLM to create comprehensive answer
↓
Return with metadata
```

**For Search Results**:
- Uses `information-synthesizer.js` → `compose()`
- Formats sources with citations

### 7. Response Enhancement

**File**: `src/workflow-engine-modular.js`

**Adds**:
- `answer` field (backward compatibility)
- `searchPerformed` boolean
- API mode information
- Timing metrics

## Example Query Flows

### Example 1: Search-Only Query

**Query**: "What are the latest features of Claude 3.5?"

```
1. workflow-engine → agent.process()
2. planNextAction → "search" (need current info)
3. search-planner → generates 5 queries
4. web-search → executes searches
5. assessProgress → 80% confidence
6. planNextAction → "synthesize"
7. synthesizeFinalAnswer → comprehensive response
8. Return with sources
```

### Example 2: Code Execution Query

**Query**: "Calculate the 30th Fibonacci number"

```
1. workflow-engine → agent.process()
2. planNextAction → "code" (calculation needed)
3. executeCode → tool calling method
4. e2b-manager → get sandbox from pool
5. Execute Python code
6. assessProgress → 90% confidence
7. synthesizeFinalAnswer → result with explanation
8. Return with execution details
```

### Example 3: Mixed Tool Query

**Query**: "Scrape https://example.com and count the words"

```
1. workflow-engine → agent.process()
2. planNextAction → "scrape" (URL provided)
3. executeScrape → fetch content
4. assessProgress → 80% confidence
5. planNextAction → "code" (need word count)
6. executeCode → count words in Python
7. assessProgress → 95% confidence
8. synthesizeFinalAnswer → complete answer
9. Return with both scrape and code results
```

## Key Decision Points

### 1. Tool Availability Check
```javascript
if (!this.hasToolAvailable(toolMapping[action.type])) {
    // Fall back to reasoning
    action.type = 'reason';
}
```

### 2. Parallel Execution
```javascript
if (action.parallelActions && action.parallelActions.length > 1) {
    // Execute multiple tools in parallel
    results = await Promise.all(promises);
}
```

### 3. Confidence Threshold
```javascript
if (confidence >= this.confidenceThreshold) {
    // Stop processing and synthesize
    break;
}
```

## Utility Functions

### Logging
- Every major step logged via `src/utils/logger.js`
- Includes requestId for tracking
- Metrics for performance monitoring

### JSON Parsing
- `src/utils/json-parser.js` handles various LLM output formats
- Extracts JSON from markdown blocks
- Fallback parsing strategies

### API Communication
- All LLM calls go through `callModel()`
- Configurable temperature and token limits
- Error handling and retries

## Environment Variables

Critical configuration that affects workflow:

```bash
USE_MODULAR_SYSTEM=true    # Use modular vs original
API_MODE=full              # none/search/full
MAX_ROUNDS=3               # Max planning iterations
CONF_THRESHOLD=0.8         # Confidence to stop
```

## Performance Optimizations

1. **Sandbox Pooling**: E2B sandboxes are reused
2. **Parallel Search**: Multiple queries executed together
3. **Early Stopping**: Confidence threshold prevents unnecessary steps
4. **Result Caching**: Search results stored in steps

## Error Handling

- Tool failures fall back to reasoning
- Scraping errors return null (handled gracefully)
- E2B failures have retry logic
- All errors logged with context

This architecture ensures robust, flexible query processing with intelligent tool selection and high-quality responses.