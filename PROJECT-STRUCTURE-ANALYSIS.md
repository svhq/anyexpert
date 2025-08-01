# AskAnyExpert Project Structure Analysis
## Essential Files for Cloud Deployment

### 1. Core Application Files

#### Main Entry Points
- `server-gpqa-ui.js` - Main UI server
- `config.js` - Configuration loader
- `.env` - Environment variables (contains API keys)

#### API Routes
- `api/gpqa-test-api-streaming.js` - Main API with real-time streaming
- `api/gpqa-test-api.js` - Standard API endpoint

### 2. Workflow Engine & Agents

#### Core Engine
- `src/workflow-engine.js` - Main orchestrator
- `src/query-analyzer.js` - Routes queries to appropriate agents
- `src/system-prompt.js` - **CRITICAL: Contains all expert personas and system instructions**

#### Specialized Agents
- `src/math-agent.js` - Handles mathematical computations with E2B
- `src/code-agent-simple.js` - Handles coding questions
- `src/information-synthesizer.js` - Complex reasoning
- `src/web-search.js` - Web search integration
- `src/search-planner.js` - Search strategy planning

#### Services
- `src/openrouter-service.js` - OpenRouter API integration

#### Utilities
- `src/utils/logger.js` - Logging utilities
- `src/utils/api-client.js` - API client utilities

### 3. E2B Code Execution Service

#### Microservice (Runs separately on port 3001)
- `microservices/run-code/server-e2b-custom.js` - E2B sandbox service
- `microservices/run-code/package.json` - E2B service dependencies
- `microservices/run-code/.env` - Can inherit from parent

#### E2B Template Files (Optional - for custom template)
- `e2b-template/Dockerfile` - Custom E2B environment
- `e2b-template/e2b.toml` - E2B configuration
- `e2b-template/test_libraries.py` - Library verification

### 4. Frontend Files
- `gpqa-testing-ui.html` - Main UI interface
- `gpqa-10-questions-clean.json` - GPQA test questions

### 5. Configuration & Dependencies
- `package.json` - Main project dependencies
- `.env` - Environment variables:
  ```
  OPENROUTER_API_KEY=your-key
  OPENROUTER_MODEL=google/gemini-2.5-flash-lite
  E2B_API_KEY=your-e2b-key
  E2B_TEMPLATE_ID=optional-custom-template
  SERPER_API_KEY=optional-web-search
  ```

### 6. Data Files
- `gpqa-10-questions-clean.json` - Test questions
- `GPQA dataset/` - Original GPQA data (optional)

## Cloud Deployment Structure

### Option 1: Monolithic Deployment
Deploy everything as a single application:
```
askanyexpert/
├── server-gpqa-ui.js
├── config.js
├── package.json
├── .env
├── api/
│   └── gpqa-test-api-streaming.js
├── src/
│   ├── workflow-engine.js
│   ├── query-analyzer.js
│   ├── system-prompt.js
│   ├── math-agent.js
│   ├── code-agent-simple.js
│   ├── information-synthesizer.js
│   ├── web-search.js
│   ├── search-planner.js
│   ├── openrouter-service.js
│   └── utils/
│       ├── logger.js
│       └── api-client.js
├── microservices/run-code/
│   ├── server-e2b-custom.js
│   └── package.json
├── gpqa-testing-ui.html
└── gpqa-10-questions-clean.json
```

### Option 2: Microservices Deployment
Deploy as separate services:

1. **Main API Service**
   - All files from Option 1 except microservices/

2. **E2B Service** (separate container/instance)
   - `microservices/run-code/` directory
   - Runs on separate port (3001)

## Critical Files Summary

### Absolutely Essential:
1. **`src/system-prompt.js`** - Contains all expert personas and instructions
2. **`.env`** - API keys (OPENROUTER_API_KEY, E2B_API_KEY)
3. **`src/workflow-engine.js`** - Core orchestration logic
4. **`src/query-analyzer.js`** - Routing logic
5. **`src/math-agent.js`** - Math handling with E2B integration
6. **`config.js`** - Configuration loader

### E2B Service (can run separately):
- **`microservices/run-code/server-e2b-custom.js`** - Code execution service
- Uses E2B cloud sandboxes, no local Docker needed

### Frontend:
- **`gpqa-testing-ui.html`** - UI interface
- **`server-gpqa-ui.js`** - Express server

## Environment Variables Required:
```bash
# Required
OPENROUTER_API_KEY=sk-or-v1-xxx
E2B_API_KEY=e2b_xxx

# Optional but recommended
OPENROUTER_MODEL=google/gemini-2.5-flash-lite
E2B_TEMPLATE_ID=your-custom-template
SERPER_API_KEY=for-web-search

# Service URLs (can be configured for cloud)
CODE_EXECUTION_URL=http://localhost:3001
```

## Deployment Notes:
1. The E2B service can be deployed separately or together
2. System prompt is embedded in the code, not a separate file
3. All expert personas are defined in `src/system-prompt.js`
4. E2B uses cloud sandboxes - no local Docker required
5. Frontend can be served statically or through Express