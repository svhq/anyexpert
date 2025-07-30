# Ask Any Expert - Web-Augmented Expert System

A cost-efficient, low-latency web-augmented reasoning system that combines an ensemble of AI experts with intelligent web search to deliver accurate, sourced answers.

## 🚀 Features

- **Expert Ensemble**: Automatically selects the most appropriate expert persona for each query
- **Intelligent Search**: Determines when web search is needed and performs up to 3 rounds of refined searches
- **Confidence-Based Stopping**: Stops searching when sufficient confidence is achieved
- **Citation Support**: Provides numbered citations for all web-sourced information
- **Structured Logging**: Complete request tracking and performance metrics
- **Cost Optimization**: Minimizes API calls while maximizing answer quality

## 📋 Prerequisites

- Node.js 20 LTS or higher
- OpenRouter API key
- Serper.dev API key

## 🛠️ Installation

1. Clone the repository:
```bash
git clone <repository-url>
cd askanyexpert
```

2. Install dependencies:
```bash
npm install
```

3. Create `.env` file with your API keys:
```env
OPENROUTER_API_KEY=your_openrouter_key
OPENROUTER_MODEL=z-ai/glm-4.5-air:free
SERPER_API_KEY=your_serper_key
```

## 🎯 Usage

### Command Line Interface
```bash
# Direct question
node index.js "What is the speed of light?"

# Question requiring web search
node index.js "What are the latest AI breakthroughs in 2025?"
```

### API Integration (Coming Soon)
```javascript
const workflowEngine = require('./src/workflow-engine');

const result = await workflowEngine.answer(
  "Your question here",
  { userId: 'user123', chatHistory: [...] }
);
```

## 🏗️ Architecture

```
┌──────────────────┐
│  workflow-engine │  ← Main orchestrator (3-round max)
└────────┬─────────┘
         ↓
┌────────────────────┐
│  query-analyzer    │  ← Determines if search needed
└────────────────────┘
         ↓ (if yes)
┌────────────────────┐
│  search-planner    │  ← Generates search queries
└────────────────────┘
         ↓
┌────────────────────┐
│  web-search        │  ← Executes Serper searches
└────────────────────┘
         ↓
┌────────────────────┐
│ info-synthesizer   │  ← Composes final response
└────────────────────┘
```

## 📊 Performance Targets

- **Latency**: ≤ 6 seconds mean response time
- **Cost**: ≤ 12 Serper API calls per complex question
- **Accuracy**: ≥ 60% exact match on benchmarks
- **Confidence**: 0.8 threshold for search completion

## 🧪 Testing

```bash
# Run all tests
npm test

# Run specific test suites
npm run test-api        # Basic API connectivity
npm run test-search     # Web search functionality
npm run test-workflow   # Full workflow integration

# Run with coverage
npm run test:coverage
```

## 📁 Project Structure

```
askanyexpert/
├── src/
│   ├── workflow-engine.js      # Main orchestrator
│   ├── query-analyzer.js       # Search necessity detection
│   ├── search-planner.js       # Query generation
│   ├── web-search.js          # Serper integration
│   ├── information-synthesizer.js # Response composition
│   ├── system-prompt.js       # Expert system prompt
│   └── utils/
│       └── logger.js          # Structured logging
├── tests/
│   └── workflow.test.js       # Jest test suite
├── .env                       # API keys (not in repo)
├── config.js                  # Configuration
└── index.js                   # CLI entry point
```

## 🔄 Development Phases

### ✅ Phase 1 (Complete)
- Core workflow engine with 3-round search
- Query analysis and search planning
- Basic web search integration
- Structured logging

### 🚧 Phase 2 (Next)
- Redis caching for search results
- Passage ranking with embeddings
- Source reliability evaluation
- Performance optimizations

### 📅 Phase 3 (Planned)
- Expert persona selection
- Context management
- Confidence scoring improvements
- Enhanced citation system

### 🔮 Phase 4 (Future)
- Contradiction detection
- Domain-specific templates
- Advanced monitoring
- UI/API server

## 🤝 Contributing

1. Follow the existing code style
2. Add tests for new features
3. Update documentation
4. Use structured logging for debugging

## 📝 License

MIT License - see LICENSE file for details