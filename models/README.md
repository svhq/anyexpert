# Model Testing System

This directory contains the model configuration and testing infrastructure for Ask Any Expert.

## Quick Start

### 1. List Available Models
```bash
npm run model:list
# or
node swap-model.js list
```

### 2. Check Current Model
```bash
npm run model:current
# or
node swap-model.js current
```

### 3. Swap to a Different Model
```bash
# Interactive mode
npm run model:swap

# Direct swap
node swap-model.js o4-mini-high
node swap-model.js claude-3.5-sonnet
node swap-model.js gpt-4o
```

### 4. Test a Model
```bash
# Test current model with quick profile (4 queries)
npm run model:test

# Test specific model
node tests/model-testing/model-test-runner.js o4-mini-high

# Run comprehensive test (12 queries)
npm run model:test-all o4-mini-high comprehensive
```

## Available Models

| Model Key | Name | Provider | Cost (per 1M tokens) | Context Window |
|-----------|------|----------|---------------------|----------------|
| glm-4.5-air | GLM-4.5 Air | Z-AI | FREE | 128K |
| o4-mini-high | OpenAI o4-mini High | OpenAI | $1.10/$4.40 | 200K |
| claude-3.5-sonnet | Claude 3.5 Sonnet | Anthropic | $3.00/$15.00 | 200K |
| gpt-4o | GPT-4 Optimized | OpenAI | $2.50/$10.00 | 128K |
| deepseek-r1 | DeepSeek R1 | DeepSeek | $0.14/$2.80 | 64K |

## Test Profiles

- **quick**: 4 essential queries across different capabilities
- **comprehensive**: 12 queries covering all features
- **tools-focus**: 6 queries requiring tools (math, search, code)
- **reasoning-focus**: 6 queries testing reasoning and expert knowledge

## Test Categories

1. **Basic**: Simple questions without tools
2. **Math**: Mathematical calculations requiring execution
3. **Search**: Current information requiring web search
4. **Code**: Programming tasks with code execution
5. **Reasoning**: Multi-step logic and problem solving
6. **Expert**: Domain-specific expert knowledge

## File Structure

```
models/
├── config.json       # Model configurations and pricing
├── model-loader.js   # Model loading and management utility
└── README.md        # This file

tests/model-testing/
├── test-queries.json      # Standardized test queries
├── model-test-runner.js   # Main test execution engine
├── compare-models.js      # Model comparison tool (TODO)
└── results/              # Test results by timestamp
```

## Usage Examples

### Test Multiple Models
```javascript
const { testModel } = require('./tests/model-testing/model-test-runner');

// Test different models with same queries
const models = ['glm-4.5-air', 'o4-mini-high', 'claude-3.5-sonnet'];

for (const model of models) {
  await testModel(model, 'quick');
}
```

### Calculate Costs
```javascript
const { ModelLoader } = require('./models/model-loader');
const loader = new ModelLoader();

// Calculate cost for 10,000 input + 2,000 output tokens
const cost = loader.calculateCost('o4-mini-high', 10000, 2000);
console.log(`Cost: $${cost.totalCost}`);
```

### Find Budget Models
```javascript
const budgetModels = loader.getModelsByPriceRange(1.0, 5.0);
console.log('Budget models:', budgetModels.map(m => m.key));
```

## Adding New Models

1. Edit `models/config.json`
2. Add model configuration with:
   - Unique key
   - OpenRouter model ID
   - Pricing information
   - Feature support flags
   - Context window limits

Example:
```json
"new-model": {
  "id": "provider/model-name",
  "name": "Display Name",
  "provider": "Provider",
  "features": {
    "tools": true,
    "search": true,
    "code": true,
    "streaming": true
  },
  "pricing": {
    "input": 1.00,
    "output": 2.00,
    "currency": "USD",
    "unit": "per_million_tokens"
  },
  "limits": {
    "contextWindow": 100000,
    "maxOutput": 4096,
    "rateLimit": "1000_rpm"
  }
}
```

## Tips

1. **Cost Management**: Always test with 'quick' profile first
2. **Rate Limits**: Add delays between tests for rate-limited models
3. **Validation**: Check test results in `results/` directory
4. **Comparison**: Use saved results to compare model performance

## Troubleshooting

- **Model not found**: Check model key in config.json
- **API errors**: Verify OpenRouter API key in .env
- **Rate limits**: Add delays or reduce test size
- **Cost concerns**: Use 'quick' profile or free models first