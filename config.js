require('dotenv').config();
const { SYSTEM_PROMPT } = require('./src/system-prompt');

const config = {
  openrouter: {
    apiKey: process.env.OPENROUTER_API_KEY,
    model: process.env.OPENROUTER_MODEL || 'z-ai/glm-4.5-air:free',
    baseUrl: process.env.OPENROUTER_BASE_URL || 'https://openrouter.ai/api/v1',
    headers: {
      'Authorization': `Bearer ${process.env.OPENROUTER_API_KEY}`,
      'Content-Type': 'application/json',
      'HTTP-Referer': 'http://localhost:3000',
      'X-Title': 'Ask Any Expert'
    }
  },
  serper: {
    apiKey: process.env.SERPER_API_KEY,
    baseUrl: 'https://google.serper.dev'
  },
  codeExecution: {
    serviceUrl: process.env.CODE_EXECUTION_URL || 'http://localhost:3001',
    timeout: parseInt(process.env.CODE_EXECUTION_TIMEOUT || '30000', 10),
    maxIterations: parseInt(process.env.CODE_MAX_ITERATIONS || '2', 10)
  },
  systemPrompt: SYSTEM_PROMPT,
  // E2B specific settings
  E2B_API_KEY: process.env.E2B_API_KEY,
  E2B_TEMPLATE_ID: process.env.E2B_TEMPLATE_ID,
  E2B_TIMEOUT: parseInt(process.env.E2B_TIMEOUT || '120000', 10),
  // Model settings
  MODEL: process.env.OPENROUTER_MODEL || 'google/gemini-2.5-flash-lite',
  OPENROUTER_API_KEY: process.env.OPENROUTER_API_KEY
};

module.exports = config;