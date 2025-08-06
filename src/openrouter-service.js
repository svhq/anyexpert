const config = require('../config');
const logger = require('./utils/logger');

/**
 * OpenRouter Service - Handles direct API calls to OpenRouter
 */
class OpenRouterService {
  constructor() {
    this.config = config;
    this.model = config.MODEL || config.openrouter?.model || 'google/gemini-2.5-flash-lite';
    this.baseUrl = config.openrouter?.baseUrl || 'https://openrouter.ai/api/v1';
    this.headers = {
      'Authorization': `Bearer ${config.OPENROUTER_API_KEY}`,
      'Content-Type': 'application/json',
      'HTTP-Referer': 'https://askanyexpert.ai',
      'X-Title': 'AskAnyExpert GPQA Testing'
    };
  }

  /**
   * Send a request to OpenRouter API
   * @param {Array} messages - Array of message objects
   * @param {Object} options - Additional options (max_tokens, temperature, etc.)
   * @returns {Promise<Object>} - Response from OpenRouter
   */
  async sendRequest(messages, options = {}) {
    const requestBody = {
      model: this.model,
      messages: messages,
      max_tokens: options.max_tokens || 16000,
      temperature: options.temperature || 0.3,
      top_p: options.top_p || 0.9,
      frequency_penalty: options.frequency_penalty || 0,
      presence_penalty: options.presence_penalty || 0
    };

    logger.info('OpenRouter request', {
      model: this.model,
      messageCount: messages.length,
      maxTokens: requestBody.max_tokens
    });

    try {
      const response = await fetch(`${this.baseUrl}/chat/completions`, {
        method: 'POST',
        headers: this.headers,
        body: JSON.stringify(requestBody)
      });

      if (!response.ok) {
        const errorData = await response.json().catch(() => ({}));
        throw new Error(`OpenRouter API error: ${response.status} - ${JSON.stringify(errorData)}`);
      }

      const data = await response.json();
      
      if (!data.choices || data.choices.length === 0) {
        throw new Error('Empty response from OpenRouter');
      }

      const content = data.choices[0].message.content;
      const usage = data.usage;

      logger.info('OpenRouter response received', {
        model: this.model,
        tokensUsed: usage?.total_tokens || 0,
        responseLength: content.length
      });

      return {
        content,
        usage,
        model: data.model || this.model
      };

    } catch (error) {
      logger.error('OpenRouter API error', {
        error: error.message,
        model: this.model
      });
      throw error;
    }
  }

  /**
   * Set the model to use
   * @param {string} model - Model identifier
   */
  setModel(model) {
    this.model = model;
    logger.info(`OpenRouter model set to: ${model}`);
  }
}

module.exports = OpenRouterService;