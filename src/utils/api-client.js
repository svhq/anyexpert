const config = require('../../config');
const logger = require('./logger');

/**
 * Shared API client with robust error handling
 */
class APIClient {
  constructor() {
    this.config = config;
  }

  /**
   * Make a safe API call to OpenRouter with proper error handling
   * @param {Array} messages - Messages array for the API
   * @param {Object} options - Additional options (temperature, max_tokens, etc.)
   * @param {Object} metadata - Request metadata for logging
   * @returns {Promise<Object>} - { content: string, tokensUsed: number }
   */
  async callOpenRouter(messages, options = {}, metadata = {}) {
    const timeout = options.timeout || 60000; // 60 seconds default
    const controller = new AbortController();
    const timeoutId = setTimeout(() => controller.abort(), timeout);

    const requestBody = {
      model: this.config.openrouter.model,
      messages,
      temperature: options.temperature ?? 0.7,
      // Don't limit max_tokens to avoid truncation unless explicitly set
      ...(options.max_tokens && { max_tokens: options.max_tokens }),
      ...options
    };

    try {
      logger.info({
        type: 'api-request',
        model: requestBody.model,
        messageCount: messages.length,
        ...metadata
      });

      const response = await fetch(`${this.config.openrouter.baseUrl}/chat/completions`, {
        method: 'POST',
        headers: this.config.openrouter.headers,
        body: JSON.stringify(requestBody),
        signal: controller.signal
      });

      clearTimeout(timeoutId);

      // Handle non-OK responses
      if (!response.ok) {
        let errorMessage = `OpenRouter API error: ${response.status}`;
        try {
          const errorText = await response.text();
          errorMessage += ` - ${errorText}`;
        } catch (e) {
          // Ignore text parsing error
        }
        throw new Error(errorMessage);
      }

      // Read response as text first to handle potential JSON errors
      const responseText = await response.text();
      
      // Check if response is empty
      if (!responseText) {
        throw new Error('Empty response from OpenRouter API');
      }

      // Try to parse JSON
      let data;
      try {
        data = JSON.parse(responseText);
      } catch (parseError) {
        logger.error({
          type: 'json-parse-error',
          error: parseError.message,
          responseLength: responseText.length,
          responsePreview: responseText.substring(0, 200),
          ...metadata
        });
        
        // Check if it might be a partial response
        if (responseText.includes('"choices"') && responseText.includes('"content"')) {
          // Try to extract content with regex as last resort
          const contentMatch = responseText.match(/"content"\s*:\s*"([^"]+)"/);
          if (contentMatch) {
            logger.warn({
              type: 'partial-response-recovered',
              extractedLength: contentMatch[1].length,
              ...metadata
            });
            return {
              content: contentMatch[1],
              tokensUsed: 0,
              partial: true
            };
          }
        }
        
        throw new Error(`Invalid JSON response from API: ${parseError.message}`);
      }

      // Validate response structure
      if (!data.choices || !data.choices[0] || !data.choices[0].message) {
        logger.error({
          type: 'invalid-response-structure',
          hasChoices: !!data.choices,
          choicesLength: data.choices?.length || 0,
          ...metadata
        });
        throw new Error('Invalid response structure from OpenRouter API');
      }

      // Extract content
      const content = data.choices[0].message.content;
      const tokensUsed = data.usage?.total_tokens || 0;

      logger.info({
        type: 'api-response',
        contentLength: content.length,
        tokensUsed,
        ...metadata
      });

      return {
        content,
        tokensUsed
      };

    } catch (error) {
      clearTimeout(timeoutId);
      
      if (error.name === 'AbortError') {
        logger.error({
          type: 'api-timeout',
          timeout,
          ...metadata
        });
        throw new Error(`Request timeout after ${timeout/1000} seconds`);
      }
      
      logger.error({
        type: 'api-error',
        error: error.message,
        ...metadata
      });
      
      throw error;
    }
  }
}

module.exports = new APIClient();