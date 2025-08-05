const unifiedAgent = require('./unified-agent');
const logger = require('./utils/logger');

// Configuration constants
const MAX_ROUNDS = parseInt(process.env.MAX_ROUNDS || '3', 10);
const CONF_THRESHOLD = parseFloat(process.env.CONF_THRESHOLD || '0.8');

/**
 * Main workflow engine that orchestrates the expert system with web search
 */
class WorkflowEngine {
  constructor(config) {
    this.config = config;
    this.maxRounds = MAX_ROUNDS;
    this.confidenceThreshold = CONF_THRESHOLD;
  }

  /**
   * Main entry point for answering user queries
   * @param {string} userQuery - The user's question
   * @param {Object} chatHistory - Previous conversation context
   * @returns {Promise<Object>} - The expert response with citations
   */
  async answer(userQuery, chatHistory = {}) {
    const startTime = Date.now();
    const requestId = this.generateRequestId();

    try {
      logger.info({ requestId, step: 'unified-agent', query: userQuery });
      
      // Use unified agent for all queries
      const response = await unifiedAgent.process(userQuery, chatHistory, { requestId });
      
      // Add backward compatibility fields
      const enhancedResponse = {
        ...response,
        answer: response.content, // For backward compatibility
        agent: 'unified',
        searchPerformed: response.metadata?.toolsUsed?.includes('search') || false
      };

      // Log final metrics for monitoring
      logger.info({ 
        type: 'request-complete',
        requestId,
        userId: chatHistory.userId || 'anonymous',
        totalDuration: Date.now() - startTime,
        ...response.metadata
      });

      return enhancedResponse;

    } catch (error) {
      logger.error({ 
        requestId, 
        error: error.message, 
        stack: error.stack,
        duration: Date.now() - startTime 
      });
      throw error;
    }
  }


  /**
   * Generate unique request ID for tracking
   */
  generateRequestId() {
    return `req-${Date.now()}-${Math.random().toString(36).substr(2, 9)}`;
  }
}

// Export singleton instance
module.exports = new WorkflowEngine(require('../config'));