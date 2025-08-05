const logger = require('./utils/logger');
const { getCurrentConfig } = require('../config/tool-config');

// Configuration constants
const MAX_ROUNDS = parseInt(process.env.MAX_ROUNDS || '3', 10);
const CONF_THRESHOLD = parseFloat(process.env.CONF_THRESHOLD || '0.8');
const USE_MODULAR_SYSTEM = process.env.USE_MODULAR_SYSTEM === 'true';

/**
 * Modular Workflow Engine - Routes to appropriate agent based on configuration
 */
class ModularWorkflowEngine {
  constructor(config) {
    this.config = config;
    this.maxRounds = MAX_ROUNDS;
    this.confidenceThreshold = CONF_THRESHOLD;
    
    // Dynamically load the appropriate agent
    if (USE_MODULAR_SYSTEM) {
      this.agent = require('./unified-agent-modular');
      this.agentType = 'modular';
    } else {
      this.agent = require('./unified-agent');
      this.agentType = 'original';
    }
    
    // Get current tool configuration
    this.toolConfig = getCurrentConfig();
    
    logger.info({
      message: 'Modular Workflow Engine initialized',
      agentType: this.agentType,
      apiMode: this.toolConfig.mode,
      useModularSystem: USE_MODULAR_SYSTEM
    });
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
      logger.info({ 
        requestId, 
        step: `${this.agentType}-workflow-start`, 
        query: userQuery,
        apiMode: this.toolConfig.mode
      });
      
      // Use the configured agent (modular or original)
      const response = await this.agent.process(userQuery, chatHistory, { requestId });
      
      // Add backward compatibility fields and API mode info
      const enhancedResponse = {
        ...response,
        answer: response.content, // For backward compatibility
        agent: this.agentType,
        apiMode: this.toolConfig.mode,
        toolsAvailable: this.toolConfig.tools.map(t => t.function.name),
        capabilities: this.toolConfig.capabilities,
        searchPerformed: response.metadata?.toolsUsed?.includes('search') || 
                        response.metadata?.toolsUsed?.includes('search_web') || false
      };

      // Log final metrics for monitoring
      logger.info({ 
        type: 'modular-workflow-complete',
        requestId,
        userId: chatHistory.userId || 'anonymous',
        totalDuration: Date.now() - startTime,
        agentType: this.agentType,
        apiMode: this.toolConfig.mode,
        ...response.metadata
      });

      return enhancedResponse;

    } catch (error) {
      logger.error({ 
        requestId, 
        error: error.message, 
        stack: error.stack,
        agentType: this.agentType,
        apiMode: this.toolConfig.mode
      });
      
      throw new Error(`Workflow processing failed: ${error.message}`);
    }
  }

  /**
   * Get current configuration information
   * @returns {Object} Current configuration details
   */
  getConfiguration() {
    return {
      agentType: this.agentType,
      apiMode: this.toolConfig.mode,
      toolsAvailable: this.toolConfig.tools.map(t => t.function.name),
      capabilities: this.toolConfig.capabilities,
      useModularSystem: USE_MODULAR_SYSTEM,
      maxRounds: this.maxRounds,
      confidenceThreshold: this.confidenceThreshold
    };
  }

  /**
   * Switch between modular and original system (for testing)
   * @param {boolean} useModular - Whether to use modular system
   */
  switchSystem(useModular) {
    if (useModular !== USE_MODULAR_SYSTEM) {
      logger.info({
        message: 'Switching agent system',
        from: this.agentType,
        to: useModular ? 'modular' : 'original'
      });
      
      // Reload the appropriate agent
      delete require.cache[require.resolve('./unified-agent-modular')];
      delete require.cache[require.resolve('./unified-agent')];
      
      if (useModular) {
        this.agent = require('./unified-agent-modular');
        this.agentType = 'modular';
      } else {
        this.agent = require('./unified-agent');
        this.agentType = 'original';
      }
      
      // Update environment variable for future instances
      process.env.USE_MODULAR_SYSTEM = useModular.toString();
    }
  }

  /**
   * Validate that the current configuration is working correctly
   * @returns {Promise<Object>} Validation results
   */
  async validateConfiguration() {
    const config = this.getConfiguration();
    const issues = [];
    const warnings = [];

    // Check if tools match API mode
    if (config.apiMode === 'none' && config.toolsAvailable.length > 0) {
      issues.push('API mode is "none" but tools are still available');
    }

    if (config.apiMode === 'search' && config.toolsAvailable.includes('run_code')) {
      issues.push('API mode is "search" but code execution is available');
    }

    if (config.apiMode === 'search' && !config.toolsAvailable.includes('search_web')) {
      issues.push('API mode is "search" but web search is not available');
    }

    if (config.apiMode === 'full' && config.toolsAvailable.length < 2) {
      warnings.push('API mode is "full" but not all expected tools are available');
    }

    // Check capabilities consistency
    if (config.capabilities.webSearch && !config.toolsAvailable.includes('search_web')) {
      issues.push('Web search capability enabled but search_web tool not available');
    }

    if (config.capabilities.codeExecution && !config.toolsAvailable.includes('run_code')) {
      issues.push('Code execution capability enabled but run_code tool not available');
    }

    return {
      valid: issues.length === 0,
      issues,
      warnings,
      configuration: config
    };
  }

  generateRequestId() {
    return Math.random().toString(36).substring(2, 15);
  }
}

module.exports = ModularWorkflowEngine;