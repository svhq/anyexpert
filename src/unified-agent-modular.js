const config = require('../config');
const logger = require('./utils/logger');
const { generateSystemPrompt } = require('./system-prompt-modular');
const { getCurrentConfig, generateLibraryDocumentation } = require('../config/tool-config');
const webSearch = require('./web-search');
const searchPlanner = require('./search-planner');
const e2bManager = require('./e2b-manager-v3');

/**
 * Modular Unified Agent - Dynamic tool loading based on API mode
 * Supports different API configurations: none, search-only, full
 */
class ModularUnifiedAgent {
  constructor() {
    this.config = config;
    this.maxSteps = 6;
    this.confidenceThreshold = 0.9;
    
    // Get current tool configuration
    this.toolConfig = getCurrentConfig();
    
    // Initialize with dynamic tool definitions
    this.toolDefinitions = this.toolConfig.tools;
    
    logger.info({
      message: 'Modular Unified Agent initialized',
      apiMode: this.toolConfig.mode,
      toolsAvailable: this.toolDefinitions.map(t => t.function.name),
      capabilities: this.toolConfig.capabilities
    });
  }

  /**
   * Check if a specific tool is available
   * @param {string} toolName - Name of the tool to check
   * @returns {boolean} Whether the tool is available
   */
  hasToolAvailable(toolName) {
    return this.toolDefinitions.some(tool => tool.function.name === toolName);
  }

  /**
   * Main entry point - process query with multi-step loop
   * @param {string} userQuery - The user's question
   * @param {Object} chatHistory - Previous conversation context
   * @param {Object} options - Processing options
   * @returns {Promise<Object>} - The final response
   */
  async process(userQuery, chatHistory = {}, options = {}) {
    const requestId = options.requestId || this.generateRequestId();
    const startTime = Date.now();
    
    logger.info({ 
      requestId, 
      step: 'modular-unified-agent-start', 
      query: userQuery,
      apiMode: this.toolConfig.mode,
      toolsAvailable: this.toolDefinitions.map(t => t.function.name)
    });

    const steps = [];
    let confidence = 0;
    let finalAnswer = null;

    try {
      // Multi-step processing loop
      for (let stepNum = 0; stepNum < this.maxSteps && confidence < this.confidenceThreshold; stepNum++) {
        logger.info({ 
          requestId, 
          stepNum, 
          previousConfidence: confidence,
          status: 'step-start' 
        });

        // Plan next action based on current state and available tools
        const action = await this.planNextAction(userQuery, steps, chatHistory);
        
        // Skip action if tool is not available in current mode
        if (action.type !== 'reason' && !this.hasToolAvailable(action.type === 'search' ? 'search_web' : 'run_code')) {
          logger.info({
            requestId,
            stepNum,
            message: `Tool ${action.type} not available in API mode ${this.toolConfig.mode}, switching to reasoning`
          });
          action.type = 'reason';
        }
        
        // Execute the action(s)
        let result;
        
        // Check for parallel execution based on model's suggestion
        if (action.parallelActions && action.parallelActions.length > 1) {
          // Filter out unavailable tools
          const availableActions = action.parallelActions.filter(actionType => {
            if (actionType === 'reason') return true;
            if (actionType === 'search') return this.hasToolAvailable('search_web');
            if (actionType === 'code') return this.hasToolAvailable('run_code');
            return false;
          });
          
          if (availableActions.length > 1) {
            logger.info({ 
              requestId, 
              stepNum,
              status: 'parallel-execution',
              actions: availableActions
            });
            
            // Execute available actions in parallel
            const parallelPromises = availableActions.map(actionType => {
              switch(actionType) {
                case 'search':
                  return this.executeSearch(userQuery, steps, requestId);
                case 'code':
                  return this.executeCode(userQuery, steps, chatHistory, requestId);
                case 'reason':
                  return this.executeReason(userQuery, steps, chatHistory);
                default:
                  throw new Error(`Unknown action type: ${actionType}`);
              }
            });
            
            const parallelResults = await Promise.all(parallelPromises);
            
            // Combine results
            result = {
              content: parallelResults.map((r, i) => 
                `${availableActions[i].toUpperCase()} Results:\n${r.content}`
              ).join('\n\n'),
              confidence: Math.max(...parallelResults.map(r => r.confidence || 0.5)),
              tokensUsed: parallelResults.reduce((sum, r) => sum + (r.tokensUsed || 0), 0),
              parallel: true
            };
            
            // Merge specific result properties
            parallelResults.forEach((r, i) => {
              if (r.sources) result.sources = r.sources;
              if (r.code) result.code = r.code;
              if (r.executionResults) result.executionResults = r.executionResults;
            });
            
            // Update action to reflect what was done
            action.type = 'parallel';
            action.executedActions = availableActions;
          } else {
            // Fall back to single action execution
            result = await this.executeAction(action, userQuery, steps, chatHistory, requestId);
          }
        } else {
          // Single action execution
          result = await this.executeAction(action, userQuery, steps, chatHistory, requestId);
        }
        
        // Store step information
        steps.push({
          stepNum,
          action,
          result,
          timestamp: Date.now()
        });

        // Assess progress and confidence
        confidence = await this.assessProgress(userQuery, steps);
        
        logger.info({ 
          requestId, 
          stepNum,
          action: action.type,
          newConfidence: confidence,
          status: 'step-complete'
        });

        // If we have a final answer with high confidence, break
        if (action.type === 'synthesize' && confidence >= this.confidenceThreshold) {
          finalAnswer = result;
          break;
        }
      }

      // If we haven't synthesized a final answer yet, do it now
      if (!finalAnswer) {
        finalAnswer = await this.synthesizeFinalAnswer(userQuery, steps, chatHistory, requestId);
      }

      // Log final metrics
      const finalMetrics = {
        requestId,
        totalDuration: Date.now() - startTime,
        totalSteps: steps.length,
        finalConfidence: confidence,
        toolsUsed: steps.map(s => s.action.type),
        tokensUsed: steps.reduce((sum, s) => sum + (s.result.tokensUsed || 0), 0),
        apiMode: this.toolConfig.mode
      };
      
      logger.info({ type: 'modular-unified-agent-complete', ...finalMetrics });

      return {
        content: finalAnswer.content,
        metadata: {
          ...finalMetrics,
          steps: steps,
          sources: this.extractSources(steps),
          finalAnswer
        }
      };

    } catch (error) {
      logger.error({ 
        requestId, 
        error: error.message, 
        stack: error.stack,
        apiMode: this.toolConfig.mode
      });
      
      throw new Error(`Processing failed: ${error.message}`);
    }
  }

  /**
   * Plan the next action based on current state and available tools
   */
  async planNextAction(userQuery, steps, chatHistory) {
    const messages = [
      {
        role: 'system',
        content: this.generateSystemMessage()
      },
      {
        role: 'user',
        content: this.buildPlanningPrompt(userQuery, steps, chatHistory)
      }
    ];

    const response = await this.callModel(messages, {
      max_tokens: 300,
      temperature: 0.3
    });

    return this.parseActionPlan(response.content);
  }

  /**
   * Generate system message with current tool configuration
   */
  generateSystemMessage() {
    const libraryDocs = this.toolConfig.capabilities.codeExecution 
      ? generateLibraryDocumentation(this.toolConfig.capabilities.pythonLibraries)
      : '';
    
    return generateSystemPrompt(this.toolConfig.systemPromptAddition, libraryDocs);
  }

  /**
   * Execute a specific action (with tool availability checks)
   */
  async executeAction(action, userQuery, steps, chatHistory, requestId) {
    switch (action.type) {
      case 'search':
        if (!this.hasToolAvailable('search_web')) {
          logger.warn({ requestId, message: 'Search requested but not available, falling back to reasoning' });
          return this.executeReason(userQuery, steps, chatHistory);
        }
        return this.executeSearch(userQuery, steps, requestId);
        
      case 'code':
        if (!this.hasToolAvailable('run_code')) {
          logger.warn({ requestId, message: 'Code execution requested but not available, falling back to reasoning' });
          return this.executeReason(userQuery, steps, chatHistory);
        }
        return this.executeCode(userQuery, steps, chatHistory, requestId);
        
      case 'reason':
        return this.executeReason(userQuery, steps, chatHistory);
        
      case 'synthesize':
        return this.synthesizeFinalAnswer(userQuery, steps, chatHistory, requestId);
        
      default:
        logger.warn({ requestId, message: `Unknown action type: ${action.type}, falling back to reasoning` });
        return this.executeReason(userQuery, steps, chatHistory);
    }
  }

  /**
   * Execute web search action
   */
  async executeSearch(userQuery, steps, requestId) {
    // Generate search query using search planner
    const searchQueries = await searchPlanner.generate(
      userQuery, 
      'Searching for relevant information',
      0,
      []
    );
    const searchQuery = Array.isArray(searchQueries) ? searchQueries[0] : searchQueries;
    
    logger.info({ 
      requestId, 
      type: 'web_search', 
      query: searchQuery 
    });

    const searchResults = await webSearch.search(searchQuery);
    
    return {
      content: this.formatSearchResults(searchResults),
      sources: searchResults.results || [],
      tokensUsed: 0,
      confidence: 0.7
    };
  }

  /**
   * Execute code execution action
   */
  async executeCode(userQuery, steps, chatHistory, requestId) {
    const messages = [
      {
        role: 'system',
        content: this.generateSystemMessage()
      },
      {
        role: 'user',
        content: this.buildCodePrompt(userQuery, steps, chatHistory)
      }
    ];

    const response = await this.callModel(messages, {
      tools: [this.toolDefinitions.find(t => t.function.name === 'run_code')],
      tool_choice: 'auto',
      max_tokens: 1000
    });

    if (response.tool_calls && response.tool_calls.length > 0) {
      const toolCall = response.tool_calls[0];
      const { code, timeout } = JSON.parse(toolCall.function.arguments);
      
      logger.info({ 
        requestId, 
        type: 'code_execution', 
        language: 'python',
        codeLength: code.length,
        success: true
      });

      const executionResult = await e2bManager.executeCode(code, { 
        timeout: timeout || 30000,
        userId: 'modular-agent'
      });
      
      return {
        content: response.content,
        code: code,
        executionResults: executionResult,
        tokensUsed: response.usage?.total_tokens || 0,
        confidence: 0.9
      };
    }

    return {
      content: response.content,
      tokensUsed: response.usage?.total_tokens || 0,
      confidence: 0.6
    };
  }

  /**
   * Execute reasoning-only action
   */
  async executeReason(userQuery, steps, chatHistory) {
    const messages = [
      {
        role: 'system',
        content: this.generateSystemMessage()
      },
      {
        role: 'user',
        content: this.buildReasoningPrompt(userQuery, steps, chatHistory)
      }
    ];

    const response = await this.callModel(messages, {
      max_tokens: 1500,
      temperature: 0.2
    });

    return {
      content: response.content,
      tokensUsed: response.usage?.total_tokens || 0,
      confidence: 0.8
    };
  }

  // ... [Rest of the methods remain the same as original unified-agent.js]
  // Including: synthesizeFinalAnswer, buildPlanningPrompt, buildCodePrompt, 
  // buildReasoningPrompt, parseActionPlan, assessProgress, callModel,
  // formatSearchResults, extractSources, generateRequestId

  /**
   * Synthesize final answer from all steps
   */
  async synthesizeFinalAnswer(userQuery, steps, chatHistory, requestId) {
    const messages = [
      {
        role: 'system',
        content: this.generateSystemMessage()
      },
      {
        role: 'user',
        content: this.buildSynthesisPrompt(userQuery, steps, chatHistory)
      }
    ];

    const response = await this.callModel(messages, {
      max_tokens: 2000,
      temperature: 0.1
    });

    return {
      content: response.content,
      tokensUsed: response.usage?.total_tokens || 0,
      confidence: 1.0
    };
  }

  buildPlanningPrompt(userQuery, steps, chatHistory) {
    let prompt = `Query: "${userQuery}"\n\n`;
    
    if (steps.length > 0) {
      prompt += `Previous steps taken:\n`;
      steps.forEach((step, i) => {
        prompt += `${i + 1}. ${step.action.type}: ${step.result.content.substring(0, 200)}...\n`;
      });
      prompt += `\n`;
    }

    prompt += `Available tools in current API mode (${this.toolConfig.mode}):\n`;
    this.toolDefinitions.forEach(tool => {
      prompt += `- ${tool.function.name}: ${tool.function.description}\n`;
    });

    prompt += `\nWhat should be the next action? Choose from: ${this.toolDefinitions.map(t => t.function.name.replace('_', '')).join(', ')}, reason, synthesize\n`;
    prompt += `Respond with just the action name and brief reasoning.`;

    return prompt;
  }

  buildCodePrompt(userQuery, steps, chatHistory) {
    let prompt = `As an expert, solve this query using Python code: "${userQuery}"\n\n`;
    
    if (steps.length > 0) {
      prompt += `Context from previous steps:\n`;
      steps.forEach((step, i) => {
        prompt += `${i + 1}. ${step.result.content.substring(0, 300)}...\n`;
      });
      prompt += `\n`;
    }

    prompt += `Write Python code to solve this problem. Use the run_code tool.`;
    return prompt;
  }

  buildReasoningPrompt(userQuery, steps, chatHistory) {
    let prompt = `As an expert, provide a comprehensive answer to: "${userQuery}"\n\n`;
    
    if (steps.length > 0) {
      prompt += `Context from previous analysis:\n`;
      steps.forEach((step, i) => {
        prompt += `${i + 1}. ${step.result.content.substring(0, 400)}...\n`;
      });
      prompt += `\n`;
    }

    prompt += `Provide your expert analysis and answer based on your knowledge and any previous context.`;
    return prompt;
  }

  buildSynthesisPrompt(userQuery, steps, chatHistory) {
    let prompt = `Synthesize a final comprehensive answer for: "${userQuery}"\n\n`;
    
    prompt += `Based on the following analysis steps:\n`;
    steps.forEach((step, i) => {
      prompt += `\nStep ${i + 1} (${step.action.type}):\n${step.result.content}\n`;
      if (step.result.sources) {
        prompt += `Sources: ${step.result.sources.map(s => s.url).join(', ')}\n`;
      }
    });

    prompt += `\nProvide a comprehensive final answer that synthesizes all the information above.`;
    return prompt;
  }

  parseActionPlan(content) {
    const lowerContent = content.toLowerCase();
    
    if (lowerContent.includes('search') && this.hasToolAvailable('search_web')) {
      return { type: 'search' };
    } else if (lowerContent.includes('code') && this.hasToolAvailable('run_code')) {
      return { type: 'code' };
    } else if (lowerContent.includes('synthesize')) {
      return { type: 'synthesize' };
    } else {
      return { type: 'reason' };
    }
  }

  async assessProgress(userQuery, steps) {
    if (steps.length === 0) return 0;
    
    // Simple confidence assessment based on step types and results
    const lastStep = steps[steps.length - 1];
    
    if (lastStep.action.type === 'synthesize') return 0.95;
    if (lastStep.action.type === 'code' && lastStep.result.executionResults) return 0.9;
    if (lastStep.action.type === 'search' && lastStep.result.sources) return 0.8;
    
    return Math.min(0.7 + (steps.length * 0.1), 0.9);
  }

  async callModel(messages, options = {}) {
    const response = await fetch(`${this.config.openrouter.baseUrl}/chat/completions`, {
      method: 'POST',
      headers: {
        'Authorization': `Bearer ${this.config.openrouter.apiKey}`,
        'Content-Type': 'application/json',
        'HTTP-Referer': this.config.openrouter.httpReferer,
        'X-Title': this.config.openrouter.xTitle
      },
      body: JSON.stringify({
        model: this.config.openrouter.model,
        messages,
        max_tokens: options.max_tokens || 1000,
        temperature: options.temperature || 0.7,
        tools: options.tools || undefined,
        tool_choice: options.tool_choice || undefined
      })
    });

    if (!response.ok) {
      throw new Error(`OpenRouter API error: ${response.status} ${response.statusText}`);
    }

    const data = await response.json();
    return data.choices[0].message;
  }

  formatSearchResults(searchResults) {
    if (!searchResults.results || searchResults.results.length === 0) {
      return 'No search results found.';
    }

    let formatted = 'Search Results:\n\n';
    searchResults.results.forEach((result, i) => {
      formatted += `${i + 1}. **${result.title}**\n`;
      formatted += `   ${result.snippet}\n`;
      formatted += `   Source: ${result.url}\n\n`;
    });

    return formatted;
  }

  extractSources(steps) {
    const sources = [];
    steps.forEach(step => {
      if (step.result.sources) {
        sources.push(...step.result.sources);
      }
    });
    return sources;
  }

  generateRequestId() {
    return Math.random().toString(36).substring(2, 15);
  }
}

module.exports = new ModularUnifiedAgent();