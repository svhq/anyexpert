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
    this.maxSteps = parseInt(process.env.MAX_PLANNING_STEPS) || 6;
    this.confidenceThreshold = parseFloat(process.env.CONFIDENCE_THRESHOLD) || 0.85;
    
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
        const toolMapping = {
          'search': 'search_web',
          'code': 'run_code',
          'scrape': 'scrape_web'
        };
        
        if (action.type !== 'reason' && action.type !== 'synthesize' && toolMapping[action.type]) {
          if (!this.hasToolAvailable(toolMapping[action.type])) {
            logger.info({
              requestId,
              stepNum,
              message: `Tool ${action.type} not available in API mode ${this.toolConfig.mode}, switching to reasoning`
            });
            action.type = 'reason';
          }
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
    // Build context from previous steps
    const context = steps.length > 0 ? 
      `\n\nPrevious steps:\n${steps.map(s => `- ${s.action.type}: ${s.result.content.substring(0, 200)}...`).join('\n')}` : '';
    
    // Dynamically determine available actions
    const availableActions = ['reason', 'synthesize'];
    if (this.hasToolAvailable('search_web')) availableActions.push('search');
    if (this.hasToolAvailable('run_code')) availableActions.push('code');
    
    // Use structured JSON planning
    const messages = [
      {
        role: 'system',
        content: this.buildPlanningSystemPrompt(availableActions)
      },
      {
        role: 'user',
        content: `Query: "${userQuery}"${context}\n\nWhat should be the next action?`
      }
    ];
    
    const response = await this.callModel(messages, {
      max_tokens: 200,
      temperature: 0.3
    });
    
    return this.parseActionResponse(response.content);
  }

  /**
   * Build planning system prompt with available actions
   */
  buildPlanningSystemPrompt(availableActions) {
    return `You are a planning assistant. Analyze the query and determine the next action(s).

Return ONLY valid JSON:
{
  "next_action": "${availableActions.join('|')}",
  "parallel_actions": ["action1", "action2"] or null,
  "rationale": "brief explanation"
}

Guidelines:
- search: Current information, facts, recent events
- code: Calculations, data analysis, algorithms
- reason: Analysis using existing knowledge
- synthesize: Final answer when sufficient information gathered
- scrape: Extract full content from specific URLs when needed

Parallel execution:
- Use when multiple independent tasks needed
- parallel_actions should include ALL actions to run in parallel
- Example: "Calculate X and find Y" → parallel_actions: ["code", "search"]
- If only one action needed → parallel_actions: null
- next_action should be the first action from parallel_actions when using parallel execution`;
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
        
      case 'scrape':
        if (!this.hasToolAvailable('scrape_web')) {
          logger.warn({ requestId, message: 'Scrape requested but not available, falling back to reasoning' });
          return this.executeReason(userQuery, steps, chatHistory);
        }
        return this.executeScrape(userQuery, steps, chatHistory, requestId);
        
      default:
        logger.warn({ requestId, message: `Unknown action type: ${action.type}, falling back to reasoning` });
        return this.executeReason(userQuery, steps, chatHistory);
    }
  }

  /**
   * Execute web search action
   */
  async executeSearch(userQuery, steps, requestId) {
    // Get previous search results if any
    const previousSearches = steps.filter(s => s.action?.type === 'search');
    const allPreviousResults = previousSearches.flatMap(s => s.result?.sources || []);

    // Generate search queries
    const round = previousSearches.length;
    const queries = await searchPlanner.generate(
      userQuery,
      'Need current information',
      round,
      allPreviousResults
    );
    
    logger.info({ 
      requestId, 
      type: 'web_search', 
      queries: queries 
    });

    // Execute searches using runBatch like original agent
    const searchResults = await webSearch.runBatch(queries);
    
    // Format results - matching original agent structure
    const sources = searchResults.map((result, idx) => ({
      number: idx + 1,
      title: result.title,
      url: result.url,
      snippet: result.snippet
    }));

    return {
      content: this.formatSearchResults(sources),
      sources,
      queries,
      confidence: Math.min(sources.length / 10, 0.7)
    };
  }

  /**
   * Execute code execution action
   */
  async executeCode(userQuery, steps, chatHistory, requestId) {
    // Support both code extraction and tool calling
    const useToolCalling = this.config.preferToolCalling ?? true;
    
    if (useToolCalling && this.hasToolAvailable('run_code')) {
      // Modern tool calling approach
      return await this.executeCodeViaTools(userQuery, steps, chatHistory, requestId);
    } else {
      // Classic code extraction approach
      return await this.executeCodeViaExtraction(userQuery, steps, chatHistory, requestId);
    }
  }

  async executeCodeViaTools(userQuery, steps, chatHistory, requestId) {
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
      max_tokens: 50000
    });
    
    if (response.tool_calls && response.tool_calls.length > 0) {
      const toolCall = response.tool_calls[0];
      const { code, timeout } = JSON.parse(toolCall.function.arguments);
      
      logger.info({ 
        requestId, 
        type: 'code_execution', 
        method: 'tool_calling',
        language: 'python',
        codeLength: code.length,
        success: true
      });
      
      const executionResult = await e2bManager.executeCode(code, { 
        timeout: timeout || 30000,
        userId: 'modular-agent'
      });
      
      // If E2B was skipped, adjust confidence
      const confidence = executionResult.skipped ? 0.6 : 0.9;
      
      return {
        content: response.content,
        code: code,
        executionResults: executionResult,
        tokensUsed: response.usage?.total_tokens || 0,
        confidence: confidence
      };
    }
    
    return {
      content: response.content,
      tokensUsed: response.usage?.total_tokens || 0,
      confidence: 0.6
    };
  }

  async executeCodeViaExtraction(userQuery, steps, chatHistory, requestId) {
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
      max_tokens: 50000,
      temperature: 0.1
    });
    
    // Extract Python code from response
    const codeMatch = response.content.match(/```python\n([\s\S]*?)\n```/);
    
    if (codeMatch && codeMatch[1]) {
      const code = codeMatch[1];
      
      logger.info({ 
        requestId, 
        type: 'code_execution', 
        method: 'code_extraction',
        language: 'python',
        codeLength: code.length,
        success: true
      });
      
      const executionResult = await e2bManager.executeCode(code, { 
        timeout: 30000,
        userId: 'modular-agent'
      });
      
      // If E2B was skipped, adjust confidence
      const confidence = executionResult.skipped ? 0.6 : 0.9;
      
      return {
        content: response.content,
        code: code,
        executionResults: executionResult,
        tokensUsed: response.usage?.total_tokens || 0,
        confidence: confidence
      };
    }
    
    // No code found
    return {
      content: response.content,
      tokensUsed: response.usage?.total_tokens || 0,
      confidence: 0.5
    };
  }

  /**
   * Execute web scraping action
   */
  async executeScrape(userQuery, steps, chatHistory, requestId) {
    const messages = [
      {
        role: 'system',
        content: this.generateSystemMessage()
      },
      {
        role: 'user',
        content: this.buildScrapePrompt(userQuery, steps, chatHistory)
      }
    ];
    
    const response = await this.callModel(messages, {
      tools: [this.toolDefinitions.find(t => t.function.name === 'scrape_web')],
      tool_choice: 'auto',
      max_tokens: 8000
    });
    
    if (response.tool_calls && response.tool_calls.length > 0) {
      const toolCall = response.tool_calls[0];
      const { url, reason } = JSON.parse(toolCall.function.arguments);
      
      logger.info({ 
        requestId, 
        type: 'web_scrape', 
        url,
        reason
      });
      
      // Execute scraping
      const scrapedContent = await webSearch.scrape(url);
      
      if (!scrapedContent) {
        return {
          content: `Unable to scrape content from ${url}. The page may be inaccessible or blocked.`,
          confidence: 0.2
        };
      }
      
      // Truncate content if too long (to avoid context overflow)
      const maxContentLength = 10000;
      let content = scrapedContent.content;
      if (content.length > maxContentLength) {
        content = content.substring(0, maxContentLength) + '\n\n[Content truncated...]';
      }
      
      return {
        content: `${response.content}\n\nScraped content from ${url}:\n\nTitle: ${scrapedContent.title}\n\n${content}`,
        source: {
          url,
          title: scrapedContent.title,
          metadata: scrapedContent.metadata
        },
        confidence: 0.95
      };
    }
    
    // No tool call made
    return {
      content: response.content,
      tokensUsed: response.usage?.total_tokens || 0,
      confidence: 0.5
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
      max_tokens: 8000,
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
      max_tokens: 8000,
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

  buildScrapePrompt(userQuery, steps, chatHistory) {
    let prompt = `As an expert, identify which URL needs to be scraped for: "${userQuery}"\n\n`;
    
    if (steps.length > 0) {
      prompt += `Context from previous steps:\n`;
      steps.forEach((step, i) => {
        if (step.result.sources) {
          prompt += `\nStep ${i + 1} found these sources:\n`;
          step.result.sources.forEach(source => {
            prompt += `- ${source.title} (${source.url})\n`;
          });
        }
      });
      prompt += `\n`;
    }

    prompt += `Use the scrape_web tool to extract full content from the most relevant URL. Choose wisely - scraping is more expensive than search.`;
    return prompt;
  }

  parseActionResponse(content) {
    try {
      // Extract JSON from markdown code blocks if present
      let jsonContent = content;
      const jsonMatch = content.match(/```json\n([\s\S]*?)\n```/);
      if (jsonMatch && jsonMatch[1]) {
        jsonContent = jsonMatch[1];
      }
      
      // Try JSON parsing
      const parsed = JSON.parse(jsonContent);
      
      // Validate and clean the response
      const action = {
        type: parsed.next_action || 'reason',
        parallelActions: Array.isArray(parsed.parallel_actions) ? parsed.parallel_actions : null,
        rationale: parsed.rationale || ''
      };
      
      // Ensure actions are available
      if (action.parallelActions) {
        action.parallelActions = action.parallelActions.filter(a => 
          this.isActionAvailable(a)
        );
        if (action.parallelActions.length === 0) {
          action.parallelActions = null;
        }
      }
      
      return action;
    } catch (e) {
      // Elegant fallback for non-JSON responses
      return this.parseActionFallback(content);
    }
  }

  parseActionFallback(content) {
    const lower = content.toLowerCase();
    
    // Check for parallel execution hints
    const hasMultipleTasks = 
      (lower.includes('both') || lower.includes('and also') || lower.includes('as well as')) &&
      (lower.includes('search') || lower.includes('find')) && 
      (lower.includes('calculate') || lower.includes('code'));
    
    if (hasMultipleTasks) {
      const actions = [];
      if (lower.includes('search') && this.hasToolAvailable('search_web')) actions.push('search');
      if (lower.includes('code') && this.hasToolAvailable('run_code')) actions.push('code');
      
      if (actions.length > 1) {
        return {
          type: actions[0],
          parallelActions: actions,
          rationale: 'Multiple independent tasks detected'
        };
      }
    }
    
    // Single action detection
    if (lower.includes('search') && this.hasToolAvailable('search_web')) {
      return { type: 'search', parallelActions: null };
    }
    if (lower.includes('code') && this.hasToolAvailable('run_code')) {
      return { type: 'code', parallelActions: null };
    }
    if (lower.includes('synthesize')) {
      return { type: 'synthesize', parallelActions: null };
    }
    
    return { type: 'reason', parallelActions: null };
  }

  isActionAvailable(actionType) {
    switch(actionType) {
      case 'search':
        return this.hasToolAvailable('search_web');
      case 'code':
        return this.hasToolAvailable('run_code');
      case 'reason':
      case 'synthesize':
        return true;
      default:
        return false;
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
        max_tokens: options.max_tokens || 16000,
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

  formatSearchResults(sources) {
    if (!sources || sources.length === 0) {
      return 'No search results found.';
    }

    let formatted = 'Search Results:\n\n';
    sources.forEach((source) => {
      formatted += `${source.number}. **${source.title}**\n`;
      formatted += `   ${source.snippet}\n`;
      formatted += `   Source: ${source.url}\n\n`;
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