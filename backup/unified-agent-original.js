const config = require('../config');
const logger = require('./utils/logger');
const { SYSTEM_PROMPT } = require('./system-prompt');
const webSearch = require('./web-search');
const searchPlanner = require('./search-planner');
const e2bManager = require('./e2b-manager-v3');

/**
 * Unified Agent - Combines all capabilities with multi-step processing
 * Supports web search, code execution, and direct reasoning in a loop
 */
class UnifiedAgent {
  constructor() {
    this.config = config;
    this.maxSteps = 6;
    this.confidenceThreshold = 0.9;
    
    // Define available tools for OpenRouter API
    this.toolDefinitions = [
      {
        type: 'function',
        function: {
          name: 'search_web',
          description: 'Search the web for current information, facts, and recent events. Returns relevant snippets and sources.',
          parameters: {
            type: 'object',
            properties: {
              query: { 
                type: 'string',
                description: 'The search query'
              },
              count: { 
                type: 'integer',
                description: 'Number of results to return',
                default: 5
              }
            },
            required: ['query']
          }
        }
      },
      {
        type: 'function',
        function: {
          name: 'run_code',
          description: 'Execute Python code in an isolated sandbox. Available libraries: numpy, scipy, sympy, mpmath, pandas, matplotlib, seaborn, scikit-learn, biopython, splicekit, deeptools, beautifulsoup4, requests.',
          parameters: {
            type: 'object',
            properties: {
              code: { 
                type: 'string',
                description: 'Python code to execute'
              },
              timeout: { 
                type: 'integer',
                description: 'Execution timeout in milliseconds',
                default: 30000
              }
            },
            required: ['code']
          }
        }
      }
    ];
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
      step: 'unified-agent-start', 
      query: userQuery 
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

        // Plan next action based on current state
        const action = await this.planNextAction(userQuery, steps, chatHistory);
        
        // Execute the action(s)
        let result;
        
        // Check for parallel execution based on model's suggestion
        if (action.parallelActions && action.parallelActions.length > 1) {
          logger.info({ 
            requestId, 
            stepNum,
            status: 'parallel-execution',
            actions: action.parallelActions
          });
          
          // Execute actions in parallel
          const parallelPromises = action.parallelActions.map(actionType => {
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
              `${action.parallelActions[i].toUpperCase()} Results:\n${r.content}`
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
          action.executedActions = action.parallelActions;
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
        tokensUsed: steps.reduce((sum, s) => sum + (s.result.tokensUsed || 0), 0)
      };
      
      logger.info({ type: 'unified-agent-complete', ...finalMetrics });

      return {
        content: finalAnswer.content,
        metadata: {
          ...finalMetrics,
          steps: steps.map(s => ({
            action: s.action,  // Return full action object, not just type
            duration: s.timestamp - startTime
          }))
        }
      };

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
   * Plan the next action based on current state
   */
  async planNextAction(userQuery, steps, chatHistory) {
    // Build context from previous steps
    let context = '';
    if (steps.length > 0) {
      context = '\n\nPrevious steps taken:\n';
      steps.forEach((step, idx) => {
        context += `${idx + 1}. ${step.action.type}: ${step.action.rationale}\n`;
        if (step.result.confidence) {
          context += `   Confidence: ${step.result.confidence}\n`;
        }
      });
    }

    const messages = [
      {
        role: 'system',
        content: `You are a planning assistant. Analyze the query and determine the next action(s).

Return ONLY valid JSON with these keys:
{
  "next_action": "search|code|reason|synthesize",
  "parallel_actions": ["search", "code"] or null,
  "rationale": "<1-2 lines explaining why>"
}

Guidelines:
- "search": Use for current info, facts, recent events, external validation
- "code": Use for calculations, verification, data analysis, algorithms
- "reason": Use for analysis without external tools
- "synthesize": Use when you have enough info to provide final answer

For parallel execution:
- If query needs BOTH search AND code that can work independently, set: parallel_actions: ["search", "code"]
- If query needs multiple searches, set: parallel_actions: ["search", "search"]
- If only one action needed, set: parallel_actions: null
- Example: "Calculate X and find Y" → parallel_actions: ["code", "search"]`
      },
      {
        role: 'user',
        content: `Query: "${userQuery}"${context}\n\nWhat should be the next action?`
      }
    ];

    try {
      const response = await this.callLLM(messages, { temperature: 0.1 });
      const planText = response.content;
      
      // Extract JSON from response
      let plan;
      const jsonMatch = planText.match(/\{[\s\S]*\}/);
      if (jsonMatch) {
        try {
          plan = JSON.parse(jsonMatch[0]);
        } catch (e) {
          logger.warn('Failed to parse planning JSON, using fallback');
          plan = this.fallbackPlan(userQuery, steps);
        }
      } else {
        plan = this.fallbackPlan(userQuery, steps);
      }

      // Validate and return plan
      const validActions = ['search', 'code', 'reason', 'synthesize'];
      if (!validActions.includes(plan.next_action)) {
        plan.next_action = 'reason';
      }

      return {
        type: plan.next_action,
        rationale: plan.rationale || 'Strategic action selection',
        parallelActions: plan.parallel_actions || null
      };

    } catch (error) {
      logger.error('Planning failed, using fallback', error);
      return this.fallbackPlan(userQuery, steps);
    }
  }

  /**
   * Fallback planning logic
   */
  fallbackPlan(userQuery, steps) {
    if (steps.length === 0) {
      // First step heuristics
      const queryLower = userQuery.toLowerCase();
      if (queryLower.includes('calculat') || queryLower.includes('comput') || 
          queryLower.includes('solve') || queryLower.includes('derive')) {
        return { type: 'code', rationale: 'Mathematical computation needed' };
      } else if (queryLower.includes('current') || queryLower.includes('latest') ||
                 queryLower.includes('recent') || queryLower.includes('2024') || 
                 queryLower.includes('2025')) {
        return { type: 'search', rationale: 'Current information needed' };
      }
      return { type: 'reason', rationale: 'Direct analysis' };
    }

    // Subsequent steps
    if (steps.length >= 2 || (steps.length >= 1 && steps[0].result.confidence >= 0.8)) {
      return { type: 'synthesize', rationale: 'Sufficient information gathered' };
    }
    
    return { type: 'reason', rationale: 'Continue analysis' };
  }

  /**
   * Execute the planned action
   */
  async executeAction(action, userQuery, steps, chatHistory, requestId) {
    switch (action.type) {
      case 'search':
        return await this.executeSearch(userQuery, steps, requestId);
      
      case 'code':
        return await this.executeCode(userQuery, steps, chatHistory, requestId);
      
      case 'reason':
        return await this.executeReasoning(userQuery, steps, chatHistory, requestId);
      
      case 'synthesize':
        return await this.synthesizeFinalAnswer(userQuery, steps, chatHistory, requestId);
      
      default:
        throw new Error(`Unknown action type: ${action.type}`);
    }
  }

  /**
   * Execute web search
   */
  async executeSearch(userQuery, steps, requestId) {
    // Get previous search results if any
    const previousSearches = steps.filter(s => s.action.type === 'search');
    const allPreviousResults = previousSearches.flatMap(s => s.result.sources || []);

    // Generate search queries
    const round = previousSearches.length;
    const queries = await searchPlanner.generate(
      userQuery,
      'Need current information',
      round,
      allPreviousResults
    );

    // Execute searches
    const searchResults = await webSearch.runBatch(queries);
    
    // Format results
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
   * Execute code
   */
  async executeCode(userQuery, steps, chatHistory, requestId) {
    // Build context from previous steps
    let context = '';
    for (const step of steps) {
      if (step.action.type === 'search' && step.result.sources?.length > 0) {
        context += '\nSearch findings:\n' + step.result.content + '\n';
      } else if (step.action.type === 'reason') {
        context += '\nAnalysis:\n' + step.result.content + '\n';
      }
    }

    const messages = [
      {
        role: 'system',
        content: `${SYSTEM_PROMPT}

Available Python libraries:
- math, decimal, fractions, mpmath (mathematical operations)
- numpy, scipy, numpy_financial (numerical computing)
- sympy (symbolic math)
- pandas, statistics, statsmodels (data analysis)
- matplotlib, seaborn, plotly (visualization)
- scikit-learn (machine learning)
- biopython (bioinformatics, sequence analysis)
- splicekit (RNA splicing analysis)
- deeptools (genomics data visualization)
- beautifulsoup4, requests (web scraping if needed)

Write Python code to solve or verify the problem. Include clear comments.`
      }
    ];

    if (context) {
      messages.push({
        role: 'user',
        content: `Context from previous analysis:\n${context}\n\nQuestion: ${userQuery}`
      });
    } else {
      messages.push({
        role: 'user',
        content: userQuery
      });
    }

    const llmResponse = await this.callLLM(messages, { temperature: 0.1 });
    const content = llmResponse.content;
    
    // Extract Python code
    let code = null;
    const codeBlockMatch = content.match(/```python\n([\s\S]*?)```/);
    if (codeBlockMatch) {
      code = codeBlockMatch[1].trim();
    }
    
    // Execute if code found
    if (code) {
      try {
        const executionResult = await this.runCode('python', code, requestId);
        
        if (executionResult.exitCode === 0) {
          return {
            content: this.formatCodeResponse(content, code, executionResult),
            code,
            executionResults: [{
              iteration: 1,
              code,
              result: executionResult
            }],
            confidence: 0.95, // Boost confidence for successful code execution
            tokensUsed: llmResponse.tokensUsed || 0
          };
        }
      } catch (error) {
        console.error('Code execution failed:', error);
      }
    }
    
    // Fallback: return explanation without execution
    return {
      content,
      code: code || 'No code extracted',
      confidence: 0.5,
      tokensUsed: llmResponse.tokensUsed || 0
    };
  }

  /**
   * Execute reasoning step
   */
  async executeReasoning(userQuery, steps, chatHistory, requestId) {
    // Build context from previous steps
    let context = '';
    for (const step of steps) {
      if (step.result.content) {
        context += `\n${step.action.type} results:\n${step.result.content}\n`;
      }
    }

    const messages = [
      {
        role: 'system',
        content: SYSTEM_PROMPT
      },
      {
        role: 'user',
        content: context ? 
          `Based on the following findings, analyze: ${userQuery}\n\nFindings:${context}` :
          userQuery
      }
    ];

    const response = await this.callLLM(messages, { temperature: 0.7 });
    
    return {
      content: response.content,
      confidence: 0.6,
      tokensUsed: response.tokensUsed || 0
    };
  }

  /**
   * Synthesize final answer from all steps
   */
  async synthesizeFinalAnswer(userQuery, steps, chatHistory, requestId) {
    // Collect all findings
    const findings = {
      searches: [],
      calculations: [],
      reasoning: []
    };

    for (const step of steps) {
      if (step.action.type === 'search' && step.result.sources) {
        findings.searches.push({
          sources: step.result.sources,
          content: step.result.content
        });
      } else if (step.action.type === 'code' && step.result.executionResults) {
        findings.calculations.push({
          code: step.result.code,
          output: step.result.executionResults[0]?.result?.stdout
        });
      } else if (step.action.type === 'reason') {
        findings.reasoning.push(step.result.content);
      }
    }

    // Build comprehensive context
    let context = '';
    
    if (findings.searches.length > 0) {
      context += '\nWeb Search Results:\n';
      findings.searches.forEach((search, idx) => {
        search.sources.forEach(source => {
          context += `[${source.number}] ${source.title}\n${source.snippet}\n\n`;
        });
      });
    }

    if (findings.calculations.length > 0) {
      context += '\nCalculations Performed:\n';
      findings.calculations.forEach(calc => {
        if (calc.output) {
          context += `Output: ${calc.output}\n\n`;
        }
      });
    }

    if (findings.reasoning.length > 0) {
      context += '\nAnalysis:\n';
      context += findings.reasoning.join('\n\n');
    }

    const messages = [
      {
        role: 'system',
        content: `${SYSTEM_PROMPT}

You have access to web search results and calculations. When using information from these sources, cite them using [1], [2], etc. format.`
      },
      {
        role: 'user',
        content: `Based on all the following research and analysis, please provide a comprehensive answer to: "${userQuery}"

${context}

Provide your final expert answer, citing sources where appropriate.`
      }
    ];

    const response = await this.callLLM(messages, { temperature: 0.7 });
    
    return {
      content: response.content,
      confidence: 0.9,
      tokensUsed: response.tokensUsed || 0,
      sources: findings.searches.flatMap(s => s.sources)
    };
  }

  /**
   * Assess progress and confidence
   */
  async assessProgress(userQuery, steps) {
    if (steps.length === 0) return 0;
    
    const lastStep = steps[steps.length - 1];
    
    // Simple heuristic based on action type and results
    let confidence = lastStep.result.confidence || 0.5;
    
    // Boost confidence if we have multiple types of verification
    const actionTypes = new Set(steps.map(s => s.action.type));
    if (actionTypes.size >= 2) {
      confidence = Math.min(confidence + 0.1, 1.0);
    }
    
    // Boost if we have both search and code
    if (actionTypes.has('search') && actionTypes.has('code')) {
      confidence = Math.min(confidence + 0.1, 1.0);
    }
    
    return confidence;
  }

  /**
   * Format search results for display
   */
  formatSearchResults(sources) {
    if (!sources || sources.length === 0) {
      return 'No search results found.';
    }

    return sources.map(source => 
      `[${source.number}] ${source.title}\n${source.snippet}`
    ).join('\n\n');
  }

  /**
   * Format code response with output
   */
  formatCodeResponse(explanation, code, executionResult) {
    // Extract the explanation part (before code)
    const parts = explanation.split('```python');
    let formatted = parts[0].trim() + '\n\n';
    
    formatted += '**Calculation:**\n```python\n' + code + '\n```\n\n';
    
    if (executionResult.stdout) {
      formatted += '**Result:**\n```\n' + executionResult.stdout.trim() + '\n```';
    }
    
    if (executionResult.stderr) {
      formatted += '\n\n**Error:**\n```\n' + executionResult.stderr + '\n```';
    }
    
    return formatted;
  }

  /**
   * Execute code via the new E2B manager
   */
  async runCode(language, source, requestId) {
    const result = await e2bManager.executeWithFallback(source, {
      language,
      timeoutMs: 120000,  // Increased to 2 minutes
      requestId,
      fallbackToLocal: true // Enable local Python fallback
    });
    
    logger.logCodeExecution(requestId, language, source.length, result.success);
    
    return {
      stdout: result.stdout || '',
      stderr: result.stderr || '',
      exitCode: result.exitCode
    };
  }

  /**
   * Call LLM with optional tool support
   */
  async callLLM(messages, options = {}) {
    const requestBody = {
      model: this.config.openrouter.model,
      messages,
      temperature: options.temperature || 0.7,
      max_tokens: options.max_tokens || 50000  // Reasonable limit
    };

    // Add tool support if requested
    if (options.tools) {
      requestBody.tools = options.tools;
      requestBody.tool_choice = options.tool_choice || 'auto';
      requestBody.parallel_tool_calls = options.parallel_tool_calls || false;
    }

    const response = await fetch(`${this.config.openrouter.baseUrl}/chat/completions`, {
      method: 'POST',
      headers: this.config.openrouter.headers,
      body: JSON.stringify(requestBody)
    });

    if (!response.ok) {
      throw new Error(`OpenRouter API error: ${response.status}`);
    }

    const data = await response.json();
    
    // Handle tool calls if present
    if (data.choices[0].message.tool_calls) {
      return {
        content: data.choices[0].message.content || '',
        tool_calls: data.choices[0].message.tool_calls,
        tokensUsed: data.usage?.total_tokens || 0
      };
    }

    return {
      content: data.choices[0].message.content,
      tokensUsed: data.usage?.total_tokens || 0
    };
  }

  /**
   * Generate unique request ID
   */
  generateRequestId() {
    return `req-${Date.now()}-${Math.random().toString(36).substr(2, 9)}`;
  }
}

// Export singleton instance
module.exports = new UnifiedAgent();