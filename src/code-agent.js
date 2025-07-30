const config = require('../config');
const logger = require('./utils/logger');
const { SYSTEM_PROMPT } = require('./system-prompt');

/**
 * Code Agent - Handles code execution and explanation queries
 * Uses structured [WRITE]/[RUN] pattern similar to existing SEARCH workflow
 */
class CodeAgent {
  constructor() {
    this.config = config;
    this.codeExecutionUrl = config.codeExecution.serviceUrl;
    this.maxIterations = config.codeExecution.maxIterations;
    this.timeout = config.codeExecution.timeout;
  }

  /**
   * Handle code-related queries with optional execution
   * @param {string} userQuery - The user's question
   * @param {Object} codingIntent - Detected coding intent from query analyzer
   * @param {Object} chatHistory - Previous conversation context
   * @param {Object} metadata - Request metadata
   * @returns {Promise<Object>} - Response with code and/or execution results
   */
  async handle(userQuery, codingIntent, chatHistory = {}, metadata = {}) {
    try {
      if (codingIntent.type === 'code-explain') {
        // For explanation-only queries, use direct reasoning
        return await this.explainCode(userQuery, chatHistory, metadata);
      } else {
        // For execution queries, use the structured [WRITE]/[RUN] approach
        return await this.executeCode(userQuery, codingIntent, chatHistory, metadata);
      }
    } catch (error) {
      console.error('Code agent error:', error);
      throw error;
    }
  }

  /**
   * Handle code explanation without execution
   */
  async explainCode(userQuery, chatHistory, metadata) {
    const messages = [
      {
        role: 'system',
        content: `${SYSTEM_PROMPT}\n\nYou are now acting as a programming expert. Provide clear explanations of code concepts, syntax, and best practices. Focus on educational value and clarity.`
      }
    ];

    // Add chat history if available
    if (chatHistory.messages && chatHistory.messages.length > 0) {
      messages.push(...chatHistory.messages.slice(-4)); // Last 2 exchanges
    }

    messages.push({
      role: 'user',
      content: userQuery
    });

    const response = await this.callLLM(messages, metadata);
    
    return {
      content: response.content,
      codeExecuted: false,
      metadata: {
        ...metadata,
        agent: 'code-explain',
        tokensUsed: response.tokensUsed || 0
      }
    };
  }

  /**
   * Handle code execution using structured [WRITE]/[RUN] pattern
   */
  async executeCode(userQuery, codingIntent, chatHistory, metadata) {
    let iteration = 0;
    let executionResults = [];
    let finalCode = '';
    
    const messages = [
      {
        role: 'system',
        content: `${SYSTEM_PROMPT}

IMPORTANT: For code execution queries, you MUST follow this EXACT format:

Step 1: Write your code between [WRITE] tags:
[WRITE ${codingIntent.language || 'python'}]
\`\`\`${codingIntent.language || 'python'}
# Your code here
\`\`\`
[/WRITE]

Step 2: Execute it with:
[RUN]
[/RUN]

You'll receive output after [RUN]. If there are errors, fix and retry (max ${this.maxIterations} times).

Example for "What does [x**2 for x in range(3)] return?":

[WRITE python]
\`\`\`python
result = [x**2 for x in range(3)]
print(result)
\`\`\`
[/WRITE]

[RUN]
[/RUN]`
      }
    ];

    // Add chat history if available
    if (chatHistory.messages && chatHistory.messages.length > 0) {
      messages.push(...chatHistory.messages.slice(-2)); // Last exchange
    }

    messages.push({
      role: 'user',
      content: `${userQuery}

Remember: You MUST use the [WRITE]...[/WRITE] and [RUN] tags to execute code. Do not just explain the answer - actually run the code to show the output.`
    });

    while (iteration < this.maxIterations) {
      iteration++;
      
      // Get LLM response with code
      const llmResponse = await this.callLLM(messages, metadata);
      const content = llmResponse.content;
      
      // Debug: log what LLM returned
      console.log('LLM Response:', content.substring(0, 500));
      
      // Extract code from [WRITE] blocks
      const codeMatch = content.match(/\[WRITE[^\]]*\]\s*```[^`]*\n([\s\S]*?)```\s*\[\/WRITE\]/);
      if (!codeMatch) {
        // No code found, return explanation
        return {
          content,
          codeExecuted: false,
          metadata: {
            ...metadata,
            agent: 'code-agent',
            iterations: iteration,
            tokensUsed: llmResponse.tokensUsed || 0
          }
        };
      }
      
      const code = codeMatch[1].trim();
      finalCode = code;
      
      // Check if execution is requested
      const hasRunTag = content.includes('[RUN]');
      if (!hasRunTag) {
        // No execution requested, return with code
        return {
          content,
          code: finalCode,
          codeExecuted: false,
          metadata: {
            ...metadata,
            agent: 'code-agent',
            iterations: iteration,
            tokensUsed: llmResponse.tokensUsed || 0
          }
        };
      }
      
      // Execute the code
      const executionResult = await this.runCode(
        codingIntent.language || 'python',
        code,
        metadata.requestId
      );
      
      executionResults.push({
        iteration,
        code,
        result: executionResult
      });
      
      // If execution successful and no errors, we're done
      if (executionResult.exitCode === 0 && !executionResult.stderr) {
        const successContent = this.formatSuccessResponse(content, code, executionResult);
        return {
          content: successContent,
          code: finalCode,
          executionResults,
          codeExecuted: true,
          metadata: {
            ...metadata,
            agent: 'code-agent',
            iterations: iteration,
            tokensUsed: llmResponse.tokensUsed || 0
          }
        };
      }
      
      // If there were errors and we have iterations left, continue
      if (iteration < this.maxIterations) {
        messages.push(
          { role: 'assistant', content },
          { 
            role: 'user', 
            content: `Execution result:\nSTDOUT: ${executionResult.stdout}\nSTDERR: ${executionResult.stderr}\nExit Code: ${executionResult.exitCode}\n\nPlease fix the issue and try again.`
          }
        );
      }
    }
    
    // Max iterations reached, return final attempt
    const finalContent = this.formatFinalResponse(messages[messages.length - 2].content, finalCode, executionResults);
    return {
      content: finalContent,
      code: finalCode,
      executionResults,
      codeExecuted: true,
      metadata: {
        ...metadata,
        agent: 'code-agent',
        iterations: this.maxIterations,
        tokensUsed: 0,
        maxIterationsReached: true
      }
    };
  }

  /**
   * Execute code via the run_code microservice
   */
  async runCode(language, source, requestId) {
    try {
      const response = await fetch(`${this.codeExecutionUrl}/run_code`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          language,
          source,
          timeout: this.timeout
        })
      });

      if (!response.ok) {
        throw new Error(`Code execution service error: ${response.status}`);
      }

      const result = await response.json();
      
      // Log execution result
      logger.logCodeExecution(requestId, language, source.length, result.exitCode === 0);
      
      return result;
    } catch (error) {
      console.error('Code execution error:', error);
      return {
        stdout: '',
        stderr: `Execution service error: ${error.message}`,
        exitCode: 1
      };
    }
  }

  /**
   * Call LLM for code generation/explanation
   */
  async callLLM(messages, metadata) {
    const response = await fetch(`${this.config.openrouter.baseUrl}/chat/completions`, {
      method: 'POST',
      headers: this.config.openrouter.headers,
      body: JSON.stringify({
        model: this.config.openrouter.model,
        messages,
        temperature: 0.1 // Low temperature for consistent code generation
      })
    });

    if (!response.ok) {
      throw new Error(`OpenRouter API error: ${response.status}`);
    }

    const data = await response.json();
    return {
      content: data.choices[0].message.content,
      tokensUsed: data.usage?.total_tokens || 0
    };
  }

  /**
   * Format successful execution response
   */
  formatSuccessResponse(originalContent, code, executionResult) {
    const codeBlock = `\`\`\`\n${code}\n\`\`\``;
    const output = executionResult.stdout ? `\n\n**Output:**\n\`\`\`\n${executionResult.stdout}\n\`\`\`` : '';
    
    // Remove the [WRITE]/[RUN] tags from display and add formatted result
    let cleanContent = originalContent
      .replace(/\[WRITE[^\]]*\][\s\S]*?\[\/WRITE\]/g, '')
      .replace(/\[RUN\][\s\S]*?\[\/RUN\]/g, '')
      .trim();
    
    return `${cleanContent}\n\n**Code:**\n${codeBlock}${output}`;
  }

  /**
   * Format final response when max iterations reached
   */
  formatFinalResponse(originalContent, code, executionResults) {
    const lastResult = executionResults[executionResults.length - 1];
    const codeBlock = `\`\`\`\n${code}\n\`\`\``;
    
    let output = '';
    if (lastResult) {
      if (lastResult.result.stdout) {
        output += `\n\n**Output:**\n\`\`\`\n${lastResult.result.stdout}\n\`\`\``;
      }
      if (lastResult.result.stderr) {
        output += `\n\n**Errors:**\n\`\`\`\n${lastResult.result.stderr}\n\`\`\``;
      }
    }
    
    let cleanContent = originalContent
      .replace(/\[WRITE[^\]]*\][\s\S]*?\[\/WRITE\]/g, '')
      .replace(/\[RUN\][\s\S]*?\[\/RUN\]/g, '')
      .trim();
    
    return `${cleanContent}\n\n**Final Code:**\n${codeBlock}${output}`;
  }
}

// Export singleton instance
module.exports = new CodeAgent();