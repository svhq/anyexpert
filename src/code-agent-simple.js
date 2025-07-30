const config = require('../config');
const logger = require('./utils/logger');
const { SYSTEM_PROMPT } = require('./system-prompt');

/**
 * Simplified Code Agent - Works with GLM-4.5-air's limitations
 * Falls back to direct code extraction when structured format fails
 */
class SimpleCodeAgent {
  constructor() {
    this.config = config;
    this.codeExecutionUrl = config.codeExecution.serviceUrl;
  }

  /**
   * Handle code-related queries
   */
  async handle(userQuery, codingIntent, chatHistory = {}, metadata = {}) {
    try {
      if (codingIntent.type === 'code-explain') {
        return await this.explainCode(userQuery, chatHistory, metadata);
      } else {
        return await this.executeCode(userQuery, codingIntent, chatHistory, metadata);
      }
    } catch (error) {
      console.error('Code agent error:', error);
      // Fallback to explanation on error
      return await this.explainCode(userQuery, chatHistory, metadata);
    }
  }

  /**
   * Handle code explanation without execution
   */
  async explainCode(userQuery, chatHistory, metadata) {
    const messages = [
      {
        role: 'system',
        content: `${SYSTEM_PROMPT}\n\nProvide clear explanations of code concepts, syntax, and best practices.`
      }
    ];

    if (chatHistory.messages && chatHistory.messages.length > 0) {
      messages.push(...chatHistory.messages.slice(-2));
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
   * Handle code execution with fallback approach
   */
  async executeCode(userQuery, codingIntent, chatHistory, metadata) {
    const messages = [
      {
        role: 'system',
        content: `${SYSTEM_PROMPT}\n\nWhen asked about code output, write the code that would produce the answer.`
      }
    ];

    if (chatHistory.messages && chatHistory.messages.length > 0) {
      messages.push(...chatHistory.messages.slice(-2));
    }

    messages.push({
      role: 'user',
      content: `${userQuery}\n\nProvide the code that would answer this question. Include the code in a code block.`
    });

    const llmResponse = await this.callLLM(messages, metadata);
    const content = llmResponse.content;
    
    // Try to extract code from various formats
    let code = null;
    let language = codingIntent.language || 'python';
    
    // Try markdown code blocks
    const codeBlockMatch = content.match(/```(\w+)?\n([\s\S]*?)```/);
    if (codeBlockMatch) {
      if (codeBlockMatch[1]) {
        language = codeBlockMatch[1];
      }
      code = codeBlockMatch[2].trim();
    } else {
      // Try to extract code heuristically
      const lines = content.split('\n');
      const codeLines = [];
      let inCode = false;
      
      for (const line of lines) {
        // Common code patterns
        if (line.includes('print(') || line.includes('console.log(') || 
            line.includes('def ') || line.includes('function ') ||
            line.includes('=') || line.includes(';')) {
          inCode = true;
        }
        
        if (inCode && line.trim()) {
          codeLines.push(line);
        }
      }
      
      if (codeLines.length > 0) {
        code = codeLines.join('\n');
      }
    }
    
    // If we found code, try to execute it
    if (code) {
      try {
        const executionResult = await this.runCode(language, code, metadata.requestId);
        
        if (executionResult.exitCode === 0) {
          // Success - format response with output
          const formattedContent = this.formatResponse(content, code, executionResult);
          return {
            content: formattedContent,
            code,
            executionResults: [{
              iteration: 1,
              code,
              result: executionResult
            }],
            codeExecuted: true,
            metadata: {
              ...metadata,
              agent: 'code-agent-simple',
              tokensUsed: llmResponse.tokensUsed || 0
            }
          };
        }
      } catch (error) {
        console.error('Execution failed:', error);
      }
    }
    
    // Fallback: return explanation without execution
    return {
      content,
      code: code || 'No code extracted',
      codeExecuted: false,
      metadata: {
        ...metadata,
        agent: 'code-agent-simple',
        tokensUsed: llmResponse.tokensUsed || 0,
        fallback: true
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
          timeout: this.config.codeExecution.timeout
        })
      });

      const result = await response.json();
      
      if (!response.ok) {
        throw new Error(`Service error: ${result.error || response.status}`);
      }
      
      logger.logCodeExecution(requestId, language, source.length, result.exitCode === 0);
      return result;
    } catch (error) {
      console.error('Code execution error:', error);
      return {
        stdout: '',
        stderr: error.message,
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
        temperature: 0.1
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
   * Format response with execution output
   */
  formatResponse(originalContent, code, executionResult) {
    // Remove the code block from content if it's there
    let cleanContent = originalContent.replace(/```[\s\S]*?```/g, '').trim();
    
    // Build formatted response
    let formatted = cleanContent + '\n\n';
    formatted += '**Code:**\n```' + '\n' + code + '\n```\n\n';
    
    if (executionResult.stdout) {
      formatted += '**Output:**\n```\n' + executionResult.stdout + '```';
    }
    
    if (executionResult.stderr) {
      formatted += '\n\n**Error:**\n```\n' + executionResult.stderr + '```';
    }
    
    return formatted;
  }
}

// Export singleton instance
module.exports = new SimpleCodeAgent();