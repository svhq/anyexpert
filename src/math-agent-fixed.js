const config = require('../config');
const logger = require('./utils/logger');
const { SYSTEM_PROMPT } = require('./system-prompt');

/**
 * Math Agent - Handles mathematical computations and explanations
 * Fixed version with better error handling for JSON parsing
 */
class MathAgent {
  constructor() {
    this.config = config;
    this.codeExecutionUrl = config.codeExecution.serviceUrl;
  }

  /**
   * Handle math-related queries
   */
  async handle(userQuery, mathIntent, chatHistory = {}, metadata = {}) {
    try {
      if (!mathIntent.needsExecution) {
        // Simple math - let the LLM handle it directly
        return await this.handleSimpleMath(userQuery, chatHistory, metadata);
      } else {
        // Complex math - use Python execution
        return await this.handleComplexMath(userQuery, mathIntent, chatHistory, metadata);
      }
    } catch (error) {
      console.error('Math agent error:', error);
      // Fallback to direct answer without code execution
      return await this.handleSimpleMath(userQuery, chatHistory, metadata);
    }
  }

  /**
   * Handle simple math without execution
   */
  async handleSimpleMath(userQuery, chatHistory, metadata) {
    const messages = [
      {
        role: 'system',
        content: `${SYSTEM_PROMPT}\n\nSolve mathematical problems step-by-step, showing your work clearly. For simple arithmetic, compute the answer directly.`
      }
    ];

    if (chatHistory.messages && chatHistory.messages.length > 0) {
      messages.push(...chatHistory.messages.slice(-4));
    }

    messages.push({ role: 'user', content: userQuery });

    const response = await this.callLLMSafe(messages, metadata);
    return response;
  }

  /**
   * Handle complex math with Python execution
   */
  async handleComplexMath(userQuery, mathIntent, chatHistory, metadata) {
    const systemPrompt = `${SYSTEM_PROMPT}

You are solving a mathematical problem that requires computation. 
First, explain your approach clearly.
Then, write Python code to solve the problem precisely.
Use numpy, scipy, or sympy as needed for mathematical operations.

Format your response as:
1. Clear explanation of the problem and approach
2. Python code block with the solution
3. Interpretation of the results`;

    const messages = [
      { role: 'system', content: systemPrompt }
    ];

    if (chatHistory.messages && chatHistory.messages.length > 0) {
      messages.push(...chatHistory.messages.slice(-2));
    }

    messages.push({ role: 'user', content: userQuery });

    try {
      // Get LLM response with mathematical solution
      const llmResponse = await this.callLLMSafe(messages, metadata);
      
      console.log('Math agent - LLM response received, length:', llmResponse.content.length);

      // Extract Python code from response
      const codeMatch = llmResponse.content.match(/```python\n([\s\S]*?)```/);
      if (!codeMatch) {
        // No code found, return as-is
        return llmResponse;
      }

      const code = codeMatch[1].trim();
      console.log('Math agent - executing python code (' + code.length + ' chars)');
      console.log('First 200 chars:', code.substring(0, 200));

      // Execute the Python code
      const executionResult = await this.executePython(code, metadata);
      
      console.log('Math agent - execution result:', { 
        exitCode: executionResult.exitCode, 
        stdoutLength: executionResult.stdout?.length || 0,
        stderrLength: executionResult.stderr?.length || 0
      });

      // Format the response with execution results
      const formattedContent = this.formatMathResponse(
        llmResponse.content,
        code,
        executionResult
      );

      return {
        content: formattedContent,
        tokensUsed: llmResponse.tokensUsed
      };

    } catch (error) {
      console.error('Complex math execution error:', error);
      // Return the LLM response without execution
      return await this.callLLMSafe(messages, metadata);
    }
  }

  /**
   * Execute Python code using the code execution service
   */
  async executePython(code, metadata) {
    try {
      const response = await fetch(this.codeExecutionUrl + '/execute', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          code,
          language: 'python',
          timeout: 30000
        })
      });

      if (!response.ok) {
        throw new Error(`Code execution service error: ${response.status}`);
      }

      return await response.json();
    } catch (error) {
      logger.error({ 
        error: error.message, 
        metadata,
        context: 'python-execution' 
      });
      
      return {
        exitCode: 1,
        stdout: '',
        stderr: `Execution error: ${error.message}`
      };
    }
  }

  /**
   * Safe LLM call with timeout and error handling
   */
  async callLLMSafe(messages, metadata) {
    const timeout = 60000; // 60 seconds timeout
    const controller = new AbortController();
    const timeoutId = setTimeout(() => controller.abort(), timeout);

    try {
      const response = await fetch(`${this.config.openrouter.baseUrl}/chat/completions`, {
        method: 'POST',
        headers: this.config.openrouter.headers,
        body: JSON.stringify({
          model: this.config.openrouter.model,
          messages,
          temperature: 0.1,
          max_tokens: 4000 // Limit response size
        }),
        signal: controller.signal
      });

      clearTimeout(timeoutId);

      if (!response.ok) {
        const errorText = await response.text().catch(() => 'Unknown error');
        throw new Error(`OpenRouter API error ${response.status}: ${errorText}`);
      }

      // Read response as text first to handle potential JSON errors
      const responseText = await response.text();
      
      let data;
      try {
        data = JSON.parse(responseText);
      } catch (parseError) {
        console.error('JSON parse error:', parseError);
        console.error('Response text (first 500 chars):', responseText.substring(0, 500));
        throw new Error(`Invalid JSON response: ${parseError.message}`);
      }

      if (!data.choices || !data.choices[0] || !data.choices[0].message) {
        throw new Error('Invalid response structure from OpenRouter');
      }

      return {
        content: data.choices[0].message.content,
        tokensUsed: data.usage?.total_tokens || 0
      };

    } catch (error) {
      clearTimeout(timeoutId);
      
      if (error.name === 'AbortError') {
        throw new Error('Request timeout after 60 seconds');
      }
      
      throw error;
    }
  }

  /**
   * Format response with mathematical output
   */
  formatMathResponse(originalContent, code, executionResult) {
    // Extract the explanation part (before code)
    const parts = originalContent.split('```python');
    let explanation = parts[0].trim();
    
    // Build formatted response
    let formatted = explanation + '\n\n';
    formatted += '**Calculation:**\n```python\n' + code + '\n```\n\n';
    
    if (executionResult.stdout) {
      formatted += '**Result:**\n```\n' + executionResult.stdout.trim() + '\n```';
      
      // Try to extract and highlight the final answer
      const output = executionResult.stdout.trim();
      const lastLine = output.split('\n').pop();
      
      // Check if the last line looks like a final answer
      if (lastLine && (lastLine.includes('=') || lastLine.match(/^[\d\.\-\+e]+$/))) {
        formatted += '\n\n**Answer: ' + lastLine.replace(/.*=\s*/, '') + '**';
      }
    }
    
    if (executionResult.stderr) {
      formatted += '\n\n**Error:**\n```\n' + executionResult.stderr + '\n```';
    }
    
    return formatted;
  }
}

module.exports = new MathAgent();