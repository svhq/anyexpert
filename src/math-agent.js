const config = require('../config');
const logger = require('./utils/logger');
const { SYSTEM_PROMPT } = require('./system-prompt');

/**
 * Math Agent - Handles mathematical computations and explanations
 * Uses E2B sandbox for high-precision calculations and symbolic math
 */
class MathAgent {
  constructor() {
    this.config = config;
    this.codeExecutionUrl = config.codeExecution.serviceUrl;
  }

  /**
   * Handle math-related queries
   * @param {string} userQuery - The user's question
   * @param {Object} mathIntent - Detected math intent from query analyzer
   * @param {Object} chatHistory - Previous conversation context
   * @param {Object} metadata - Request metadata
   * @returns {Promise<Object>} - Response with calculation and/or explanation
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
      // Fallback to explanation without execution
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
        content: SYSTEM_PROMPT
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
      mathExecuted: false,
      metadata: {
        ...metadata,
        agent: 'math-simple',
        tokensUsed: response.tokensUsed || 0
      }
    };
  }

  /**
   * Handle complex math with Python execution
   */
  async handleComplexMath(userQuery, mathIntent, chatHistory, metadata) {
    const messages = [
      {
        role: 'system',
        content: `${SYSTEM_PROMPT}

Available Python libraries for calculations:
- math, decimal, fractions, mpmath (mathematical operations)
- numpy, scipy (numerical computing)
- sympy (symbolic math)
- pandas, statistics, statsmodels (data analysis)
- matplotlib, seaborn, plotly (visualization)
- scikit-learn (machine learning)
- rdkit, rdchiral (chemistry, molecular structures, reactions)`
      }
    ];

    if (chatHistory.messages && chatHistory.messages.length > 0) {
      messages.push(...chatHistory.messages.slice(-2));
    }

    messages.push({
      role: 'user',
      content: userQuery
    });

    const llmResponse = await this.callLLM(messages, metadata);
    const content = llmResponse.content;
    
    console.log('Math agent - LLM response received, length:', content.length);
    
    // Extract Python code
    let code = null;
    const codeBlockMatch = content.match(/```python\n([\s\S]*?)```/);
    if (codeBlockMatch) {
      code = codeBlockMatch[1].trim();
    }
    
    // If we found code, execute it
    if (code) {
      try {
        // Ensure imports for common math libraries
        const enhancedCode = this.enhanceMathCode(code);
        
        const executionResult = await this.runCode('python', enhancedCode, metadata.requestId);
        
        if (executionResult.exitCode === 0) {
          // Success - format response with output
          const formattedContent = this.formatMathResponse(content, code, executionResult);
          return {
            content: formattedContent,
            code,
            executionResults: [{
              iteration: 1,
              code: enhancedCode,
              result: executionResult
            }],
            mathExecuted: true,
            metadata: {
              ...metadata,
              agent: 'math-complex',
              tokensUsed: llmResponse.tokensUsed || 0
            }
          };
        }
      } catch (error) {
        console.error('Math execution failed:', error);
      }
    }
    
    // Fallback: return explanation without execution
    return {
      content,
      code: code || 'No code extracted',
      mathExecuted: false,
      metadata: {
        ...metadata,
        agent: 'math-complex',
        tokensUsed: llmResponse.tokensUsed || 0,
        fallback: true
      }
    };
  }

  /**
   * Enhance math code with common imports
   */
  enhanceMathCode(code) {
    // Check what libraries are used and add imports if missing
    const imports = [];
    
    // Always import these for math operations
    if (!code.includes('import math') && !code.includes('from math')) {
      imports.push('import math');
    }
    
    // Check for specific library usage
    if (code.includes('np.') && !code.includes('import numpy')) {
      imports.push('import numpy as np');
    }
    
    if (code.includes('sp.') && !code.includes('import sympy')) {
      imports.push('import sympy as sp');
    }
    
    if (code.includes('plt.') && !code.includes('import matplotlib')) {
      imports.push('import matplotlib.pyplot as plt');
    }
    
    if (code.includes('Decimal') && !code.includes('from decimal')) {
      imports.push('from decimal import Decimal, getcontext');
      imports.push('getcontext().prec = 100  # Set high precision');
    }
    
    if (code.includes('Fraction') && !code.includes('from fractions')) {
      imports.push('from fractions import Fraction');
    }
    
    // Combine imports with original code
    if (imports.length > 0) {
      return imports.join('\n') + '\n\n' + code;
    }
    
    return code;
  }

  /**
   * Execute code via the run_code microservice
   */
  async runCode(language, source, requestId) {
    try {
      console.log(`Math agent - executing ${language} code (${source.length} chars)`);
      console.log('First 200 chars:', source.substring(0, 200));
      
      const response = await fetch(`${this.codeExecutionUrl}/run_code`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          language,
          source,
          timeout: 60000  // 60 seconds for complex mathematical calculations
        })
      });

      const result = await response.json();
      
      console.log('Math agent - execution result:', {
        exitCode: result.exitCode,
        stdoutLength: result.stdout ? result.stdout.length : 0,
        stderrLength: result.stderr ? result.stderr.length : 0
      });
      
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
   * Call LLM for mathematical analysis
   */
  async callLLM(messages, metadata) {
    const response = await fetch(`${this.config.openrouter.baseUrl}/chat/completions`, {
      method: 'POST',
      headers: this.config.openrouter.headers,
      body: JSON.stringify({
        model: this.config.openrouter.model,
        messages,
        temperature: 0.1,  // Low temperature for consistent mathematical reasoning
        max_tokens: 1000000  // Very high limit to capture full responses
      })
    });

    if (!response.ok) {
      throw new Error(`OpenRouter API error: ${response.status}`);
    }

    let data;
    try {
      let text = await response.text();
      if (!text) {
        throw new Error('Empty response from OpenRouter');
      }
      
      // GLM-4.5-air sometimes pads responses with whitespace
      text = text.trim();
      
      // Find the start of JSON (GLM sometimes prepends whitespace)
      const jsonStart = text.indexOf('{');
      if (jsonStart > 0) {
        console.log(`Trimming ${jsonStart} bytes of padding from response`);
        text = text.substring(jsonStart);
      }
      
      data = JSON.parse(text);
    } catch (parseError) {
      console.error('Failed to parse response:', parseError.message);
      console.error('Response text length:', text ? text.length : 0);
      if (text) {
        console.error('First 500 chars:', text.substring(0, 500));
      }
      throw new Error(`JSON parse error: ${parseError.message}`);
    }
    
    if (!data.choices || !data.choices[0]) {
      throw new Error('Invalid response structure');
    }
    
    return {
      content: data.choices[0].message.content,
      tokensUsed: data.usage?.total_tokens || 0
    };
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

// Export singleton instance
module.exports = new MathAgent();