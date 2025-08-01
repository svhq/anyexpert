const config = require('../config');
const logger = require('./utils/logger');
const { SYSTEM_PROMPT } = require('./system-prompt');

/**
 * Math Agent - Handles mathematical computations without numpy
 * Uses only built-in Python libraries available in E2B
 */
class MathAgentNoNumpy {
  constructor() {
    this.config = config;
    this.codeExecutionUrl = config.codeExecution.serviceUrl;
  }

  async handle(userQuery, mathIntent, chatHistory = {}, metadata = {}) {
    try {
      if (!mathIntent.needsExecution) {
        return await this.handleSimpleMath(userQuery, chatHistory, metadata);
      } else {
        return await this.handleComplexMath(userQuery, mathIntent, chatHistory, metadata);
      }
    } catch (error) {
      console.error('Math agent error:', error);
      return await this.handleSimpleMath(userQuery, chatHistory, metadata);
    }
  }

  async handleSimpleMath(userQuery, chatHistory, metadata) {
    const messages = [
      {
        role: 'system',
        content: `${SYSTEM_PROMPT}\n\nSolve mathematical problems step-by-step, showing your work clearly. For simple arithmetic, compute the answer directly.`
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

  async handleComplexMath(userQuery, mathIntent, chatHistory, metadata) {
    const messages = [
      {
        role: 'system',
        content: `${SYSTEM_PROMPT}

For mathematical queries requiring precision or complex calculations, write Python code using ONLY these built-in libraries:
- math - for mathematical functions
- decimal - for high precision decimal arithmetic
- fractions - for exact fraction arithmetic
- statistics - for statistical calculations (mean, median, stdev, etc.)

DO NOT use numpy, scipy, sympy, or any external libraries as they are not available.

Always:
1. Explain the mathematical approach
2. Write clean Python code using only built-in libraries
3. Show the exact result
4. Provide interpretation of the result`
      }
    ];

    if (chatHistory.messages && chatHistory.messages.length > 0) {
      messages.push(...chatHistory.messages.slice(-2));
    }

    messages.push({
      role: 'user',
      content: `${userQuery}\n\nProvide the mathematical solution with Python code using only built-in libraries (no numpy/scipy).`
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
    
    if (code) {
      try {
        // Clean the code to ensure no numpy imports
        const cleanedCode = this.cleanCodeForE2B(code);
        
        const executionResult = await this.runCode('python', cleanedCode, metadata.requestId);
        
        if (executionResult.exitCode === 0) {
          const formattedContent = this.formatMathResponse(content, cleanedCode, executionResult);
          return {
            content: formattedContent,
            code: cleanedCode,
            executionResults: [{
              iteration: 1,
              code: cleanedCode,
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
   * Clean code to ensure it works in E2B environment
   */
  cleanCodeForE2B(code) {
    // Remove numpy imports and replace with math equivalents
    let cleaned = code;
    
    // Remove numpy imports
    cleaned = cleaned.replace(/import numpy as np\n?/g, '');
    cleaned = cleaned.replace(/from numpy import.*\n?/g, '');
    cleaned = cleaned.replace(/import numpy\n?/g, '');
    
    // Replace common numpy functions with math equivalents
    cleaned = cleaned.replace(/np\.pi/g, 'math.pi');
    cleaned = cleaned.replace(/np\.e/g, 'math.e');
    cleaned = cleaned.replace(/np\.sqrt/g, 'math.sqrt');
    cleaned = cleaned.replace(/np\.exp/g, 'math.exp');
    cleaned = cleaned.replace(/np\.log/g, 'math.log');
    cleaned = cleaned.replace(/np\.sin/g, 'math.sin');
    cleaned = cleaned.replace(/np\.cos/g, 'math.cos');
    cleaned = cleaned.replace(/np\.tan/g, 'math.tan');
    cleaned = cleaned.replace(/np\.power/g, 'math.pow');
    cleaned = cleaned.replace(/np\.abs/g, 'abs');
    
    // Ensure math is imported if needed
    if (cleaned.includes('math.') && !cleaned.includes('import math')) {
      cleaned = 'import math\n' + cleaned;
    }
    
    return cleaned;
  }

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
          timeout: 60000
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

  formatMathResponse(originalContent, code, executionResult) {
    const parts = originalContent.split('```python');
    let explanation = parts[0].trim();
    
    let formatted = explanation + '\n\n';
    formatted += '**Calculation:**\n```python\n' + code + '\n```\n\n';
    
    if (executionResult.stdout) {
      formatted += '**Result:**\n```\n' + executionResult.stdout.trim() + '\n```';
      
      const output = executionResult.stdout.trim();
      const lastLine = output.split('\n').pop();
      
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

module.exports = new MathAgentNoNumpy();