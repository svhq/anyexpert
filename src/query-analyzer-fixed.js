const config = require('../config');
const logger = require('./utils/logger');
const { SYSTEM_PROMPT } = require('./system-prompt');

/**
 * Query Analyzer Module - Fixed version with better math detection
 * Determines if a query needs web search or can be answered directly
 */
class QueryAnalyzer {
  constructor() {
    this.config = config;
  }

  /**
   * Detect mathematical intent in user query - IMPROVED VERSION
   * @param {string} userQuery - The user's question
   * @returns {Object} - { isMath: boolean, complexity: string, needsExecution: boolean }
   */
  detectMathIntent(userQuery) {
    const queryLower = userQuery.toLowerCase();
    
    // Skip if it looks like a multiple choice question
    if (queryLower.includes('which of the following') || 
        queryLower.includes('select the correct answer') ||
        /\b[A-E]\)\s/.test(userQuery)) {
      return { isMath: false, complexity: 'none', needsExecution: false };
    }
    
    // Skip if it's asking about concepts rather than calculations
    if (queryLower.includes('statement') && queryLower.includes('true') ||
        queryLower.includes('which') && queryLower.includes('is') ||
        queryLower.includes('definition') ||
        queryLower.includes('property') ||
        queryLower.includes('theorem')) {
      return { isMath: false, complexity: 'none', needsExecution: false };
    }
    
    // Simple math patterns (can be answered directly)
    const simplePatterns = [
      /^\s*\d+\s*[\+\-\*\/\%]\s*\d+\s*$/,  // Simple arithmetic
      /what is \d+ (plus|minus|times|divided by) \d+/i,
      /^[\d\s\+\-\*\/\(\)\.]+$/,  // Pure arithmetic expression
      /\b(add|subtract|multiply|divide)\s+\d+\s+(and|by|from|to)\s+\d+/i
    ];
    
    // Complex math patterns (needs execution) - MORE SPECIFIC
    const complexPatterns = [
      /\b(factorial of \d+|permutation|combination)\b/i,
      /\b(integrate|derivative|differentiate)\s+.*\b(function|equation)\b/i,
      /\b(solve)\s+.*\b(equation|for [a-z])\b/i,
      /\bcalculate.*\b(integral|derivative|limit)\b/i,
      /\b(eigenvalue|eigenvector|determinant)\s+of\s+matrix/i,
      /\b(taylor series|fourier|laplace transform)\s+of/i,
      /\bcalculate.*\b(standard deviation|variance|correlation)\b/i,
      /\bto \d+ (decimal places|digits|significant figures)\b/i,
      /\b(plot|graph|visualize).*\b(function|equation|y\s*=)/i,
      /\b(monte carlo|simulation|numerical)\s+(method|calculation)/i
    ];
    
    // Check if it's simple arithmetic
    const isSimple = simplePatterns.some(pattern => pattern.test(queryLower));
    
    // Check if it's complex math WITH specific calculation request
    const isComplex = complexPatterns.some(pattern => pattern.test(queryLower));
    
    // More restrictive math keyword check
    const calculationKeywords = ['calculate', 'compute', 'evaluate', 'find the value'];
    const hasCalculationRequest = calculationKeywords.some(keyword => 
      queryLower.includes(keyword) && /\d/.test(queryLower)
    );
    
    // Determine if it's math and complexity
    const isMath = isSimple || (isComplex && hasCalculationRequest);
    let complexity = 'none';
    let needsExecution = false;
    
    if (isMath) {
      if (isSimple) {
        complexity = 'simple';
        needsExecution = false;
      } else {
        complexity = 'complex';
        needsExecution = true;
      }
    }
    
    return {
      isMath,
      complexity,
      needsExecution
    };
  }

  /**
   * Detect coding intent in user query
   * @param {string} userQuery - The user's question
   * @returns {Object} - { isCoding: boolean, language: string|null, type: string }
   */
  detectCodingIntent(userQuery) {
    const queryLower = userQuery.toLowerCase();
    
    // Skip if it looks like a multiple choice question
    if (queryLower.includes('which of the following') || 
        queryLower.includes('select the correct answer') ||
        /\b[A-E]\)\s/.test(userQuery)) {
      return { isCoding: false, language: null, type: 'none' };
    }
    
    // Language detection patterns
    const languages = {
      python: /\b(python|py|pip|pandas|numpy|matplotlib|jupyter|\.py\b)/i,
      javascript: /\b(javascript|js|node|npm|react|vue|angular|\.js\b)/i,
      java: /\b(java|\.java\b|spring|maven|gradle)/i,
      cpp: /\b(c\+\+|cpp|\.cpp\b|\.h\b)/i,
      csharp: /\b(c#|csharp|\.cs\b|\.net|dotnet)/i,
      go: /\b(golang|go\b|\.go\b)/i,
      rust: /\b(rust|cargo|\.rs\b)/i,
      php: /\b(php|\.php\b|laravel|composer)/i,
      ruby: /\b(ruby|rails|gem|\.rb\b)/i,
      bash: /\b(bash|shell|terminal|\.sh\b|command line)/i
    };
    
    // Coding activity patterns - more specific
    const codingPatterns = [
      /\b(debug|debugging this code|error in my code|traceback|stack trace)\b/i,
      /\b(output of this code|result of this program|run this|execute this code)\b/i,
      /\b(fix this bug|optimize this code|refactor this|code review)\b/i,
      /\bwrite a (function|method|class|program) to\b/i,
      /\b(syntax error|runtime error|compile error|build failed)\b/i,
      /```[\s\S]*```/,  // Code blocks
      /\bprint\s*\(|\bconsole\.log\(|\bSystem\.out/i,  // Print statements
    ];
    
    // Check for language indicators
    let detectedLanguage = null;
    for (const [lang, pattern] of Object.entries(languages)) {
      if (pattern.test(userQuery)) {
        detectedLanguage = lang;
        break;
      }
    }
    
    // Check for coding patterns
    const hasCodingPattern = codingPatterns.some(pattern => pattern.test(userQuery));
    
    // Determine coding type
    let codingType = 'none';
    if (hasCodingPattern || detectedLanguage) {
      if (userQuery.includes('debug') || userQuery.includes('error') || userQuery.includes('fix')) {
        codingType = 'code-debug';
      } else if (userQuery.includes('explain') || userQuery.includes('what does') || userQuery.includes('how does')) {
        codingType = 'code-explain';
      } else if (userQuery.includes('write') || userQuery.includes('create') || userQuery.includes('implement')) {
        codingType = 'code-write';
      } else if (/```[\s\S]*```/.test(userQuery) || userQuery.includes('output') || userQuery.includes('result')) {
        codingType = 'code-exec';
      } else {
        codingType = 'code-general';
      }
    }
    
    return {
      isCoding: !!(hasCodingPattern || detectedLanguage),
      language: detectedLanguage,
      type: codingType
    };
  }

  /**
   * Assess if a query needs web search or code execution
   * @param {string} userQuery - The user's question
   * @param {Object} chatHistory - Previous conversation context
   * @returns {Promise<Object>} - { needsSearch: boolean, rationale: string, codingIntent?: Object }
   */
  async assess(userQuery, chatHistory = {}) {
    try {
      // First check for math intent
      const mathIntent = this.detectMathIntent(userQuery);
      
      // If it's math-related, return early with math decision
      if (mathIntent.isMath) {
        const result = {
          needsSearch: false,  // Math queries don't need web search
          rationale: `Detected ${mathIntent.complexity} math query${mathIntent.needsExecution ? ' requiring execution' : ''}`,
          mathIntent
        };
        
        logger.logSearchDecision('query-analyzer', userQuery, false, result.rationale, 'math-detection');
        return result;
      }
      
      // Then check for coding intent
      const codingIntent = this.detectCodingIntent(userQuery);
      
      // If it's coding-related, return early with coding decision
      if (codingIntent.isCoding) {
        const result = {
          needsSearch: false,  // Code queries don't need web search
          rationale: `Detected ${codingIntent.type} query for ${codingIntent.language || 'code'}`,
          codingIntent
        };
        
        logger.logSearchDecision('query-analyzer', userQuery, false, result.rationale, 'coding-detection');
        return result;
      }

      // For general queries, use LLM assessment
      const messages = [
        {
          role: 'system',
          content: `${SYSTEM_PROMPT}

You are assessing whether a user query requires web search or can be answered from general knowledge.

Return a JSON object with:
{
  "needsSearch": "Y" or "N",
  "rationale": "Brief explanation of why search is/isn't needed"
}

Guidelines:
- "Y" if the query needs current information, specific facts, recent events, or specialized knowledge
- "N" if the query can be answered from general knowledge, is conversational, or is asking for explanations of well-known concepts`
        },
        { role: 'user', content: userQuery }
      ];

      const response = await fetch(`${this.config.openrouter.baseUrl}/chat/completions`, {
        method: 'POST',
        headers: this.config.openrouter.headers,
        body: JSON.stringify({
          model: this.config.openrouter.model,
          messages,
          temperature: 0,
          max_tokens: 200
        })
      });

      if (!response.ok) {
        throw new Error(`OpenRouter API error: ${response.status}`);
      }

      const data = await response.json();
      const content = data.choices[0].message.content;
      
      // Parse the JSON response
      try {
        const analysis = JSON.parse(content);
        const result = {
          needsSearch: analysis.needsSearch === 'Y',
          rationale: analysis.rationale
        };
        
        logger.logSearchDecision('query-analyzer', userQuery, result.needsSearch, result.rationale, 'json');
        return result;
      } catch (parseError) {
        // Fallback: try to extract from non-JSON response
        const jsonMatch = content.match(/\{[\s\S]*\}/);
        if (jsonMatch) {
          try {
            const extractedJson = JSON.parse(jsonMatch[0]);
            return {
              needsSearch: extractedJson.needsSearch === 'Y',
              rationale: extractedJson.rationale || 'Assessment completed'
            };
          } catch (e) {
            // If JSON parsing fails, use text analysis
          }
        }
        
        // Ultimate fallback: analyze the text response
        const needsSearch = content.toLowerCase().includes('"y"') || 
                          content.toLowerCase().includes('yes') ||
                          content.toLowerCase().includes('search needed');
        
        logger.logSearchDecision('query-analyzer', userQuery, needsSearch, 'Parsed from text response', 'text');
        
        return {
          needsSearch,
          rationale: content.substring(0, 200)
        };
      }
    } catch (error) {
      console.error('Query analysis error:', error);
      // Default to no search on error
      return {
        needsSearch: false,
        rationale: 'Analysis error - defaulting to direct answer'
      };
    }
  }
}

module.exports = new QueryAnalyzer();