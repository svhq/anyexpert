const config = require('../config');
const logger = require('./utils/logger');

/**
 * Query Analyzer - Determines if web search is needed for a query
 * Uses LLM with few-shot examples for consistent decision making
 */
class QueryAnalyzer {
  constructor() {
    this.config = config;
    this.analysisPrompt = this.buildAnalysisPrompt();
  }

  /**
   * Detect math intent in user query
   * @param {string} userQuery - The user's question
   * @returns {Object} - { isMath: boolean, complexity: string, needsExecution: boolean }
   */
  detectMathIntent(userQuery) {
    const queryLower = userQuery.toLowerCase();
    
    // Simple arithmetic patterns (LLM can handle directly)
    const simplePatterns = [
      /^\s*\d+\s*[\+\-\*\/]\s*\d+\s*$/,  // Simple 2-operand arithmetic
      /^what is \d+ (plus|minus|times|divided by) \d+/i,
      /^calculate \d+\s*[\+\-\*\/]\s*\d+$/i,
      /^compute \d+\s*[\+\-\*\/]\s*\d+$/i,
      /^\d+\s*[\+\-\*\/]\s*\d+\s*[\+\-\*\/]\s*\d+$/,  // 3-operand simple math
      /^what is \d+\s*[\+\-\*\/]\s*\d+$/i
    ];
    
    // Complex math patterns (needs execution)
    const complexPatterns = [
      /\b(factorial|permutation|combination)\b/i,
      /\b(integral|integrate|derivative|differentiate|∫|∂)\b/i,
      /\b(matrix|matrices|eigenvalue|eigenvector|determinant)\b/i,
      /\b(solve.*equation|solve for [a-z])\b/i,
      /\b(limit|lim\s)/i,
      /\b(sum|product|∑|∏)\s*(from|of)/i,
      /\b(taylor series|fourier|laplace transform)\b/i,
      /\b(standard deviation|variance|correlation|regression)\b/i,
      /\b(probability|distribution|statistical)\b/i,
      /\bto \d+ (decimal places|digits|significant figures)\b/i,
      /\b(exact|precise|high precision)\b/i,
      /\b(plot|graph|visualize|chart)\b.*\b(function|equation|data)\b/i,
      /\b(monte carlo|simulation|numerical)\b/i,
      /\b(differential equation|ODE|PDE)\b/i
    ];
    
    // Math operation keywords
    const mathKeywords = [
      'calculate', 'compute', 'evaluate', 'simplify', 'factor',
      'expand', 'solve', 'find', 'determine', 'derive',
      'prove', 'verify', 'integrate', 'differentiate'
    ];
    
    // Check if it's simple arithmetic
    const isSimple = simplePatterns.some(pattern => pattern.test(queryLower));
    
    // Check if it's complex math
    const isComplex = complexPatterns.some(pattern => pattern.test(queryLower));
    
    // Check for math keywords
    const hasMathKeyword = mathKeywords.some(keyword => 
      queryLower.includes(keyword) && 
      (queryLower.includes('number') || queryLower.includes('value') || 
       queryLower.includes('result') || /\d/.test(queryLower))
    );
    
    // Determine if it's math and complexity
    const isMath = isSimple || isComplex || hasMathKeyword;
    let complexity = 'none';
    let needsExecution = false;
    
    if (isMath) {
      if (isSimple && !queryLower.includes('exact') && !queryLower.includes('precise')) {
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
    
    // Coding activity patterns
    const codingPatterns = [
      /\b(debug|debugging|error|traceback|stack trace|exception)\b/i,
      /\b(output of|result of|run this|execute|what does.*return)\b/i,
      /\b(fix.*bug|optimize|refactor|code review)\b/i,
      /\b(algorithm|function|method|class|variable)\b/i,
      /\b(syntax error|runtime error|compile|build)\b/i,
      /```[\s\S]*```/,  // Code blocks
      /\bprint\s*\(|\bconsole\.log\(|\bSystem\.out/i,  // Print statements
      /\bdef\s+\w+|\bfunction\s+\w+|\bclass\s+\w+/i,  // Function/class definitions
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
    let codingType = null;
    if (hasCodingPattern || detectedLanguage) {
      if (/\b(output|result|return|execute|run)\b/i.test(queryLower)) {
        codingType = 'code-exec';  // Needs execution
      } else if (/\b(explain|how|what is|difference|compare)\b/i.test(queryLower)) {
        codingType = 'code-explain';  // Needs explanation only
      } else {
        codingType = 'code-exec';  // Default to execution for ambiguous cases
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
      
      // For non-coding queries, proceed with normal search analysis
      const messages = [
        {
          role: 'system',
          content: this.analysisPrompt
        },
        {
          role: 'user',
          content: `Query: "${userQuery}"\n\nContext: ${JSON.stringify(chatHistory.recentTopics || [])}\n\nAnalyze if this query needs web search.`
        }
      ];

      const response = await fetch(`${this.config.openrouter.baseUrl}/chat/completions`, {
        method: 'POST',
        headers: this.config.openrouter.headers,
        body: JSON.stringify({
          model: this.config.openrouter.model,
          messages,
          temperature: 0 // Deterministic for consistent decisions
        })
      });

      if (!response.ok) {
        throw new Error(`OpenRouter API error: ${response.status}`);
      }

      const data = await response.json();
      const content = data.choices[0].message.content;
      
      try {
        const analysis = JSON.parse(content);
        const result = {
          needsSearch: analysis.needsSearch === 'Y',
          rationale: analysis.rationale || 'No rationale provided'
        };
        logger.logSearchDecision('query-analyzer', userQuery, result.needsSearch, result.rationale, 'json');
        return result;
      } catch (parseError) {
        logger.logJsonParsingFailure('query-analyzer', 'query-analyzer', content, parseError);
        // Fallback: try to extract from text if JSON parsing fails
        console.warn('Failed to parse JSON response, using enhanced fallback logic');
        
        // Try to extract JSON from within the response
        const jsonMatch = content.match(/\{[^{}]*"needsSearch"[^{}]*\}/i);
        if (jsonMatch) {
          try {
            const extractedJson = JSON.parse(jsonMatch[0]);
            return {
              needsSearch: extractedJson.needsSearch === 'Y',
              rationale: extractedJson.rationale || 'Extracted from partial JSON'
            };
          } catch (e) {
            // Continue to heuristic analysis
          }
        }
        
        // Enhanced heuristic analysis
        const responseIndicatesSearch = 
          content.toLowerCase().includes('"needssearch": "y"') || 
          content.toLowerCase().includes('needssearch: y') ||
          content.toLowerCase().includes('need search') ||
          content.toLowerCase().includes('requires search') ||
          content.toLowerCase().includes('current information needed') ||
          content.toLowerCase().includes('web search') ||
          content.toLowerCase().includes('external information');
          
        // Enhanced query content analysis for temporal/current info needs
        const queryLower = userQuery.toLowerCase();
        const queryNeedsSearch = 
          queryLower.includes('latest') ||
          queryLower.includes('current') ||
          queryLower.includes('2025') ||
          queryLower.includes('2024') ||
          queryLower.includes('today') ||
          queryLower.includes('recent') ||
          queryLower.includes('news') ||
          queryLower.includes('update') ||
          queryLower.includes('breakthrough') ||
          queryLower.includes('development') ||
          queryLower.includes('trend') ||
          queryLower.includes('regulation') ||
          queryLower.includes('guideline') ||
          queryLower.includes('best practice') ||
          queryLower.includes('best practices') ||
          queryLower.includes('modern') ||
          queryLower.includes('now') ||
          queryLower.includes('this year') ||
          queryLower.includes('how to implement') ||
          queryLower.includes('implementation') ||
          queryLower.includes('tutorial') ||
          queryLower.includes('guide') ||
          queryLower.includes('setup') ||
          queryLower.includes('configure') ||
          queryLower.includes('optimization') ||
          queryLower.includes('performance') ||
          queryLower.includes('troubleshooting') ||
          queryLower.includes('solutions') ||
          queryLower.includes('approaches') ||
          queryLower.includes('strategies') ||
          queryLower.includes('techniques') ||
          queryLower.includes('methods') ||
          queryLower.includes('frameworks') ||
          queryLower.includes('libraries') ||
          queryLower.includes('tools') ||
          queryLower.includes('recommendations');
          
        const needsSearch = responseIndicatesSearch || queryNeedsSearch;
        
        const result = {
          needsSearch,
          rationale: needsSearch ? 'Query contains temporal/current info request (heuristic analysis)' : 'General knowledge query (heuristic analysis)'
        };
        
        logger.logSearchDecision('query-analyzer', userQuery, result.needsSearch, result.rationale, 'heuristic');
        return result;
      }
    } catch (error) {
      console.error('Query analysis error:', error);
      // Default to searching on error
      return {
        needsSearch: true,
        rationale: 'Error in analysis, defaulting to search'
      };
    }
  }

  /**
   * Build the analysis prompt with few-shot examples
   */
  buildAnalysisPrompt() {
    return `You are a query analyzer that determines if a user's question requires web search for current/external information.

IMPORTANT: You MUST respond with ONLY a JSON object in this exact format:
{"needsSearch": "Y", "rationale": "explanation"} 
or 
{"needsSearch": "N", "rationale": "explanation"}

Few-shot examples:

Query: "What is photosynthesis?"
Output: {"needsSearch": "N", "rationale": "Basic scientific concept that can be explained from general knowledge"}

Query: "Latest news about OpenAI?"
Output: {"needsSearch": "Y", "rationale": "Requires current events and recent updates"}

Query: "How to implement quicksort in Python?"
Output: {"needsSearch": "N", "rationale": "Standard algorithm that can be explained from programming knowledge"}

Query: "Current stock price of AAPL"
Output: {"needsSearch": "Y", "rationale": "Requires real-time financial data"}

Query: "Who won the 2024 Super Bowl?"
Output: {"needsSearch": "Y", "rationale": "Specific recent event information needed"}

Query: "Explain the theory of relativity"
Output: {"needsSearch": "N", "rationale": "Established scientific theory from general knowledge"}

Query: "Best restaurants in Tokyo 2025"
Output: {"needsSearch": "Y", "rationale": "Current recommendations and recent reviews needed"}

Decision criteria:
- Need search (Y): Current events, real-time data, recent updates, location-specific current info, latest research/developments
- No search (N): General knowledge, established facts, standard procedures, historical events (pre-2024), basic explanations

Analyze the given query and respond with the JSON format.`;
  }
}

// Export singleton instance
module.exports = new QueryAnalyzer();