const config = require('../config');
const logger = require('./utils/logger');
const FlexibleJSONParser = require('./utils/json-parser');

/**
 * Search Planner - Generates optimized search queries for each round
 */
class SearchPlanner {
  constructor() {
    this.config = config;
    this.queryTemplates = this.initializeTemplates();
  }

  /**
   * Generate search queries based on round and previous results
   * @param {string} userQuery - Original user query
   * @param {string} rationale - Why search is needed
   * @param {number} round - Current search round (0-based)
   * @param {Array} previousResults - Results from previous rounds
   * @returns {Promise<Array>} - Array of search queries
   */
  async generate(userQuery, rationale, round = 0, previousResults = []) {
    try {
      const prompt = this.buildPlanningPrompt(userQuery, rationale, round, previousResults);
      
      const messages = [
        {
          role: 'system',
          content: 'You are a search query optimizer that generates diverse, high-quality search queries to find comprehensive information.'
        },
        {
          role: 'user',
          content: prompt
        }
      ];

      const response = await fetch(`${this.config.openrouter.baseUrl}/chat/completions`, {
        method: 'POST',
        headers: this.config.openrouter.headers,
        body: JSON.stringify({
          model: this.config.openrouter.model,
          messages,
          temperature: 0.3 // Some creativity but mostly consistent
        })
      });

      if (!response.ok) {
        throw new Error(`OpenRouter API error: ${response.status}`);
      }

      const data = await response.json();
      const content = data.choices[0].message.content;
      
      // Use flexible parser
      const parseResult = FlexibleJSONParser.parse(content);
      
      if (parseResult.success && parseResult.data) {
        const result = parseResult.data;
        // Ensure we return an array of queries
        return Array.isArray(result.queries) ? result.queries : [userQuery];
      } else {
        logger.logJsonParsingFailure('search-planner', 'search-planner', content, parseResult.error);
        // Don't log warning in production - fallback is working as designed
        logger.debug('Using fallback query extraction');
        
        // Try to extract queries from response text
        const queries = this.extractQueriesFromText(content, userQuery, round);
        if (queries.length > 0) {
          return queries;
        }
        
        // Use intelligent fallback
        return this.generateFallbackQueries(userQuery, round);
      }
    } catch (error) {
      console.error('Search planning error:', error);
      // Fallback to simple query variations
      return this.generateFallbackQueries(userQuery, round);
    }
  }

  /**
   * Build the planning prompt based on round
   */
  buildPlanningPrompt(userQuery, rationale, round, previousResults) {
    let prompt = `Generate search queries for: "${userQuery}"\n`;
    prompt += `Search needed because: ${rationale}\n\n`;

    if (round === 0) {
      prompt += `This is the FIRST search round. Generate 3-5 diverse queries that:
1. Include the original query (slightly rephrased)
2. Add a more specific/technical version
3. Add a broader/context version
4. Include temporal qualifiers if relevant (2024, 2025, latest, current)
5. Use different phrasings to maximize result diversity\n`;
    } else {
      prompt += `This is search round ${round + 1}. Previous searches found ${previousResults.length} results.\n`;
      
      // Extract what we already know
      const knownDomains = [...new Set(previousResults.map(r => new URL(r.url).hostname))];
      prompt += `Already searched domains: ${knownDomains.slice(0, 5).join(', ')}\n\n`;
      
      prompt += `Generate 3-4 NEW queries that:
1. Fill information gaps not covered in previous results
2. Seek different perspectives or sources
3. Validate or expand on partial information found
4. Use MORE SPECIFIC terms based on what was discovered\n`;
      
      if (round >= 2) {
        prompt += `5. Focus on fact-checking and authoritative sources\n`;
      }
    }

    prompt += `\nOutput format:
{
  "queries": [
    "search query 1",
    "search query 2",
    "search query 3"
  ]
}`;

    return prompt;
  }

  /**
   * Extract queries from LLM response text when JSON parsing fails
   */
  extractQueriesFromText(content, userQuery, round) {
    const queries = [];
    
    // Look for quoted strings that might be queries
    const quotedMatches = content.match(/"([^"]+)"/g);
    if (quotedMatches) {
      for (const match of quotedMatches) {
        const query = match.replace(/"/g, '').trim();
        if (query.length > 10 && query.length < 200) {
          queries.push(query);
        }
      }
    }
    
    // Look for numbered list items
    const numberedMatches = content.match(/\d+\.\s*([^\n]+)/g);
    if (numberedMatches) {
      for (const match of numberedMatches) {
        const query = match.replace(/^\d+\.\s*/, '').trim();
        if (query.length > 10 && query.length < 200) {
          queries.push(query);
        }
      }
    }
    
    // Look for bullet points
    const bulletMatches = content.match(/[-*]\s*([^\n]+)/g);
    if (bulletMatches) {
      for (const match of bulletMatches) {
        const query = match.replace(/^[-*]\s*/, '').trim();
        if (query.length > 10 && query.length < 200) {
          queries.push(query);
        }
      }
    }
    
    // Deduplicate and limit
    const uniqueQueries = [...new Set(queries)];
    return uniqueQueries.slice(0, 4);
  }

  /**
   * Fallback query generation if LLM fails
   */
  generateFallbackQueries(userQuery, round) {
    const queries = [userQuery];
    
    if (round === 0) {
      // Add variations
      queries.push(`${userQuery} explained`);
      queries.push(`latest ${userQuery} 2025`);
      queries.push(`${userQuery} comprehensive guide`);
    } else {
      // Add refinements
      queries.push(`${userQuery} detailed analysis`);
      queries.push(`${userQuery} expert opinion`);
      queries.push(`${userQuery} authoritative sources`);
    }
    
    return queries.slice(0, 4); // Limit to 4 queries
  }

  /**
   * Initialize query templates for different domains
   * (Will be expanded in Phase 4)
   */
  initializeTemplates() {
    return {
      technical: {
        patterns: [
          '{topic} implementation best practices',
          '{topic} technical documentation',
          '{topic} tutorial guide 2025',
          '{topic} common issues solutions'
        ]
      },
      research: {
        patterns: [
          '{topic} latest research 2025',
          '{topic} academic papers',
          '{topic} scientific studies',
          '{topic} peer reviewed'
        ]
      },
      news: {
        patterns: [
          '{topic} latest news today',
          '{topic} breaking updates',
          '{topic} current developments',
          '{topic} recent announcements'
        ]
      },
      general: {
        patterns: [
          'what is {topic}',
          '{topic} explained simply',
          '{topic} comprehensive overview',
          '{topic} key facts'
        ]
      }
    };
  }
}

// Export singleton instance
module.exports = new SearchPlanner();