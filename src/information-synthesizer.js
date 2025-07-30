const config = require('../config');
const { SYSTEM_PROMPT } = require('./system-prompt');

/**
 * Information Synthesizer - Combines search results with expert knowledge
 */
class InformationSynthesizer {
  constructor() {
    this.config = config;
  }

  /**
   * Generate direct answer without search
   */
  async answerDirect(userQuery, chatHistory, metadata) {
    try {
      const messages = [
        {
          role: 'system',
          content: SYSTEM_PROMPT
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

      const response = await fetch(`${this.config.openrouter.baseUrl}/chat/completions`, {
        method: 'POST',
        headers: this.config.openrouter.headers,
        body: JSON.stringify({
          model: this.config.openrouter.model,
          messages,
          temperature: 0.7
        })
      });

      if (!response.ok) {
        throw new Error(`OpenRouter API error: ${response.status}`);
      }

      const data = await response.json();
      
      const content = data.choices[0].message.content;
      
      // Log response quality
      const logger = require('./utils/logger');
      logger.logResponseQuality(
        metadata?.requestId || 'unknown',
        content.length,
        false, // Direct answers aren't truncated since we removed token limits
        0 // No citations for direct answers
      );
      
      return {
        content,
        searchPerformed: false,
        metadata: {
          ...metadata,
          model: this.config.openrouter.model,
          tokensUsed: data.usage?.total_tokens || 0
        }
      };
    } catch (error) {
      console.error('Direct answer generation error:', error);
      throw error;
    }
  }

  /**
   * Compose answer using search results
   */
  async compose(userQuery, searchResults, chatHistory, metadata) {
    try {
      // Prepare search evidence
      const searchEvidence = this.formatSearchResults(searchResults);
      
      const messages = [
        {
          role: 'system',
          content: `${SYSTEM_PROMPT}\n\nYou have access to web search results. When using information from these sources, cite them using [1], [2], etc. format.`
        }
      ];

      // Add relevant chat history
      if (chatHistory.messages && chatHistory.messages.length > 0) {
        messages.push(...chatHistory.messages.slice(-2)); // Last exchange
      }

      // Add search results and query
      messages.push({
        role: 'user',
        content: `Based on the following search results, please answer this question as the most appropriate expert: "${userQuery}"

Search Results:
${searchEvidence.text}

Please provide a comprehensive answer using the search results and your expertise. Include citations [n] for specific facts from the sources.`
      });

      const response = await fetch(`${this.config.openrouter.baseUrl}/chat/completions`, {
        method: 'POST',
        headers: this.config.openrouter.headers,
        body: JSON.stringify({
          model: this.config.openrouter.model,
          messages,
          temperature: 0.7
        })
      });

      if (!response.ok) {
        throw new Error(`OpenRouter API error: ${response.status}`);
      }

      const data = await response.json();
      const content = data.choices[0].message.content;
      
      // Log response quality for search-based answers
      const logger = require('./utils/logger');
      logger.logResponseQuality(
        metadata?.requestId || 'unknown',
        content.length,
        false, // No truncation with removed token limits
        searchEvidence.sources.length
      );
      
      return {
        content,
        searchPerformed: true,
        sources: searchEvidence.sources,
        metadata: {
          ...metadata,
          model: this.config.openrouter.model,
          tokensUsed: data.usage?.total_tokens || 0,
          sourcesUsed: searchEvidence.sources.length
        }
      };
    } catch (error) {
      console.error('Response composition error:', error);
      throw error;
    }
  }

  /**
   * Format search results for LLM consumption
   */
  formatSearchResults(searchResults) {
    if (!searchResults || searchResults.length === 0) {
      return { text: 'No search results available.', sources: [] };
    }

    const sources = [];
    const texts = [];
    
    // Deduplicate by URL
    const uniqueResults = [];
    const seenUrls = new Set();
    
    for (const result of searchResults) {
      const url = result.url || result.link;
      if (!seenUrls.has(url)) {
        seenUrls.add(url);
        uniqueResults.push(result);
      }
    }

    // Format top results (limit to avoid context overflow)
    const topResults = uniqueResults.slice(0, 10);
    
    topResults.forEach((result, index) => {
      const num = index + 1;
      const source = {
        number: num,
        title: result.title,
        url: result.url || result.link,
        snippet: result.snippet
      };
      sources.push(source);
      
      texts.push(`[${num}] ${source.title}
URL: ${source.url}
${source.snippet}
`);
    });

    return {
      text: texts.join('\n---\n'),
      sources
    };
  }
}

// Export singleton instance
module.exports = new InformationSynthesizer();