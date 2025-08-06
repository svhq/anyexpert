const config = require('../config');

class WebSearch {
  constructor() {
    this.apiKey = config.serper.apiKey;
    this.baseUrl = config.serper.baseUrl;
  }

  async search(query, options = {}) {
    const searchOptions = {
      q: query,
      gl: options.country || 'us',
      hl: options.language || 'en',
      num: options.numResults || 10,
      ...options
    };

    try {
      const response = await fetch(`${this.baseUrl}/search`, {
        method: 'POST',
        headers: {
          'X-API-KEY': this.apiKey,
          'Content-Type': 'application/json'
        },
        body: JSON.stringify(searchOptions)
      });

      if (!response.ok) {
        throw new Error(`Serper API error: ${response.status}`);
      }

      const data = await response.json();
      
      // Log raw response for debugging if needed
      if (process.env.NODE_ENV === 'development') {
        console.log('Serper raw response structure:', {
          hasOrganic: !!data.organic,
          organicCount: data.organic?.length || 0,
          hasAnswerBox: !!data.answerBox,
          hasKnowledgeGraph: !!data.knowledgeGraph
        });
      }
      
      return this.formatResults(data);
    } catch (error) {
      console.error('Web search error:', error);
      throw error;
    }
  }

  formatResults(data) {
    const results = {
      query: data.searchParameters?.q || '',
      organic: [],
      answerBox: null,
      knowledgeGraph: null,
      relatedSearches: []
    };

    // Organic search results
    if (data.organic) {
      results.organic = data.organic.map(item => ({
        title: item.title || 'No title',
        link: item.link,
        url: item.link, // Ensure both link and url are available
        snippet: item.snippet || 'No snippet available',
        position: item.position
      }));
    }

    // Answer box if available
    if (data.answerBox) {
      results.answerBox = {
        answer: data.answerBox.answer,
        snippet: data.answerBox.snippet,
        title: data.answerBox.title,
        link: data.answerBox.link
      };
    }

    // Knowledge graph if available
    if (data.knowledgeGraph) {
      results.knowledgeGraph = {
        title: data.knowledgeGraph.title,
        type: data.knowledgeGraph.type,
        description: data.knowledgeGraph.description,
        attributes: data.knowledgeGraph.attributes
      };
    }

    // Related searches
    if (data.relatedSearches) {
      results.relatedSearches = data.relatedSearches.map(item => item.query);
    }

    return results;
  }

  async searchNews(query, options = {}) {
    const newsOptions = {
      q: query,
      type: 'news',
      ...options
    };

    return this.search(query, newsOptions);
  }

  async searchImages(query, options = {}) {
    const imageOptions = {
      q: query,
      type: 'images',
      ...options
    };

    return this.search(query, imageOptions);
  }

  /**
   * Run multiple searches in parallel
   * @param {Array<string>} queries - Array of search queries
   * @param {Object} options - Search options
   * @returns {Promise<Array>} - Array of search results
   */
  async runBatch(queries, options = {}) {
    try {
      // Execute searches in parallel
      const searchPromises = queries.map(query => 
        this.search(query, options)
          .then(results => ({
            query,
            success: true,
            results: results.organic || []
          }))
          .catch(error => ({
            query,
            success: false,
            error: error.message,
            results: []
          }))
      );

      const batchResults = await Promise.all(searchPromises);
      
      // Flatten results while maintaining source query info
      const flatResults = [];
      for (const batch of batchResults) {
        if (batch.success && batch.results.length > 0) {
          // Add source query to each result and ensure URL consistency
          const resultsWithQuery = batch.results.map(result => ({
            ...result,
            url: result.url || result.link, // Ensure url field is always available
            sourceQuery: batch.query
          }));
          flatResults.push(...resultsWithQuery);
        }
      }

      return flatResults;
    } catch (error) {
      console.error('Batch search error:', error);
      throw error;
    }
  }

  /**
   * Scrape full content from a URL using Serper's scrape endpoint
   * @param {string} url - URL to scrape
   * @param {Object} options - Additional options
   * @returns {Promise<Object>} - Scraped content
   */
  async scrape(url, options = {}) {
    try {
      const response = await fetch(`${this.baseUrl}/scrape`, {
        method: 'POST',
        headers: {
          'X-API-KEY': this.apiKey,
          'Content-Type': 'application/json'
        },
        body: JSON.stringify({ url })
      });

      if (!response.ok) {
        const errorText = await response.text();
        throw new Error(`Serper scrape error: ${response.status} - ${errorText}`);
      }

      const data = await response.json();
      
      // Log for debugging if needed
      if (process.env.NODE_ENV === 'development') {
        console.log('Serper scrape response:', {
          hasText: !!data.text,
          hasMarkdown: !!data.markdown,
          textLength: data.text?.length || 0
        });
      }
      
      return {
        url,
        title: data.title || 'No title',
        content: data.markdown || data.text || '',
        metadata: {
          description: data.description,
          language: data.language,
          ...data.metadata
        }
      };
    } catch (error) {
      console.error('Scraping error for URL:', url, error);
      // Return null to indicate failure, let agent handle gracefully
      return null;
    }
  }
}

// Export singleton instance
module.exports = new WebSearch();