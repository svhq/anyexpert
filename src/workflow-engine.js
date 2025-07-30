const queryAnalyzer = require('./query-analyzer');
const searchPlanner = require('./search-planner');
const webSearch = require('./web-search');
const informationSynthesizer = require('./information-synthesizer');
const codeAgent = require('./code-agent-simple');
const mathAgent = require('./math-agent');
const logger = require('./utils/logger');

// Configuration constants
const MAX_ROUNDS = parseInt(process.env.MAX_ROUNDS || '3', 10);
const CONF_THRESHOLD = parseFloat(process.env.CONF_THRESHOLD || '0.8');

/**
 * Main workflow engine that orchestrates the expert system with web search
 */
class WorkflowEngine {
  constructor(config) {
    this.config = config;
    this.maxRounds = MAX_ROUNDS;
    this.confidenceThreshold = CONF_THRESHOLD;
  }

  /**
   * Main entry point for answering user queries
   * @param {string} userQuery - The user's question
   * @param {Object} chatHistory - Previous conversation context
   * @returns {Promise<Object>} - The expert response with citations
   */
  async answer(userQuery, chatHistory = {}) {
    const startTime = Date.now();
    const requestId = this.generateRequestId();
    const logs = [];

    try {
      // Step 1: Analyze query type (search, code, or direct)
      logger.info({ requestId, step: 'query-analysis', query: userQuery });
      const analysis = await queryAnalyzer.assess(userQuery, chatHistory);
      
      // Step 2: Route based on query type
      if (analysis.mathIntent && analysis.mathIntent.isMath) {
        // Handle math queries with math agent
        logger.info({ requestId, queryType: 'math', complexity: analysis.mathIntent.complexity });
        return await mathAgent.handle(
          userQuery,
          analysis.mathIntent,
          chatHistory,
          { requestId }
        );
      } else if (analysis.codingIntent && analysis.codingIntent.isCoding) {
        // Handle coding queries with code agent
        logger.info({ requestId, queryType: 'coding', codingType: analysis.codingIntent.type });
        return await codeAgent.handle(
          userQuery,
          analysis.codingIntent,
          chatHistory,
          { requestId }
        );
      } else if (!analysis.needsSearch) {
        // Handle direct reasoning queries
        logger.info({ requestId, needsSearch: false, reason: analysis.rationale });
        return await informationSynthesizer.answerDirect(
          userQuery, 
          chatHistory,
          { requestId, searchPerformed: false }
        );
      }

      // Step 3: Execute search loop for search-required queries
      let confidence = 0;
      let allSearchResults = [];
      let round = 0;

      while (this.shouldContinue(confidence, round)) {
        logger.info({ 
          requestId, 
          round, 
          previousConfidence: confidence,
          status: 'search-round-start' 
        });

        // Generate search queries
        const queries = await searchPlanner.generate(
          userQuery,
          analysis.rationale,
          round,
          allSearchResults
        );

        // Execute searches in parallel
        const roundResults = await webSearch.runBatch(queries);
        
        // For Phase 1, we'll use simple result compilation
        // Phase 2 will add passage ranking
        const processedResults = this.processSearchResults(roundResults);
        
        // Simple confidence estimation for Phase 1
        // Phase 3 will add expert-selector for proper confidence scoring
        confidence = this.estimateConfidence(processedResults, userQuery);

        // Log round details
        const roundLog = {
          round,
          queries,
          resultCount: processedResults.length,
          topUrls: processedResults.slice(0, 3).map(r => r.url),
          confidence,
          duration: Date.now() - startTime
        };
        logs.push(roundLog);
        logger.info({ requestId, ...roundLog });

        // Accumulate results
        allSearchResults = allSearchResults.concat(processedResults);
        round++;
      }

      // Step 3: Generate final response
      const response = await informationSynthesizer.compose(
        userQuery,
        allSearchResults,
        chatHistory,
        {
          requestId,
          rounds: logs.length,
          finalConfidence: confidence
        }
      );

      // Log final metrics
      const finalMetrics = {
        requestId,
        userId: chatHistory.userId || 'anonymous',
        totalDuration: Date.now() - startTime,
        rounds: logs.length,
        finalConfidence: confidence,
        serperCalls: logs.reduce((sum, log) => sum + log.queries.length, 0),
        tokensUsed: response.tokensUsed || 0
      };
      logger.info({ type: 'request-complete', ...finalMetrics });

      return {
        ...response,
        metadata: {
          ...finalMetrics,
          searchLogs: logs
        }
      };

    } catch (error) {
      logger.error({ 
        requestId, 
        error: error.message, 
        stack: error.stack,
        duration: Date.now() - startTime 
      });
      throw error;
    }
  }

  /**
   * Determine if search should continue
   */
  shouldContinue(confidence, round) {
    return confidence < this.confidenceThreshold && round < this.maxRounds;
  }

  /**
   * Process raw search results (Phase 1 simple implementation)
   */
  processSearchResults(results) {
    return results.map(result => ({
      title: result.title,
      url: result.url,
      snippet: result.snippet,
      position: result.position
    }));
  }

  /**
   * Simple confidence estimation (Phase 1)
   * Phase 3 will replace this with expert-selector
   */
  estimateConfidence(results, query) {
    if (results.length === 0) return 0;
    
    // Simple heuristic: more results = higher confidence
    // Capped at 0.9 for Phase 1 (never reaches threshold on first round)
    const baseConfidence = Math.min(results.length / 10, 0.9);
    
    // Boost confidence if we have authoritative sources
    const hasAuthoritativeSources = results.some(r => 
      r.url && (
        r.url.includes('.edu') || 
        r.url.includes('.gov') ||
        r.url.includes('wikipedia.org')
      )
    );
    
    return hasAuthoritativeSources ? 
      Math.min(baseConfidence + 0.2, 1.0) : 
      baseConfidence;
  }

  /**
   * Generate unique request ID for tracking
   */
  generateRequestId() {
    return `req-${Date.now()}-${Math.random().toString(36).substr(2, 9)}`;
  }
}

// Export singleton instance
module.exports = new WorkflowEngine(require('../config'));