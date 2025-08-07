const logger = require('./utils/logger');
const SemanticContextManager = require('./semantic-context-manager');
const ModernContextManager = require('./modern-context-manager');

/**
 * ContextManager - Intelligently selects and manages conversation context
 * Scores relevance, manages token limits, and creates contextual prompts
 */
class ContextManager {
  constructor(options = {}) {
    this.maxContextTokens = options.maxContextTokens || 
      parseInt(process.env.MAX_CONTEXT_TOKENS || '2000');
    this.relevanceThreshold = options.relevanceThreshold || 
      parseFloat(process.env.CONTEXT_RELEVANCE_THRESHOLD || '0.3');
    this.enableSummary = options.enableSummary !== false && 
      process.env.ENABLE_CONTEXT_SUMMARY !== 'false';
    
    // Choose context analysis mode
    this.useSemanticMode = process.env.USE_SEMANTIC_CONTEXT !== 'false';
    this.useModernMode = process.env.USE_MODERN_CONTEXT === 'true';
    
    // Initialize context managers based on configuration
    if (this.useModernMode) {
      this.modernManager = new ModernContextManager({ 
        maxContextTokens: this.maxContextTokens 
      });
    } else if (this.useSemanticMode) {
      this.semanticManager = new SemanticContextManager();
    }
    
    // Approximate tokens per character (rough estimate)
    this.tokensPerChar = 0.25;
    
    logger.info({
      message: 'ContextManager initialized',
      mode: this.useModernMode ? 'modern-2025' : (this.useSemanticMode ? 'semantic' : 'pattern-based'),
      maxContextTokens: this.maxContextTokens,
      relevanceThreshold: this.relevanceThreshold,
      enableSummary: this.enableSummary
    });
  }

  /**
   * Select relevant context from history for current query
   * @param {string} query - Current user query
   * @param {Object} history - Conversation history
   * @param {number} maxTokens - Maximum tokens for context
   * @returns {Object} Selected context with metadata
   */
  async selectRelevant(query, history, maxTokens = this.maxContextTokens) {
    if (!history || !history.messages || history.messages.length === 0) {
      return {
        messages: [],
        summary: null,
        hasContext: false,
        continuityType: 'new'
      };
    }

    // Use modern 2025 approach if enabled
    if (this.useModernMode && this.modernManager) {
      try {
        const analysis = await this.modernManager.analyzeContextNeed(query, history);
        return {
          messages: analysis.selectedMessages,
          summary: null,
          hasContext: analysis.needsContext,
          continuityType: analysis.metadata?.contextType || 'modern',
          confidence: analysis.confidence,
          reasoning: analysis.reasoning,
          topics: history.metadata?.topics || [],
          entities: history.metadata?.entities || [],
          expertsUsed: history.metadata?.expertsUsed || [],
          modernData: analysis.metadata
        };
      } catch (error) {
        logger.warn('Modern context analysis failed, falling back to semantic', { error: error.message });
      }
    }

    // Use semantic mode if available
    if (this.useSemanticMode && this.semanticManager) {
      try {
        return await this.selectRelevantSemantic(query, history, maxTokens);
      } catch (error) {
        logger.warn('Semantic analysis failed, falling back to pattern-based', { error: error.message });
        // Fall through to pattern-based approach
      }
    }

    // Pattern-based approach (original)
    return this.selectRelevantPatternBased(query, history, maxTokens);
  }

  /**
   * Semantic-based context selection
   */
  async selectRelevantSemantic(query, history, maxTokens) {
    const contextAnalysis = await this.semanticManager.analyzeContextNeed(query, history);
    
    return {
      messages: contextAnalysis.selectedMessages,
      summary: this.enableSummary ? this.generateSemanticSummary(contextAnalysis) : null,
      hasContext: contextAnalysis.needsContext,
      continuityType: contextAnalysis.contextType,
      confidence: contextAnalysis.confidence,
      reasoning: contextAnalysis.reasoning,
      topics: history.metadata?.topics || [],
      entities: history.metadata?.entities || [],
      expertsUsed: history.metadata?.expertsUsed || [],
      semanticData: {
        intent: contextAnalysis.intent,
        adaptiveThreshold: contextAnalysis.metadata?.adaptiveThreshold
      }
    };
  }

  /**
   * Pattern-based context selection (original implementation)
   */
  selectRelevantPatternBased(query, history, maxTokens) {
    // Score and rank messages
    const scoredMessages = this.scoreMessages(query, history.messages);
    
    // Detect conversation continuity type
    const continuityType = this.detectContinuityType(query, history);
    
    // Select messages within token budget
    const selectedMessages = this.selectWithinBudget(scoredMessages, maxTokens);
    
    // Generate summary if enabled and helpful
    const summary = this.enableSummary ? 
      this.generateContextSummary(history, selectedMessages) : null;
    
    return {
      messages: selectedMessages,
      summary,
      hasContext: selectedMessages.length > 0,
      continuityType,
      topics: history.metadata?.topics || [],
      entities: history.metadata?.entities || [],
      expertsUsed: history.metadata?.expertsUsed || []
    };
  }

  /**
   * Score messages based on relevance to current query
   * @param {string} query - Current query
   * @param {Array} messages - Previous messages
   * @returns {Array} Messages with relevance scores
   */
  scoreMessages(query, messages) {
    const queryLower = query.toLowerCase();
    const queryWords = new Set(queryLower.split(/\s+/));
    
    return messages.map((msg, index) => {
      let score = 0;
      const msgLower = msg.content.toLowerCase();
      
      // Recency score (exponential decay)
      const recencyScore = Math.exp(-0.2 * (messages.length - index - 1));
      score += recencyScore * 0.3;
      
      // Word overlap score
      const msgWords = new Set(msgLower.split(/\s+/));
      const overlap = [...queryWords].filter(w => msgWords.has(w)).length;
      const overlapScore = overlap / Math.max(queryWords.size, 1);
      score += overlapScore * 0.4;
      
      // Entity continuity score
      if (msg.metadata?.entities) {
        const entityMatch = this.checkEntityOverlap(query, msg.metadata.entities);
        score += entityMatch * 0.2;
      }
      
      // Topic relevance score
      if (msg.metadata?.topics) {
        const topicMatch = this.checkTopicRelevance(query, msg.metadata.topics);
        score += topicMatch * 0.1;
      }
      
      // Boost score for assistant messages with high confidence
      if (msg.role === 'assistant' && msg.metadata?.confidence > 0.8) {
        score *= 1.2;
      }
      
      return {
        ...msg,
        relevanceScore: score,
        includeInContext: score >= this.relevanceThreshold
      };
    }).sort((a, b) => b.relevanceScore - a.relevanceScore);
  }

  /**
   * Detect the type of conversation continuity
   * @param {string} query - Current query
   * @param {Object} history - Conversation history
   * @returns {string} Continuity type
   */
  detectContinuityType(query, history) {
    const queryLower = query.toLowerCase();
    
    // Check for explicit follow-up patterns
    const followUpPatterns = [
      /^(what|how|why|when|where|who)\s+(about|regarding|concerning)/i,
      /^(tell|explain|describe)\s+me\s+more/i,
      /^(can|could|would)\s+you\s+(explain|elaborate|clarify)/i,
      /^(also|additionally|furthermore|moreover)/i,
      /^(but|however|although|though)\s+/i
    ];
    
    if (followUpPatterns.some(pattern => pattern.test(query))) {
      return 'explicit-followup';
    }
    
    // Check for pronoun references
    const pronounPatterns = [
      /\b(it|its|it's)\b/i,
      /\b(this|that|these|those)\b/i,
      /\b(they|them|their)\b/i,
      /\bthe (previous|last|above|mentioned)\b/i
    ];
    
    if (pronounPatterns.some(pattern => pattern.test(queryLower))) {
      return 'pronoun-reference';
    }
    
    // Check for topic continuation
    if (history.metadata?.topics && history.metadata.topics.length > 0) {
      const topicContinuation = history.metadata.topics.some(topic => 
        queryLower.includes(topic.toLowerCase())
      );
      if (topicContinuation) {
        return 'topic-continuation';
      }
    }
    
    // Check for entity references
    if (history.metadata?.entities && history.metadata.entities.length > 0) {
      const entityReference = history.metadata.entities.some(entity => 
        queryLower.includes(entity.toLowerCase())
      );
      if (entityReference) {
        return 'entity-reference';
      }
    }
    
    return 'new';
  }

  /**
   * Select messages within token budget
   * @param {Array} scoredMessages - Messages with scores
   * @param {number} maxTokens - Token budget
   * @returns {Array} Selected messages
   */
  selectWithinBudget(scoredMessages, maxTokens) {
    const selected = [];
    let tokenCount = 0;
    
    // Always try to include the most recent exchange if relevant
    const recentExchange = scoredMessages.slice(0, 2)
      .filter(msg => msg.includeInContext);
    
    for (const msg of recentExchange) {
      const msgTokens = this.estimateTokens(msg.content);
      if (tokenCount + msgTokens <= maxTokens) {
        selected.push(msg);
        tokenCount += msgTokens;
      }
    }
    
    // Add other relevant messages by score
    const otherMessages = scoredMessages.slice(2)
      .filter(msg => msg.includeInContext && !selected.includes(msg));
    
    for (const msg of otherMessages) {
      const msgTokens = this.estimateTokens(msg.content);
      if (tokenCount + msgTokens <= maxTokens) {
        selected.push(msg);
        tokenCount += msgTokens;
      } else {
        // Try to include a truncated version
        const truncated = this.truncateMessage(msg, maxTokens - tokenCount);
        if (truncated) {
          selected.push(truncated);
          break;
        }
      }
    }
    
    // Sort by timestamp to maintain conversation flow
    return selected.sort((a, b) => a.timestamp - b.timestamp);
  }

  /**
   * Generate a summary of the conversation context
   * @param {Object} history - Full conversation history
   * @param {Array} selectedMessages - Messages selected for context
   * @returns {string|null} Context summary
   */
  generateContextSummary(history, selectedMessages) {
    if (!history.messages || history.messages.length < 3) {
      return null;
    }
    
    const topics = history.metadata?.topics || [];
    const entities = history.metadata?.entities || [];
    const experts = history.metadata?.expertsUsed || [];
    
    let summary = [];
    
    if (topics.length > 0) {
      summary.push(`Topics discussed: ${topics.slice(0, 3).join(', ')}`);
    }
    
    if (entities.length > 0) {
      summary.push(`Key entities: ${entities.slice(0, 5).join(', ')}`);
    }
    
    if (experts.length > 0) {
      summary.push(`Experts consulted: ${experts.slice(0, 3).join(', ')}`);
    }
    
    // Add key findings from high-confidence responses
    const keyFindings = history.messages
      .filter(msg => msg.role === 'assistant' && msg.metadata?.confidence > 0.85)
      .slice(-2)
      .map(msg => this.extractKeyPoint(msg.content));
    
    if (keyFindings.length > 0) {
      summary.push(`Key findings: ${keyFindings.join('; ')}`);
    }
    
    return summary.length > 0 ? summary.join('. ') : null;
  }

  /**
   * Build a contextual prompt with conversation history
   * @param {string} query - Current query
   * @param {Object} context - Selected context
   * @param {string} basePrompt - Base prompt to enhance
   * @returns {string} Enhanced prompt with context
   */
  buildContextualPrompt(query, context, basePrompt = '') {
    // Use modern prompt building if available
    if (this.useModernMode && this.modernManager && context.modernData) {
      return this.modernManager.buildContextualPrompt(query, context, basePrompt);
    }
    
    // Use semantic prompt building if available and semantic data present
    if (this.useSemanticMode && this.semanticManager && context.semanticData) {
      return this.semanticManager.buildContextualPrompt(query, context, basePrompt);
    }
    
    // Original pattern-based prompt building
    return this.buildPatternBasedPrompt(query, context, basePrompt);
  }

  /**
   * Generate semantic context summary
   */
  generateSemanticSummary(contextAnalysis) {
    if (!contextAnalysis || !contextAnalysis.selectedMessages) return null;
    
    const summary = [];
    
    if (contextAnalysis.intent) {
      summary.push(`Intent: ${contextAnalysis.intent.type}`);
    }
    
    if (contextAnalysis.contextType !== 'general') {
      summary.push(`Query type: ${contextAnalysis.contextType}`);
    }
    
    if (contextAnalysis.confidence) {
      summary.push(`Confidence: ${(contextAnalysis.confidence * 100).toFixed(0)}%`);
    }
    
    return summary.length > 0 ? summary.join(', ') : null;
  }

  /**
   * Original pattern-based prompt building
   */
  buildPatternBasedPrompt(query, context, basePrompt = '') {
    let prompt = '';
    
    // Add conversation summary if available
    if (context.summary) {
      prompt += `[Conversation Context: ${context.summary}]\n\n`;
    }
    
    // Add continuity indicator
    if (context.continuityType !== 'new') {
      const continuityHints = {
        'explicit-followup': 'This is a follow-up question.',
        'pronoun-reference': 'This question references previous discussion.',
        'topic-continuation': 'This continues the previous topic.',
        'entity-reference': 'This refers to previously mentioned entities.'
      };
      prompt += `[Note: ${continuityHints[context.continuityType]}]\n\n`;
    }
    
    // Add relevant previous messages
    if (context.messages && context.messages.length > 0) {
      prompt += 'Previous relevant exchanges:\n';
      context.messages.forEach(msg => {
        const role = msg.role === 'user' ? 'User' : 'Assistant';
        const content = msg.content.length > 500 ? 
          msg.content.substring(0, 500) + '...' : msg.content;
        prompt += `${role}: ${content}\n`;
      });
      prompt += '\n';
    }
    
    // Add current query
    prompt += `Current question: ${query}\n\n`;
    
    // Add base prompt if provided
    if (basePrompt) {
      prompt += basePrompt;
    }
    
    // Add context-aware instructions
    if (context.hasContext) {
      prompt += '\nConsider the conversation context and maintain continuity with previous responses. ';
      if (context.expertsUsed && context.expertsUsed.length > 0) {
        prompt += `Continue as ${context.expertsUsed[context.expertsUsed.length - 1]} if appropriate for this question.`;
      }
    }
    
    return prompt;
  }

  /**
   * Estimate token count for text
   * @param {string} text - Text to estimate
   * @returns {number} Estimated token count
   */
  estimateTokens(text) {
    return Math.ceil(text.length * this.tokensPerChar);
  }

  /**
   * Truncate a message to fit within token budget
   * @param {Object} message - Message to truncate
   * @param {number} maxTokens - Maximum tokens allowed
   * @returns {Object|null} Truncated message or null
   */
  truncateMessage(message, maxTokens) {
    const targetLength = Math.floor(maxTokens / this.tokensPerChar);
    if (targetLength < 50) return null; // Too short to be useful
    
    return {
      ...message,
      content: message.content.substring(0, targetLength) + '... [truncated]',
      truncated: true
    };
  }

  /**
   * Check entity overlap between query and message entities
   * @param {string} query - Current query
   * @param {Array} entities - Entities from message
   * @returns {number} Overlap score (0-1)
   */
  checkEntityOverlap(query, entities) {
    if (!entities || entities.length === 0) return 0;
    
    const queryLower = query.toLowerCase();
    const matches = entities.filter(entity => 
      queryLower.includes(entity.toLowerCase())
    );
    
    return matches.length / entities.length;
  }

  /**
   * Check topic relevance between query and message topics
   * @param {string} query - Current query
   * @param {Array} topics - Topics from message
   * @returns {number} Relevance score (0-1)
   */
  checkTopicRelevance(query, topics) {
    if (!topics || topics.length === 0) return 0;
    
    const queryLower = query.toLowerCase();
    const relevantTopics = topics.filter(topic => 
      queryLower.includes(topic.toLowerCase()) ||
      this.areTopicsRelated(query, topic)
    );
    
    return relevantTopics.length / topics.length;
  }

  /**
   * Check if topics are related (simple heuristic)
   * @param {string} query - Query text
   * @param {string} topic - Topic to check
   * @returns {boolean} Whether topics are related
   */
  areTopicsRelated(query, topic) {
    // Simple related terms mapping (could be enhanced)
    const relatedTerms = {
      'programming': ['code', 'software', 'development', 'algorithm'],
      'mathematics': ['math', 'calculation', 'equation', 'formula'],
      'science': ['research', 'experiment', 'theory', 'hypothesis'],
      'ai': ['machine learning', 'neural', 'model', 'artificial intelligence']
    };
    
    const queryLower = query.toLowerCase();
    const topicLower = topic.toLowerCase();
    
    for (const [key, terms] of Object.entries(relatedTerms)) {
      if (topicLower.includes(key) && terms.some(term => queryLower.includes(term))) {
        return true;
      }
    }
    
    return false;
  }

  /**
   * Extract key point from a message (first significant sentence)
   * @param {string} content - Message content
   * @returns {string} Key point
   */
  extractKeyPoint(content) {
    const sentences = content.match(/[^.!?]+[.!?]+/g) || [];
    const keypoint = sentences.find(s => s.length > 30 && s.length < 200) || sentences[0];
    return keypoint ? keypoint.trim().substring(0, 150) : '';
  }
}

// Export singleton instance
module.exports = new ContextManager();