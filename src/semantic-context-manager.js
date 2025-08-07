const logger = require('./utils/logger');

/**
 * SemanticContextManager - Modern LLM-inspired context understanding
 * Uses embeddings and attention mechanisms instead of pattern matching
 */
class SemanticContextManager {
  constructor() {
    this.embeddingModel = null;
    this.queryClassifier = null;
    this.attentionScorer = null;
    this.memoryManager = null;
    
    // Configuration
    this.maxContextTokens = parseInt(process.env.MAX_CONTEXT_TOKENS || '2000');
    this.baseThreshold = parseFloat(process.env.SEMANTIC_THRESHOLD || '0.4');
    this.confidenceThreshold = parseFloat(process.env.CONFIDENCE_THRESHOLD || '0.7');
    
    // Pre-computed query type embeddings (will be loaded lazily)
    this.queryTypeEmbeddings = new Map();
    this.intentEmbeddings = new Map();
    
    // Cache for performance
    this.embeddingCache = new Map();
    this.maxCacheSize = 1000;
    
    // Initialize components
    this.initializeAsync();
    
    logger.info({
      message: 'SemanticContextManager initialized',
      maxContextTokens: this.maxContextTokens,
      baseThreshold: this.baseThreshold
    });
  }

  /**
   * Initialize embedding model and other components asynchronously
   */
  async initializeAsync() {
    try {
      // Use a fallback approach if transformers aren't available
      this.embeddingModel = await this.loadEmbeddingModel();
      this.queryClassifier = new QueryClassifier(this);
      this.attentionScorer = new AttentionScorer();
      this.memoryManager = new AdaptiveMemoryManager(this);
      
      // Pre-compute query type embeddings
      await this.precomputeQueryTypes();
      
      logger.info('Semantic components initialized successfully');
    } catch (error) {
      logger.warn('Failed to initialize transformers, using fallback semantic analysis', { error: error.message });
      this.initializeFallbackMode();
    }
  }

  /**
   * Load embedding model with fallback
   */
  async loadEmbeddingModel() {
    try {
      // Try to load transformers
      const { pipeline } = await import('@xenova/transformers');
      return await pipeline('feature-extraction', 'Xenova/all-MiniLM-L6-v2');
    } catch (error) {
      // Fallback to simple text analysis
      return new FallbackEmbeddingModel();
    }
  }

  /**
   * Initialize fallback mode without transformers
   */
  initializeFallbackMode() {
    this.embeddingModel = new FallbackEmbeddingModel();
    this.queryClassifier = new QueryClassifier(this);
    this.attentionScorer = new AttentionScorer();
    this.memoryManager = new AdaptiveMemoryManager(this);
  }

  /**
   * Pre-compute embeddings for common query types and intents
   */
  async precomputeQueryTypes() {
    const queryTypes = {
      'follow-up': 'continuing previous topic, referring to it, building on earlier discussion',
      'new-topic': 'starting new subject, different domain, unrelated to previous',
      'comparison': 'comparing, contrasting, differences between, versus, better than',
      'elaboration': 'more details, explanation, examples, expand on, tell me more',
      'clarification': 'what do you mean, unclear, confusing, explain again',
      'application': 'how to use, practical example, real world application',
      'definition': 'what is, define, meaning of, explain the concept'
    };

    const intents = {
      'reference': 'it, this, that, the previous, mentioned before',
      'continuation': 'also, additionally, furthermore, next, then',
      'contrast': 'but, however, although, different, instead',
      'causation': 'because, since, due to, results in, leads to'
    };

    for (const [type, description] of Object.entries(queryTypes)) {
      this.queryTypeEmbeddings.set(type, await this.getEmbedding(description));
    }

    for (const [intent, description] of Object.entries(intents)) {
      this.intentEmbeddings.set(intent, await this.getEmbedding(description));
    }
  }

  /**
   * Main context analysis method - replaces pattern-based approach
   */
  async analyzeContextNeed(query, history) {
    if (!history || !history.messages || history.messages.length === 0) {
      return {
        needsContext: false,
        contextType: 'new-conversation',
        confidence: 1.0,
        selectedMessages: [],
        reasoning: 'No conversation history available'
      };
    }

    try {
      // 1. Get query embedding
      const queryEmbedding = await this.getEmbedding(query);
      
      // 2. Classify query type and intent
      const queryType = await this.queryClassifier.classify(queryEmbedding, query);
      const queryIntent = await this.classifyIntent(queryEmbedding, query);
      
      // 3. Get embeddings for conversation history
      const messageEmbeddings = await this.getHistoryEmbeddings(history.messages);
      
      // 4. Calculate semantic similarities
      const similarities = messageEmbeddings.map(msgEmbed => 
        this.cosineSimilarity(queryEmbedding, msgEmbed.embedding)
      );
      
      // 5. Use attention mechanism for relevance scoring
      const attentionScores = this.attentionScorer.score(
        queryEmbedding,
        messageEmbeddings,
        queryIntent,
        queryType
      );
      
      // 6. Determine context need with adaptive thresholds
      const contextAnalysis = this.determineContextNeed(
        query,
        queryType,
        queryIntent,
        similarities,
        attentionScores,
        history
      );
      
      // 7. Select relevant messages
      const selectedMessages = this.selectRelevantMessages(
        history.messages,
        attentionScores,
        contextAnalysis.threshold
      );

      return {
        needsContext: contextAnalysis.needsContext,
        contextType: queryType,
        intent: queryIntent,
        confidence: contextAnalysis.confidence,
        selectedMessages,
        reasoning: contextAnalysis.reasoning,
        metadata: {
          similarities,
          attentionScores,
          adaptiveThreshold: contextAnalysis.threshold
        }
      };

    } catch (error) {
      logger.error('Error in semantic context analysis', { error: error.message });
      // Fallback to simple analysis
      return this.fallbackContextAnalysis(query, history);
    }
  }

  /**
   * Get embedding for text with caching
   */
  async getEmbedding(text) {
    // Check cache first
    const cacheKey = this.hashString(text);
    if (this.embeddingCache.has(cacheKey)) {
      return this.embeddingCache.get(cacheKey);
    }

    let embedding;
    if (this.embeddingModel && typeof this.embeddingModel === 'function') {
      // Transformers model
      const result = await this.embeddingModel(text);
      embedding = result.data || result[0] || result;
    } else {
      // Fallback model
      embedding = this.embeddingModel.encode(text);
    }

    // Cache the result
    if (this.embeddingCache.size >= this.maxCacheSize) {
      // Remove oldest entry
      const firstKey = this.embeddingCache.keys().next().value;
      this.embeddingCache.delete(firstKey);
    }
    this.embeddingCache.set(cacheKey, embedding);

    return embedding;
  }

  /**
   * Get embeddings for conversation history
   */
  async getHistoryEmbeddings(messages) {
    const embeddings = [];
    
    for (let i = 0; i < messages.length; i++) {
      const msg = messages[i];
      const embedding = await this.getEmbedding(msg.content);
      embeddings.push({
        embedding,
        message: msg,
        index: i
      });
    }
    
    return embeddings;
  }

  /**
   * Classify query intent using semantic similarity
   */
  async classifyIntent(queryEmbedding, query) {
    let bestIntent = 'general';
    let bestScore = 0;

    for (const [intent, intentEmbedding] of this.intentEmbeddings) {
      const similarity = this.cosineSimilarity(queryEmbedding, intentEmbedding);
      if (similarity > bestScore) {
        bestIntent = intent;
        bestScore = similarity;
      }
    }

    // Enhanced reference detection for implicit follow-ups
    const queryLower = query.toLowerCase();
    const referenceWords = ['it', 'this', 'that', 'they', 'them', 'the previous', 'above', 'these', 'those'];
    const implicitFollowUps = [
      /^what are the main/i,
      /^how does.*work/i, 
      /^what.*types/i,
      /^give.*example/i,
      /^which.*better/i,
      /^how.*compare/i
    ];
    
    const hasExplicitReference = referenceWords.some(word => queryLower.includes(word));
    const isImplicitFollowUp = implicitFollowUps.some(pattern => pattern.test(query));
    
    if (hasExplicitReference && bestIntent !== 'reference') {
      bestIntent = 'reference';
      bestScore = Math.max(bestScore, 0.8);
    } else if (isImplicitFollowUp && bestIntent === 'general') {
      bestIntent = 'continuation';
      bestScore = Math.max(bestScore, 0.6);
    }

    return { type: bestIntent, confidence: bestScore };
  }

  /**
   * Determine context need with adaptive thresholds
   */
  determineContextNeed(query, queryType, queryIntent, similarities, attentionScores, history) {
    // Calculate adaptive threshold based on conversation characteristics
    const conversationLength = history.messages.length;
    const maxSimilarity = Math.max(...similarities);
    const maxAttention = Math.max(...attentionScores);
    
    // Start with base threshold
    let threshold = this.baseThreshold;
    
    // Adjust based on query type (more aggressive for context detection)
    const queryTypeAdjustments = {
      'follow-up': -0.2,       // Lower threshold for follow-ups
      'comparison': -0.15,     // Lower for comparisons
      'elaboration': -0.15,    // Lower for elaborations
      'clarification': -0.25,  // Much lower for clarifications
      'new-topic': 0.05,       // Slightly higher for new topics
      'definition': 0.0        // Neutral for definitions
    };
    
    threshold += queryTypeAdjustments[queryType] || 0;
    
    // Adjust based on intent (more aggressive for references)
    if (queryIntent.type === 'reference') {
      threshold -= 0.25; // Much lower for reference queries
    }
    if (queryIntent.type === 'continuation') {
      threshold -= 0.15; // Lower for continuation queries
    }
    
    // Adjust based on conversation length
    if (conversationLength < 4) {
      threshold -= 0.1; // More aggressive in short conversations
    }
    
    // Ensure threshold stays in reasonable range
    threshold = Math.max(0.1, Math.min(0.8, threshold));
    
    // Determine if context is needed
    const needsContext = maxSimilarity > threshold || maxAttention > threshold;
    
    // Calculate confidence
    const confidence = Math.max(maxSimilarity, maxAttention);
    
    // Generate reasoning
    let reasoning = `Query type: ${queryType}, Intent: ${queryIntent.type}. `;
    reasoning += `Max similarity: ${maxSimilarity.toFixed(3)}, threshold: ${threshold.toFixed(3)}. `;
    reasoning += needsContext ? 'Context needed.' : 'No context needed.';

    return {
      needsContext,
      threshold,
      confidence,
      reasoning
    };
  }

  /**
   * Select relevant messages based on attention scores
   */
  selectRelevantMessages(messages, attentionScores, threshold) {
    const selected = [];
    let tokenCount = 0;
    const tokensPerChar = 0.25; // Rough estimate
    
    // Create array of message-score pairs
    const messageScores = messages.map((msg, index) => ({
      message: msg,
      score: attentionScores[index],
      index
    }));
    
    // Sort by score (highest first)
    messageScores.sort((a, b) => b.score - a.score);
    
    // Select messages above threshold, respecting token budget
    for (const item of messageScores) {
      if (item.score > threshold) {
        const msgTokens = Math.ceil(item.message.content.length * tokensPerChar);
        if (tokenCount + msgTokens <= this.maxContextTokens) {
          selected.push({
            ...item.message,
            relevanceScore: item.score,
            originalIndex: item.index
          });
          tokenCount += msgTokens;
        }
      }
    }
    
    // Sort selected messages by original order for coherent context
    selected.sort((a, b) => a.originalIndex - b.originalIndex);
    
    return selected;
  }

  /**
   * Calculate cosine similarity between two embeddings
   */
  cosineSimilarity(a, b) {
    if (!a || !b || a.length !== b.length) return 0;
    
    let dotProduct = 0;
    let normA = 0;
    let normB = 0;
    
    for (let i = 0; i < a.length; i++) {
      dotProduct += a[i] * b[i];
      normA += a[i] * a[i];
      normB += b[i] * b[i];
    }
    
    const magnitude = Math.sqrt(normA) * Math.sqrt(normB);
    return magnitude === 0 ? 0 : dotProduct / magnitude;
  }

  /**
   * Build contextual prompt using semantic understanding
   */
  buildContextualPrompt(query, contextAnalysis, basePrompt = '') {
    if (!contextAnalysis.needsContext || contextAnalysis.selectedMessages.length === 0) {
      return basePrompt + `\n\nCurrent question: ${query}`;
    }

    let prompt = '';
    
    // Add context summary
    const contextSummary = this.generateContextSummary(contextAnalysis);
    if (contextSummary) {
      prompt += `[Conversation Context: ${contextSummary}]\n\n`;
    }
    
    // Add continuity indicator based on semantic analysis
    const continuityHints = {
      'follow-up': 'This question builds upon previous discussion.',
      'comparison': 'This question asks for comparison with previous topics.',
      'elaboration': 'This question seeks more detail about previous topics.',
      'clarification': 'This question asks for clarification of previous responses.',
      'reference': 'This question refers to previously mentioned items.'
    };
    
    if (continuityHints[contextAnalysis.contextType]) {
      prompt += `[Note: ${continuityHints[contextAnalysis.contextType]}]\n\n`;
    }
    
    // Add relevant previous messages
    if (contextAnalysis.selectedMessages.length > 0) {
      prompt += 'Relevant conversation history:\n';
      contextAnalysis.selectedMessages.forEach((msg, index) => {
        const role = msg.role === 'user' ? 'User' : 'Assistant';
        const content = msg.content.length > 400 ? 
          msg.content.substring(0, 400) + '...' : msg.content;
        prompt += `${role}: ${content}\n`;
      });
      prompt += '\n';
    }
    
    // Add current query
    prompt += `Current question: ${query}\n\n`;
    
    // Add base prompt
    if (basePrompt) {
      prompt += basePrompt;
    }
    
    return prompt;
  }

  /**
   * Generate semantic context summary
   */
  generateContextSummary(contextAnalysis) {
    if (!contextAnalysis.selectedMessages || contextAnalysis.selectedMessages.length === 0) {
      return null;
    }
    
    const topics = new Set();
    const entities = new Set();
    
    contextAnalysis.selectedMessages.forEach(msg => {
      // Extract key terms from high-relevance messages
      if (msg.relevanceScore > 0.6) {
        const words = msg.content.split(/\s+/)
          .filter(word => word.length > 3)
          .filter(word => /^[A-Za-z]/.test(word));
        
        words.slice(0, 3).forEach(word => topics.add(word.toLowerCase()));
      }
    });
    
    const summary = [];
    if (topics.size > 0) {
      summary.push(`Topics: ${Array.from(topics).slice(0, 3).join(', ')}`);
    }
    
    return summary.join('. ');
  }

  /**
   * Fallback context analysis when semantic analysis fails
   */
  fallbackContextAnalysis(query, history) {
    const queryLower = query.toLowerCase();
    const referenceWords = ['it', 'this', 'that', 'they', 'them'];
    const hasReference = referenceWords.some(word => queryLower.includes(word));
    
    return {
      needsContext: hasReference,
      contextType: hasReference ? 'reference' : 'new-topic',
      confidence: hasReference ? 0.8 : 0.2,
      selectedMessages: hasReference ? history.messages.slice(-2) : [],
      reasoning: 'Fallback analysis: ' + (hasReference ? 'Reference detected' : 'No reference detected')
    };
  }

  /**
   * Simple hash function for caching
   */
  hashString(str) {
    let hash = 0;
    if (str.length === 0) return hash;
    for (let i = 0; i < str.length; i++) {
      const char = str.charCodeAt(i);
      hash = ((hash << 5) - hash) + char;
      hash = hash & hash; // Convert to 32-bit integer
    }
    return hash;
  }
}

/**
 * Query Classifier using semantic similarity
 */
class QueryClassifier {
  constructor(semanticManager) {
    this.semanticManager = semanticManager;
  }

  async classify(queryEmbedding, queryText) {
    let bestType = 'general';
    let bestScore = 0;

    // Wait for query type embeddings to be available
    if (this.semanticManager.queryTypeEmbeddings.size === 0) {
      await new Promise(resolve => setTimeout(resolve, 100)); // Wait briefly
    }

    for (const [type, typeEmbedding] of this.semanticManager.queryTypeEmbeddings) {
      const similarity = this.semanticManager.cosineSimilarity(queryEmbedding, typeEmbedding);
      if (similarity > bestScore) {
        bestType = type;
        bestScore = similarity;
      }
    }

    // Override with explicit patterns if very clear
    const queryLower = queryText.toLowerCase();
    if (queryLower.startsWith('what is') || queryLower.startsWith('define')) {
      return 'definition';
    }
    if (queryLower.includes('example') || queryLower.includes('show me')) {
      return 'application';
    }
    if (queryLower.startsWith('what are the') || queryLower.startsWith('how does') || 
        queryLower.startsWith('which') || queryLower.includes('give me an example')) {
      return 'follow-up';
    }
    if (queryLower.includes('compare') || queryLower.includes('difference') || queryLower.includes('versus')) {
      return 'comparison';
    }

    return bestType;
  }
}

/**
 * Attention Scorer for dynamic relevance assessment
 */
class AttentionScorer {
  score(queryEmbedding, messageEmbeddings, queryIntent, queryType) {
    return messageEmbeddings.map((msgData, index) => {
      const msgEmbedding = msgData.embedding;
      const message = msgData.message;
      
      // Semantic similarity (main signal)
      const semanticScore = this.cosineSimilarity(queryEmbedding, msgEmbedding);
      
      // Recency score with adaptive decay
      const position = messageEmbeddings.length - index - 1;
      const recencyScore = Math.exp(-0.1 * position);
      
      // Intent alignment bonus
      const intentBonus = this.calculateIntentBonus(message, queryIntent);
      
      // Query type bonus
      const typeBonus = this.calculateTypeBonus(message, queryType);
      
      // Weighted combination
      return (
        semanticScore * 0.6 +
        recencyScore * 0.2 +
        intentBonus * 0.1 +
        typeBonus * 0.1
      );
    });
  }

  cosineSimilarity(a, b) {
    if (!a || !b || a.length !== b.length) return 0;
    
    let dotProduct = 0;
    let normA = 0;
    let normB = 0;
    
    for (let i = 0; i < a.length; i++) {
      dotProduct += a[i] * b[i];
      normA += a[i] * a[i];
      normB += b[i] * b[i];
    }
    
    const magnitude = Math.sqrt(normA) * Math.sqrt(normB);
    return magnitude === 0 ? 0 : dotProduct / magnitude;
  }

  calculateIntentBonus(message, queryIntent) {
    if (!queryIntent || !queryIntent.type) return 0;
    
    const content = message.content.toLowerCase();
    const intentBonuses = {
      'reference': content.includes('this') || content.includes('that') ? 0.1 : 0,
      'continuation': content.includes('also') || content.includes('additionally') ? 0.1 : 0,
      'contrast': content.includes('however') || content.includes('but') ? 0.1 : 0
    };
    
    return intentBonuses[queryIntent.type] || 0;
  }

  calculateTypeBonus(message, queryType) {
    const content = message.content.toLowerCase();
    const typeBonuses = {
      'definition': content.includes('is defined') || content.includes('means') ? 0.1 : 0,
      'comparison': content.includes('compared') || content.includes('versus') ? 0.1 : 0,
      'application': content.includes('example') || content.includes('used') ? 0.1 : 0
    };
    
    return typeBonuses[queryType] || 0;
  }
}

/**
 * Adaptive Memory Manager with semantic chunking
 */
class AdaptiveMemoryManager {
  constructor(semanticManager) {
    this.semanticManager = semanticManager;
    this.chunks = [];
    this.chunkingThreshold = 0.6;
  }

  async addMessage(message) {
    const messageEmbedding = await this.semanticManager.getEmbedding(message.content);
    
    if (this.chunks.length === 0) {
      this.chunks.push(new ConversationChunk([message], messageEmbedding));
      return;
    }

    // Find most similar chunk
    let bestChunk = null;
    let bestSimilarity = 0;

    for (const chunk of this.chunks) {
      const similarity = this.semanticManager.cosineSimilarity(messageEmbedding, chunk.centroidEmbedding);
      if (similarity > bestSimilarity) {
        bestSimilarity = similarity;
        bestChunk = chunk;
      }
    }

    if (bestSimilarity > this.chunkingThreshold) {
      bestChunk.addMessage(message, messageEmbedding);
    } else {
      this.chunks.push(new ConversationChunk([message], messageEmbedding));
    }

    this.pruneChunks();
  }

  async getRelevantChunks(queryEmbedding, maxChunks = 3) {
    const chunkScores = this.chunks.map(chunk => ({
      chunk,
      score: this.semanticManager.cosineSimilarity(queryEmbedding, chunk.centroidEmbedding)
    }));

    return chunkScores
      .sort((a, b) => b.score - a.score)
      .slice(0, maxChunks)
      .map(item => item.chunk);
  }

  pruneChunks() {
    const maxChunks = 10;
    if (this.chunks.length > maxChunks) {
      this.chunks = this.chunks.slice(-maxChunks);
    }
  }
}

/**
 * Conversation Chunk for semantic memory management
 */
class ConversationChunk {
  constructor(messages, initialEmbedding) {
    this.messages = [...messages];
    this.centroidEmbedding = initialEmbedding;
    this.createdAt = Date.now();
  }

  addMessage(message, messageEmbedding) {
    this.messages.push(message);
    this.updateCentroid(messageEmbedding);
  }

  updateCentroid(newEmbedding) {
    if (!this.centroidEmbedding || !newEmbedding) return;
    
    // Update centroid as weighted average
    const weight = 1 / this.messages.length;
    for (let i = 0; i < this.centroidEmbedding.length; i++) {
      this.centroidEmbedding[i] = this.centroidEmbedding[i] * (1 - weight) + newEmbedding[i] * weight;
    }
  }
}

/**
 * Fallback Embedding Model when transformers aren't available
 */
class FallbackEmbeddingModel {
  constructor() {
    this.dimension = 128;
    this.vocabulary = new Map();
    this.idfScores = new Map();
  }

  encode(text) {
    // Simple TF-IDF-like embedding
    const words = text.toLowerCase().split(/\s+/).filter(w => w.length > 2);
    const embedding = new Array(this.dimension).fill(0);
    
    words.forEach((word, index) => {
      const hash = this.hashWord(word);
      const position = Math.abs(hash) % this.dimension;
      embedding[position] += 1 / (1 + index * 0.1); // Position-weighted
    });
    
    // Normalize
    const norm = Math.sqrt(embedding.reduce((sum, val) => sum + val * val, 0));
    return norm > 0 ? embedding.map(val => val / norm) : embedding;
  }

  hashWord(word) {
    let hash = 0;
    for (let i = 0; i < word.length; i++) {
      const char = word.charCodeAt(i);
      hash = ((hash << 5) - hash) + char;
      hash = hash & hash;
    }
    return hash;
  }
}

module.exports = SemanticContextManager;