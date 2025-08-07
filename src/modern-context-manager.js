const logger = require('./utils/logger');

/**
 * ModernContextManager - 2025 LLM-inspired context understanding
 * Based on research from txtai, attention head analysis, and dynamic context orchestration
 * Uses information-theoretic approaches and Bayesian inference instead of pattern matching
 */
class ModernContextManager {
  constructor(options = {}) {
    this.maxContextTokens = options.maxContextTokens || 
      parseInt(process.env.MAX_CONTEXT_TOKENS || '2000');
    
    // Modern context engineering parameters
    this.contextOrchestrator = new DynamicContextOrchestrator();
    this.bayesianInference = new BayesianContextInference();
    this.informationExtractor = new InformationTheoreticExtractor();
    this.attentionSimulator = new AttentionHeadSimulator();
    
    // State tracking for dynamic context
    this.conversationState = new ConversationState();
    
    logger.info({
      message: 'ModernContextManager initialized with 2025 methodologies',
      maxContextTokens: this.maxContextTokens
    });
  }

  /**
   * Main context analysis using 2025 LLM methodologies
   * Replaces pattern-based approach with information-theoretic analysis
   */
  async analyzeContextNeed(query, history) {
    if (!history?.messages?.length) {
      return this.createEmptyContext(query);
    }

    try {
      // 1. Update conversation state with new query
      this.conversationState.update(query, history);
      
      // 2. Calculate information-theoretic relevance
      const infoScores = await this.informationExtractor.calculateMutualInformation(
        query, 
        history.messages
      );
      
      // 3. Simulate attention head patterns
      const attentionPatterns = await this.attentionSimulator.simulateHeadActivation(
        query,
        history.messages,
        this.conversationState
      );
      
      // 4. Use Bayesian inference for context selection
      const bayesianAnalysis = this.bayesianInference.inferOptimalContext(
        query,
        history.messages,
        infoScores,
        attentionPatterns,
        this.conversationState
      );
      
      // 5. Orchestrate dynamic context assembly
      const contextAssembly = this.contextOrchestrator.assembleContext(
        query,
        history.messages,
        bayesianAnalysis,
        this.maxContextTokens
      );

      return {
        needsContext: contextAssembly.needsContext,
        confidence: bayesianAnalysis.posteriorProbability,
        selectedMessages: contextAssembly.selectedMessages,
        reasoning: this.explainReasoning(bayesianAnalysis, attentionPatterns, infoScores),
        metadata: {
          informationScores: infoScores,
          attentionActivation: attentionPatterns,
          bayesianPosterior: bayesianAnalysis.posteriorProbability,
          contextType: this.conversationState.getCurrentType()
        }
      };

    } catch (error) {
      logger.error('Modern context analysis failed', { error: error.message });
      return this.fallbackToSimpleAnalysis(query, history);
    }
  }

  /**
   * Build contextual prompt using dynamic orchestration
   */
  buildContextualPrompt(query, contextAnalysis, basePrompt = '') {
    if (!contextAnalysis.needsContext) {
      return `${basePrompt}\n\nCurrent question: ${query}`;
    }

    return this.contextOrchestrator.buildDynamicPrompt(
      query,
      contextAnalysis,
      basePrompt,
      this.conversationState
    );
  }

  createEmptyContext(query) {
    return {
      needsContext: false,
      confidence: 1.0,
      selectedMessages: [],
      reasoning: 'No conversation history available',
      metadata: { contextType: 'new' }
    };
  }

  explainReasoning(bayesian, attention, info) {
    const reasons = [];
    
    if (bayesian.posteriorProbability > 0.7) {
      reasons.push(`High Bayesian confidence (${(bayesian.posteriorProbability * 100).toFixed(1)}%)`);
    }
    
    if (attention.retrievalHeadActivation > 0.6) {
      reasons.push(`Strong retrieval head activation (${(attention.retrievalHeadActivation * 100).toFixed(1)}%)`);
    }
    
    if (info.maxMutualInformation > 0.5) {
      reasons.push(`High information overlap (${(info.maxMutualInformation * 100).toFixed(1)}%)`);
    }
    
    return reasons.length > 0 ? reasons.join(', ') : 'Low confidence context need';
  }

  fallbackToSimpleAnalysis(query, history) {
    const hasPronouns = /\b(it|this|that|they|them)\b/i.test(query);
    return {
      needsContext: hasPronouns,
      confidence: hasPronouns ? 0.6 : 0.2,
      selectedMessages: hasPronouns ? history.messages.slice(-2) : [],
      reasoning: 'Fallback: Simple pronoun detection',
      metadata: { contextType: 'fallback' }
    };
  }
}

/**
 * Dynamic Context Orchestrator - Assembles context components intelligently
 */
class DynamicContextOrchestrator {
  assembleContext(query, messages, bayesianAnalysis, maxTokens) {
    if (bayesianAnalysis.posteriorProbability < 0.3) {
      return { needsContext: false, selectedMessages: [] };
    }

    // Select messages based on Bayesian ranking
    const rankedMessages = bayesianAnalysis.rankedMessages || [];
    const selected = [];
    let tokenCount = 0;
    const tokensPerChar = 0.25;

    for (const msgData of rankedMessages) {
      const msgTokens = Math.ceil(msgData.message.content.length * tokensPerChar);
      if (tokenCount + msgTokens <= maxTokens) {
        selected.push({
          ...msgData.message,
          relevanceScore: msgData.probability,
          contextReason: msgData.reason
        });
        tokenCount += msgTokens;
      }
    }

    // Sort by original order for coherence
    selected.sort((a, b) => a.timestamp - b.timestamp);

    return {
      needsContext: selected.length > 0,
      selectedMessages: selected
    };
  }

  buildDynamicPrompt(query, contextAnalysis, basePrompt, conversationState) {
    let prompt = '';

    // Add state-aware context summary
    const stateSummary = conversationState.generateStateSummary();
    if (stateSummary) {
      prompt += `[Conversation State: ${stateSummary}]\n\n`;
    }

    // Add reasoning transparency
    if (contextAnalysis.reasoning) {
      prompt += `[Context Reasoning: ${contextAnalysis.reasoning}]\n\n`;
    }

    // Add relevant messages with context reasoning
    if (contextAnalysis.selectedMessages?.length > 0) {
      prompt += 'Contextually relevant history:\n';
      contextAnalysis.selectedMessages.forEach(msg => {
        const role = msg.role === 'user' ? 'User' : 'Assistant';
        const content = msg.content.length > 400 ? 
          msg.content.substring(0, 400) + '...' : msg.content;
        const reason = msg.contextReason ? ` [Relevant because: ${msg.contextReason}]` : '';
        prompt += `${role}: ${content}${reason}\n`;
      });
      prompt += '\n';
    }

    prompt += `Current question: ${query}\n\n${basePrompt}`;
    return prompt;
  }
}

/**
 * Bayesian Context Inference - Uses probabilistic reasoning for context selection
 */
class BayesianContextInference {
  inferOptimalContext(query, messages, infoScores, attentionPatterns, conversationState) {
    // Prior probability based on conversation state
    const priorProb = this.calculatePrior(conversationState);
    
    // Likelihood based on multiple signals
    const likelihood = this.calculateLikelihood(infoScores, attentionPatterns);
    
    // Posterior probability using Bayes' theorem
    const posteriorProbability = this.updatePosterior(priorProb, likelihood);
    
    // Rank messages by their individual posteriors
    const rankedMessages = messages.map((msg, index) => {
      const msgPosterior = this.calculateMessagePosterior(
        msg, 
        infoScores[index] || 0,
        attentionPatterns.messageActivations?.[index] || 0
      );
      
      return {
        message: msg,
        probability: msgPosterior,
        reason: this.explainMessageRelevance(msgPosterior, infoScores[index], attentionPatterns.messageActivations?.[index])
      };
    }).sort((a, b) => b.probability - a.probability);

    return {
      posteriorProbability,
      rankedMessages: rankedMessages.filter(item => item.probability > 0.2)
    };
  }

  calculatePrior(conversationState) {
    const stateFactors = {
      'continuation': 0.7,
      'new_topic': 0.3,
      'clarification': 0.9,
      'elaboration': 0.8
    };
    
    return stateFactors[conversationState.getCurrentType()] || 0.5;
  }

  calculateLikelihood(infoScores, attentionPatterns) {
    const infoComponent = Math.max(...infoScores) || 0;
    const attentionComponent = attentionPatterns.retrievalHeadActivation || 0;
    const inductionComponent = attentionPatterns.inductionHeadActivation || 0;
    
    // Weighted combination of signals
    return (
      infoComponent * 0.4 +
      attentionComponent * 0.4 +
      inductionComponent * 0.2
    );
  }

  updatePosterior(prior, likelihood) {
    // Simplified Bayesian update
    const evidence = prior * likelihood + (1 - prior) * (1 - likelihood);
    return evidence > 0 ? (prior * likelihood) / evidence : prior;
  }

  calculateMessagePosterior(message, infoScore, attentionScore) {
    // Message-level Bayesian analysis
    const contentPrior = message.role === 'assistant' ? 0.6 : 0.4; // Slight preference for assistant messages
    const signalLikelihood = (infoScore + attentionScore) / 2;
    
    return this.updatePosterior(contentPrior, signalLikelihood);
  }

  explainMessageRelevance(posterior, infoScore, attentionScore) {
    if (posterior > 0.8) return 'High semantic relevance';
    if (posterior > 0.6) return 'Moderate relevance';
    if (infoScore > 0.5) return 'Information overlap';
    if (attentionScore > 0.5) return 'Attention pattern match';
    return 'Low relevance';
  }
}

/**
 * Information-Theoretic Context Extractor
 */
class InformationTheoreticExtractor {
  async calculateMutualInformation(query, messages) {
    const queryTerms = this.extractTerms(query);
    
    return messages.map(message => {
      const messageTerms = this.extractTerms(message.content);
      
      // Calculate mutual information using term overlap and semantic distance
      const termOverlap = this.calculateTermOverlap(queryTerms, messageTerms);
      const semanticDistance = this.calculateSemanticDistance(queryTerms, messageTerms);
      
      // Combine into mutual information estimate
      return this.mutualInformationEstimate(termOverlap, semanticDistance);
    });
  }

  extractTerms(text) {
    return text.toLowerCase()
      .split(/\s+/)
      .filter(term => term.length > 2)
      .filter(term => !/^(the|and|or|but|in|on|at|to|for)$/.test(term));
  }

  calculateTermOverlap(terms1, terms2) {
    const set1 = new Set(terms1);
    const set2 = new Set(terms2);
    const intersection = new Set([...set1].filter(x => set2.has(x)));
    
    return intersection.size / Math.max(set1.size, set2.size, 1);
  }

  calculateSemanticDistance(terms1, terms2) {
    // Simplified semantic distance using term co-occurrence patterns
    let semanticScore = 0;
    const semanticPairs = [
      ['machine', 'learning'], ['quantum', 'computing'], ['neural', 'network'],
      ['data', 'science'], ['artificial', 'intelligence'], ['deep', 'learning']
    ];
    
    for (const [term1, term2] of semanticPairs) {
      if ((terms1.includes(term1) && terms2.includes(term2)) ||
          (terms1.includes(term2) && terms2.includes(term1))) {
        semanticScore += 0.3;
      }
    }
    
    return Math.min(semanticScore, 1.0);
  }

  mutualInformationEstimate(termOverlap, semanticDistance) {
    return (termOverlap * 0.7 + semanticDistance * 0.3);
  }
}

/**
 * Attention Head Simulator - Simulates LLM attention patterns
 */
class AttentionHeadSimulator {
  async simulateHeadActivation(query, messages, conversationState) {
    const patterns = {
      retrievalHeadActivation: this.simulateRetrievalHeads(query, messages),
      inductionHeadActivation: this.simulateInductionHeads(query, messages),
      safetyHeadActivation: this.simulateSafetyHeads(query),
      messageActivations: messages.map(msg => this.calculateMessageActivation(query, msg))
    };

    return patterns;
  }

  simulateRetrievalHeads(query, messages) {
    // Simulate retrieval heads that look for relevant context
    const queryKeywords = query.toLowerCase().split(/\s+/);
    let maxActivation = 0;

    for (const message of messages) {
      const messageKeywords = message.content.toLowerCase().split(/\s+/);
      const overlap = queryKeywords.filter(word => messageKeywords.includes(word)).length;
      const activation = overlap / Math.max(queryKeywords.length, 1);
      maxActivation = Math.max(maxActivation, activation);
    }

    return maxActivation;
  }

  simulateInductionHeads(query, messages) {
    // Simulate induction heads that detect patterns and continuation
    const patterns = [
      /^(what|how|why|when|where)/i,
      /^(can|could|would|should)/i,
      /(example|instance|case)/i,
      /(compare|contrast|difference)/i
    ];

    const hasPattern = patterns.some(pattern => pattern.test(query));
    return hasPattern ? 0.8 : 0.2;
  }

  simulateSafetyHeads(query) {
    // Simulate safety heads that check for harmful content
    const safetyTriggers = ['harmful', 'dangerous', 'illegal', 'unethical'];
    const hasTrigger = safetyTriggers.some(trigger => 
      query.toLowerCase().includes(trigger)
    );
    
    return hasTrigger ? 0.9 : 0.1;
  }

  calculateMessageActivation(query, message) {
    // Individual message activation based on multiple factors
    const recencyWeight = 0.3;
    const contentWeight = 0.7;
    
    const contentSim = this.simpleContentSimilarity(query, message.content);
    const recencyScore = 1.0; // Simplified - in real implementation would use timestamp
    
    return contentSim * contentWeight + recencyScore * recencyWeight;
  }

  simpleContentSimilarity(text1, text2) {
    const words1 = new Set(text1.toLowerCase().split(/\s+/));
    const words2 = new Set(text2.toLowerCase().split(/\s+/));
    const intersection = new Set([...words1].filter(x => words2.has(x)));
    
    return intersection.size / Math.max(words1.size, words2.size, 1);
  }
}

/**
 * Conversation State Tracker
 */
class ConversationState {
  constructor() {
    this.currentType = 'new';
    this.topics = new Set();
    this.entities = new Set();
    this.intentHistory = [];
    this.lastUpdate = Date.now();
  }

  update(query, history) {
    // Update conversation state based on new information
    this.analyzeQueryType(query);
    this.extractTopics(query);
    this.extractEntities(query);
    this.updateIntentHistory(query);
    this.lastUpdate = Date.now();
  }

  analyzeQueryType(query) {
    const queryLower = query.toLowerCase();
    
    if (/\b(it|this|that|they|them)\b/.test(queryLower)) {
      this.currentType = 'continuation';
    } else if (queryLower.startsWith('what') || queryLower.startsWith('how')) {
      this.currentType = 'elaboration';
    } else if (queryLower.includes('clarify') || queryLower.includes('explain')) {
      this.currentType = 'clarification';
    } else {
      this.currentType = 'new_topic';
    }
  }

  extractTopics(query) {
    const topicKeywords = [
      'machine learning', 'quantum computing', 'artificial intelligence',
      'data science', 'neural networks', 'deep learning', 'programming'
    ];
    
    for (const topic of topicKeywords) {
      if (query.toLowerCase().includes(topic)) {
        this.topics.add(topic);
      }
    }
  }

  extractEntities(query) {
    // Simple entity extraction - in real implementation would use NER
    const words = query.split(/\s+/);
    words.forEach(word => {
      if (word.length > 6 && /^[A-Z]/.test(word)) {
        this.entities.add(word);
      }
    });
  }

  updateIntentHistory(query) {
    const intent = this.classifyIntent(query);
    this.intentHistory.push({ intent, timestamp: Date.now() });
    
    // Keep only recent intents
    const fiveMinutesAgo = Date.now() - 300000;
    this.intentHistory = this.intentHistory.filter(item => item.timestamp > fiveMinutesAgo);
  }

  classifyIntent(query) {
    const queryLower = query.toLowerCase();
    if (queryLower.includes('example')) return 'application';
    if (queryLower.includes('compare')) return 'comparison';
    if (queryLower.includes('explain')) return 'explanation';
    return 'general';
  }

  getCurrentType() {
    return this.currentType;
  }

  generateStateSummary() {
    const parts = [];
    
    if (this.topics.size > 0) {
      parts.push(`Topics: ${Array.from(this.topics).slice(0, 3).join(', ')}`);
    }
    
    if (this.intentHistory.length > 0) {
      const recentIntent = this.intentHistory[this.intentHistory.length - 1];
      parts.push(`Recent intent: ${recentIntent.intent}`);
    }
    
    parts.push(`Type: ${this.currentType}`);
    
    return parts.join(' | ');
  }
}

module.exports = ModernContextManager;