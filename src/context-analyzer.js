const logger = require('./utils/logger');

/**
 * ContextAnalyzer - Analyzes queries and responses for context extraction
 * Detects follow-ups, extracts entities/topics, identifies expert personas
 */
class ContextAnalyzer {
  constructor() {
    // Common follow-up patterns
    this.followUpPatterns = [
      { pattern: /^(what|how|why|when|where|who)\s+(about|regarding)/i, type: 'question-followup' },
      { pattern: /^(tell|explain|show)\s+me\s+more/i, type: 'elaboration-request' },
      { pattern: /^(can|could|would)\s+you\s+(explain|elaborate|clarify)/i, type: 'clarification' },
      { pattern: /^(also|additionally|furthermore|moreover)/i, type: 'addition' },
      { pattern: /^(but|however|although|though)\s+/i, type: 'contrast' },
      { pattern: /\b(it|its|it's)\s+(performance|capabilities|features|benefits)/i, type: 'pronoun-attribute' },
      { pattern: /\bcompare(d?)\s+(it|this|that|them)\s+to/i, type: 'comparison' },
      { pattern: /\b(the|that|this)\s+(previous|last|above|mentioned)/i, type: 'explicit-reference' },
      { pattern: /^(yes|no|correct|wrong|right),?\s+but/i, type: 'correction' },
      { pattern: /^wait,?\s+/i, type: 'interruption' }
    ];

    // Pronoun patterns that indicate context dependency
    this.pronounPatterns = [
      { pattern: /\b(it|its|it's)\b(?!\s+is\s+\w+ing)/i, referenceType: 'singular-entity' },
      { pattern: /\b(they|them|their|theirs)\b/i, referenceType: 'plural-entity' },
      { pattern: /\b(this|that)\b(?!\s+is\s+\w+ing)/i, referenceType: 'demonstrative' },
      { pattern: /\b(these|those)\b/i, referenceType: 'plural-demonstrative' },
      { pattern: /\bthe\s+(system|algorithm|method|approach|solution|problem)\b/i, referenceType: 'definite-concept' }
    ];

    // Expert persona patterns
    this.expertPatterns = [
      { domain: 'mathematics', keywords: ['equation', 'formula', 'calculate', 'theorem', 'proof', 'derivative', 'integral'] },
      { domain: 'physics', keywords: ['quantum', 'relativity', 'particle', 'force', 'energy', 'momentum', 'wave'] },
      { domain: 'chemistry', keywords: ['molecule', 'reaction', 'compound', 'element', 'bond', 'oxidation', 'catalyst'] },
      { domain: 'biology', keywords: ['cell', 'dna', 'protein', 'evolution', 'species', 'organism', 'gene'] },
      { domain: 'computer-science', keywords: ['algorithm', 'code', 'program', 'software', 'database', 'api', 'framework'] },
      { domain: 'medicine', keywords: ['diagnosis', 'treatment', 'symptom', 'disease', 'patient', 'therapy', 'medication'] },
      { domain: 'history', keywords: ['century', 'war', 'civilization', 'empire', 'revolution', 'historical', 'period'] },
      { domain: 'economics', keywords: ['market', 'economy', 'inflation', 'gdp', 'trade', 'finance', 'investment'] },
      { domain: 'philosophy', keywords: ['ethics', 'morality', 'consciousness', 'existence', 'knowledge', 'truth', 'logic'] },
      { domain: 'engineering', keywords: ['design', 'build', 'structure', 'system', 'optimize', 'efficiency', 'mechanism'] }
    ];
  }

  /**
   * Analyze a query for context dependencies
   * @param {string} query - User query
   * @param {Object} previousContext - Previous conversation context
   * @returns {Object} Analysis results
   */
  analyzeQuery(query, previousContext = null) {
    const analysis = {
      isFollowUp: false,
      followUpType: null,
      pronounReferences: [],
      requiredContext: [],
      suggestedDomain: null,
      entities: [],
      topics: [],
      contextDependency: 'low' // low, medium, high
    };

    // Check for follow-up patterns
    const followUp = this.detectFollowUp(query);
    if (followUp) {
      analysis.isFollowUp = true;
      analysis.followUpType = followUp.type;
      analysis.contextDependency = 'high';
    }

    // Check for pronoun references
    const pronouns = this.detectPronouns(query);
    if (pronouns.length > 0) {
      analysis.pronounReferences = pronouns;
      analysis.contextDependency = analysis.contextDependency === 'low' ? 'medium' : 'high';
      
      // Determine what context is needed
      pronouns.forEach(p => {
        if (p.referenceType === 'singular-entity' || p.referenceType === 'demonstrative') {
          analysis.requiredContext.push('last-mentioned-entity');
        } else if (p.referenceType === 'plural-entity' || p.referenceType === 'plural-demonstrative') {
          analysis.requiredContext.push('recent-entities-list');
        } else if (p.referenceType === 'definite-concept') {
          analysis.requiredContext.push('discussed-concept');
        }
      });
    }

    // Extract entities and topics
    analysis.entities = this.extractEntities(query);
    analysis.topics = this.extractTopics(query);

    // Determine suggested domain
    analysis.suggestedDomain = this.determineDomain(query);

    // Check if query is incomplete without context
    if (this.isIncompleteWithoutContext(query)) {
      analysis.contextDependency = 'high';
      analysis.requiredContext.push('previous-question');
    }

    return analysis;
  }

  /**
   * Analyze a response for metadata extraction
   * @param {string} response - Assistant response
   * @param {string} query - Original query
   * @returns {Object} Extracted metadata
   */
  analyzeResponse(response, query) {
    const metadata = {
      entities: [],
      topics: [],
      keyFacts: [],
      expertDomain: null,
      expertPersona: null,
      confidence: 0.5
    };

    // Extract entities from response
    metadata.entities = this.extractEntities(response);
    
    // Extract topics
    metadata.topics = this.extractTopics(response);
    
    // Extract key facts (sentences with specific patterns)
    metadata.keyFacts = this.extractKeyFacts(response);
    
    // Determine expert domain based on content
    metadata.expertDomain = this.determineDomain(response);
    
    // Extract expert persona if mentioned
    metadata.expertPersona = this.extractExpertPersona(response);
    
    // Estimate confidence based on response characteristics
    metadata.confidence = this.estimateConfidence(response);

    return metadata;
  }

  /**
   * Detect follow-up patterns in query
   * @param {string} query - User query
   * @returns {Object|null} Follow-up detection result
   */
  detectFollowUp(query) {
    for (const followUp of this.followUpPatterns) {
      if (followUp.pattern.test(query)) {
        return {
          type: followUp.type,
          matched: query.match(followUp.pattern)[0]
        };
      }
    }
    return null;
  }

  /**
   * Detect pronoun references in query
   * @param {string} query - User query
   * @returns {Array} Detected pronouns with types
   */
  detectPronouns(query) {
    const detected = [];
    
    for (const pronounDef of this.pronounPatterns) {
      const matches = query.match(new RegExp(pronounDef.pattern, 'gi'));
      if (matches) {
        matches.forEach(match => {
          detected.push({
            pronoun: match,
            referenceType: pronounDef.referenceType,
            position: query.indexOf(match)
          });
        });
      }
    }
    
    return detected.sort((a, b) => a.position - b.position);
  }

  /**
   * Extract entities from text (simple NER)
   * @param {string} text - Text to analyze
   * @returns {Array} Extracted entities
   */
  extractEntities(text) {
    const entities = [];
    
    // Extract capitalized words (proper nouns)
    const properNouns = text.match(/\b[A-Z][a-z]+(?:\s+[A-Z][a-z]+)*/g) || [];
    entities.push(...properNouns.filter(n => n.length > 2));
    
    // Extract acronyms
    const acronyms = text.match(/\b[A-Z]{2,}\b/g) || [];
    entities.push(...acronyms);
    
    // Extract quoted terms
    const quoted = text.match(/"([^"]+)"/g) || [];
    entities.push(...quoted.map(q => q.replace(/"/g, '')));
    
    // Extract technical terms (simple heuristic)
    const technicalPatterns = [
      /\b\w+(?:ation|ology|ometry|ography)\b/gi,
      /\b\w+(?:\.js|\.py|\.java|\.cpp)\b/gi,
      /\b(?:API|SDK|GUI|CLI|CPU|GPU|RAM)\b/g
    ];
    
    technicalPatterns.forEach(pattern => {
      const matches = text.match(pattern) || [];
      entities.push(...matches);
    });
    
    // Remove duplicates and filter
    return [...new Set(entities)]
      .filter(e => e.length > 2 && e.length < 50)
      .slice(0, 10); // Limit to top 10
  }

  /**
   * Extract topics from text
   * @param {string} text - Text to analyze
   * @returns {Array} Extracted topics
   */
  extractTopics(text) {
    const topics = [];
    const textLower = text.toLowerCase();
    
    // Check against known domains
    for (const expert of this.expertPatterns) {
      const matchCount = expert.keywords.filter(keyword => 
        textLower.includes(keyword)
      ).length;
      
      if (matchCount >= 2) {
        topics.push(expert.domain);
      }
    }
    
    // Extract noun phrases (simplified)
    const nounPhrases = text.match(/\b(?:the\s+)?(?:\w+\s+){0,2}\w+(?:ing|tion|ment|ness|ity)\b/gi) || [];
    topics.push(...nounPhrases.slice(0, 3).map(np => np.toLowerCase().trim()));
    
    return [...new Set(topics)].slice(0, 5);
  }

  /**
   * Extract key facts from response
   * @param {string} response - Assistant response
   * @returns {Array} Key facts
   */
  extractKeyFacts(response) {
    const facts = [];
    const sentences = response.match(/[^.!?]+[.!?]+/g) || [];
    
    const factPatterns = [
      /\b(?:is|are|was|were)\s+(?:defined|known|called|considered)\s+/i,
      /\b(?:means|indicates|suggests|shows|demonstrates)\s+that\s+/i,
      /\b(?:approximately|roughly|about|around)\s+[\d.]+/i,
      /\b(?:discovered|invented|created|developed)\s+(?:by|in)\s+/i,
      /\b(?:consists?\s+of|comprises?|includes?)\s+/i
    ];
    
    sentences.forEach(sentence => {
      if (sentence.length > 30 && sentence.length < 200) {
        for (const pattern of factPatterns) {
          if (pattern.test(sentence)) {
            facts.push(sentence.trim());
            break;
          }
        }
      }
    });
    
    return facts.slice(0, 5);
  }

  /**
   * Determine domain from text content
   * @param {string} text - Text to analyze
   * @returns {string|null} Determined domain
   */
  determineDomain(text) {
    const textLower = text.toLowerCase();
    let bestMatch = null;
    let bestScore = 0;
    
    for (const expert of this.expertPatterns) {
      const score = expert.keywords.filter(keyword => 
        textLower.includes(keyword)
      ).length;
      
      if (score > bestScore) {
        bestScore = score;
        bestMatch = expert.domain;
      }
    }
    
    return bestScore >= 2 ? bestMatch : null;
  }

  /**
   * Extract expert persona from response
   * @param {string} response - Assistant response
   * @returns {string|null} Expert persona
   */
  extractExpertPersona(response) {
    // Look for persona introduction patterns
    const personaPatterns = [
      /(?:I'm|I am)\s+(?:Dr\.|Professor|Prof\.)?\s*([A-Z][a-z]+(?:\s+[A-Z][a-z]+)*),?\s+(?:a|an)?\s*([^.]+)/i,
      /As\s+(?:Dr\.|Professor|Prof\.)?\s*([A-Z][a-z]+(?:\s+[A-Z][a-z]+)*),?\s+(?:a|an)?\s*([^,]+)/i,
      /\*(?:Dr\.|Professor|Prof\.)?\s*([A-Z][a-z]+(?:\s+[A-Z][a-z]+)*)[^*]+\*/
    ];
    
    for (const pattern of personaPatterns) {
      const match = response.match(pattern);
      if (match) {
        return match[0].substring(0, 100); // Return the full matched persona
      }
    }
    
    return null;
  }

  /**
   * Check if query is incomplete without context
   * @param {string} query - User query
   * @returns {boolean} Whether query needs context
   */
  isIncompleteWithoutContext(query) {
    // Very short queries with pronouns
    if (query.length < 20 && /\b(it|this|that|they)\b/i.test(query)) {
      return true;
    }
    
    // Questions starting with certain patterns
    const incompletePatterns = [
      /^(and|but|or|so)\s+/i,
      /^(what|how|why)\s+about\b/i,
      /^(yes|no|okay|sure),?\s*$/i,
      /^\.\.\./,
      /^(more|less|better|worse)\s+than\s+what/i
    ];
    
    return incompletePatterns.some(pattern => pattern.test(query));
  }

  /**
   * Estimate confidence of response
   * @param {string} response - Assistant response
   * @returns {number} Confidence score (0-1)
   */
  estimateConfidence(response) {
    let confidence = 0.5;
    
    // Increase confidence for detailed responses
    if (response.length > 500) confidence += 0.1;
    if (response.length > 1000) confidence += 0.1;
    
    // Check for confidence indicators
    const highConfidenceTerms = [
      /\b(?:definitely|certainly|absolutely|clearly|obviously)\b/i,
      /\b(?:proven|established|confirmed|verified)\b/i,
      /\b(?:exactly|precisely|specifically)\b/i
    ];
    
    const lowConfidenceTerms = [
      /\b(?:might|maybe|perhaps|possibly|probably)\b/i,
      /\b(?:unclear|uncertain|unknown|debated)\b/i,
      /\b(?:approximately|roughly|about|around)\b/i
    ];
    
    highConfidenceTerms.forEach(term => {
      if (term.test(response)) confidence += 0.05;
    });
    
    lowConfidenceTerms.forEach(term => {
      if (term.test(response)) confidence -= 0.05;
    });
    
    // Check for citations or sources
    if (/\[\d+\]|\(\d{4}\)|https?:\/\//i.test(response)) {
      confidence += 0.1;
    }
    
    // Ensure confidence is within bounds
    return Math.max(0.1, Math.min(1.0, confidence));
  }

  /**
   * Determine if expert persona should change
   * @param {string} newQuery - New query
   * @param {string} previousDomain - Previous expert domain
   * @returns {Object} Switch recommendation
   */
  shouldSwitchExpert(newQuery, previousDomain) {
    const newDomain = this.determineDomain(newQuery);
    
    if (!newDomain || !previousDomain) {
      return {
        shouldSwitch: false,
        reason: 'insufficient-domain-signals'
      };
    }
    
    if (newDomain === previousDomain) {
      return {
        shouldSwitch: false,
        reason: 'same-domain'
      };
    }
    
    // Check if domains are related
    const relatedDomains = {
      'mathematics': ['physics', 'engineering', 'computer-science'],
      'physics': ['mathematics', 'engineering', 'chemistry'],
      'chemistry': ['biology', 'physics', 'medicine'],
      'biology': ['chemistry', 'medicine'],
      'computer-science': ['mathematics', 'engineering'],
      'medicine': ['biology', 'chemistry']
    };
    
    const isRelated = relatedDomains[previousDomain]?.includes(newDomain);
    
    return {
      shouldSwitch: true,
      newDomain,
      isRelated,
      reason: isRelated ? 'related-domain-shift' : 'domain-change'
    };
  }
}

// Export singleton instance
module.exports = new ContextAnalyzer();