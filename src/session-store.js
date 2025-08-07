const crypto = require('crypto');
const logger = require('./utils/logger');

/**
 * SessionStore - Manages conversation sessions with TTL and cleanup
 * Simple in-memory implementation, easily upgradeable to Redis
 */
class SessionStore {
  constructor(options = {}) {
    this.sessions = new Map();
    this.ttl = options.ttl || parseInt(process.env.SESSION_TTL || '3600000'); // 1 hour default
    this.maxSessions = options.maxSessions || 1000;
    this.cleanupInterval = options.cleanupInterval || 300000; // 5 minutes
    
    // Start cleanup timer
    this.startCleanupTimer();
    
    logger.info({
      message: 'SessionStore initialized',
      ttl: this.ttl,
      maxSessions: this.maxSessions
    });
  }

  /**
   * Generate a unique session ID
   * @returns {string} Session ID
   */
  generateId() {
    return `session-${Date.now()}-${crypto.randomBytes(8).toString('hex')}`;
  }

  /**
   * Get or create a session
   * @param {string} sessionId - Optional session ID
   * @returns {Object} Session object
   */
  getOrCreate(sessionId) {
    // If no sessionId provided or session doesn't exist, create new
    if (!sessionId || !this.sessions.has(sessionId)) {
      return this.create();
    }
    
    const session = this.sessions.get(sessionId);
    
    // Check if session expired
    if (this.isExpired(session)) {
      this.delete(sessionId);
      return this.create();
    }
    
    // Update last accessed time
    session.lastAccessed = Date.now();
    return session;
  }

  /**
   * Create a new session
   * @returns {Object} New session object
   */
  create() {
    // Check if we've reached max sessions
    if (this.sessions.size >= this.maxSessions) {
      this.cleanupOldest();
    }
    
    const sessionId = this.generateId();
    const session = {
      id: sessionId,
      messages: [],
      created: Date.now(),
      lastAccessed: Date.now(),
      metadata: {
        messageCount: 0,
        topics: new Set(),
        entities: new Set(),
        expertsUsed: new Set()
      }
    };
    
    this.sessions.set(sessionId, session);
    
    logger.info({
      message: 'New session created',
      sessionId,
      totalSessions: this.sessions.size
    });
    
    return session;
  }

  /**
   * Get a session by ID
   * @param {string} sessionId - Session ID
   * @returns {Object|null} Session object or null if not found/expired
   */
  get(sessionId) {
    if (!sessionId || !this.sessions.has(sessionId)) {
      return null;
    }
    
    const session = this.sessions.get(sessionId);
    
    if (this.isExpired(session)) {
      this.delete(sessionId);
      return null;
    }
    
    session.lastAccessed = Date.now();
    return session;
  }

  /**
   * Add a message to a session
   * @param {string} sessionId - Session ID
   * @param {string} role - 'user' or 'assistant'
   * @param {string} content - Message content
   * @param {Object} metadata - Additional metadata
   */
  addMessage(sessionId, role, content, metadata = {}) {
    const session = this.get(sessionId);
    if (!session) {
      logger.warn({ sessionId, message: 'Attempted to add message to non-existent session' });
      return false;
    }
    
    const message = {
      role,
      content,
      timestamp: Date.now(),
      metadata
    };
    
    session.messages.push(message);
    session.metadata.messageCount++;
    
    // Update session metadata
    if (metadata.topics) {
      metadata.topics.forEach(topic => session.metadata.topics.add(topic));
    }
    if (metadata.entities) {
      metadata.entities.forEach(entity => session.metadata.entities.add(entity));
    }
    if (metadata.expertPersona) {
      session.metadata.expertsUsed.add(metadata.expertPersona);
    }
    
    return true;
  }

  /**
   * Get conversation history for a session
   * @param {string} sessionId - Session ID
   * @returns {Object} Formatted history object
   */
  getHistory(sessionId) {
    const session = this.get(sessionId);
    if (!session) {
      return {
        messages: [],
        metadata: {}
      };
    }
    
    return {
      sessionId: session.id,
      messages: session.messages,
      metadata: {
        ...session.metadata,
        topics: Array.from(session.metadata.topics),
        entities: Array.from(session.metadata.entities),
        expertsUsed: Array.from(session.metadata.expertsUsed)
      },
      created: session.created,
      lastAccessed: session.lastAccessed
    };
  }

  /**
   * Delete a session
   * @param {string} sessionId - Session ID
   */
  delete(sessionId) {
    const deleted = this.sessions.delete(sessionId);
    if (deleted) {
      logger.info({
        message: 'Session deleted',
        sessionId,
        remainingSessions: this.sessions.size
      });
    }
    return deleted;
  }

  /**
   * Check if a session is expired
   * @param {Object} session - Session object
   * @returns {boolean} True if expired
   */
  isExpired(session) {
    return Date.now() - session.lastAccessed > this.ttl;
  }

  /**
   * Clean up expired sessions
   */
  cleanup() {
    let cleaned = 0;
    for (const [sessionId, session] of this.sessions.entries()) {
      if (this.isExpired(session)) {
        this.sessions.delete(sessionId);
        cleaned++;
      }
    }
    
    if (cleaned > 0) {
      logger.info({
        message: 'Session cleanup completed',
        sessionsRemoved: cleaned,
        remainingSessions: this.sessions.size
      });
    }
  }

  /**
   * Remove oldest sessions when at capacity
   */
  cleanupOldest() {
    const sortedSessions = Array.from(this.sessions.entries())
      .sort((a, b) => a[1].lastAccessed - b[1].lastAccessed);
    
    // Remove oldest 10% of sessions
    const toRemove = Math.ceil(this.sessions.size * 0.1);
    for (let i = 0; i < toRemove && i < sortedSessions.length; i++) {
      this.sessions.delete(sortedSessions[i][0]);
    }
    
    logger.info({
      message: 'Removed oldest sessions',
      removed: toRemove,
      remainingSessions: this.sessions.size
    });
  }

  /**
   * Start the cleanup timer
   */
  startCleanupTimer() {
    this.cleanupTimer = setInterval(() => {
      this.cleanup();
    }, this.cleanupInterval);
  }

  /**
   * Stop the cleanup timer (for graceful shutdown)
   */
  stopCleanupTimer() {
    if (this.cleanupTimer) {
      clearInterval(this.cleanupTimer);
      this.cleanupTimer = null;
    }
  }

  /**
   * Get statistics about sessions
   * @returns {Object} Session statistics
   */
  getStats() {
    const now = Date.now();
    const sessions = Array.from(this.sessions.values());
    
    return {
      totalSessions: this.sessions.size,
      activeSessions: sessions.filter(s => now - s.lastAccessed < 300000).length, // Active in last 5 min
      oldestSession: sessions.reduce((oldest, s) => 
        !oldest || s.created < oldest.created ? s : oldest, null),
      averageMessageCount: sessions.reduce((sum, s) => 
        sum + s.metadata.messageCount, 0) / (sessions.length || 1)
    };
  }
}

// Export singleton instance
module.exports = new SessionStore();