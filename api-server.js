const express = require("express");
const cors = require("cors");
// Use modular workflow engine if configured
const WorkflowEngine = require("./src/workflow-engine-modular");
const workflowEngine = new WorkflowEngine(require('./config'));
const sessionStore = require("./src/session-store");
const contextAnalyzer = require("./src/context-analyzer");
const logger = require("./src/utils/logger");

const app = express();
const PORT = process.env.PORT || 8080;

app.use(cors());
app.use(express.json());

// Health check
app.get("/health", (req, res) => {
  res.json({ status: "healthy", service: "Ask Any Expert API" });
});

// API Health check (for frontend compatibility)
app.get("/api/health", (req, res) => {
  res.json({ status: "healthy", service: "Ask Any Expert API" });
});

// Main ask endpoint with session support
app.post("/api/ask", async (req, res) => {
  const { question, sessionId } = req.body;
  
  if (!question) {
    return res.status(400).json({ error: "Question is required" });
  }
  
  // Get or create session
  const session = sessionStore.getOrCreate(sessionId);
  
  logger.info({ 
    question, 
    sessionId: session.id,
    isNewSession: !sessionId || sessionId !== session.id
  }, "Processing API query");
  
  try {
    // Analyze the query for context dependencies
    const queryAnalysis = contextAnalyzer.analyzeQuery(question, session);
    
    // Get conversation history
    const history = sessionStore.getHistory(session.id);
    
    // Process the query with context
    const result = await workflowEngine.answer(question, history);
    
    // Analyze the response for metadata extraction
    const responseMetadata = contextAnalyzer.analyzeResponse(result.content, question);
    
    // Add messages to session
    sessionStore.addMessage(session.id, 'user', question, {
      ...queryAnalysis,
      timestamp: Date.now()
    });
    
    sessionStore.addMessage(session.id, 'assistant', result.content, {
      ...responseMetadata,
      toolsUsed: result.metadata?.toolsUsed || [],
      sources: result.sources || [],
      timestamp: Date.now()
    });
    
    // Prepare response with session info
    res.json({
      answer: result.content,
      sources: result.sources,
      searchPerformed: result.searchPerformed,
      agent: result.metadata?.agent || "unknown",
      confidence: responseMetadata.confidence || result.metadata?.finalConfidence,
      metadata: {
        ...result.metadata,
        contextUsed: queryAnalysis.contextDependency !== 'low',
        isFollowUp: queryAnalysis.isFollowUp,
        expertPersona: responseMetadata.expertPersona,
        sessionInfo: {
          sessionId: session.id,
          messageCount: session.metadata.messageCount,
          topics: Array.from(session.metadata.topics || []),
          expertsUsed: Array.from(session.metadata.expertsUsed || [])
        }
      },
      sessionId: session.id
    });
    
  } catch (error) {
    logger.error({ error: error.message, sessionId: session.id }, "Query processing failed");
    res.status(500).json({ 
      error: "Failed to process question",
      message: error.message,
      sessionId: session.id
    });
  }
});

// Session info endpoint
app.get("/api/session/:sessionId", (req, res) => {
  const { sessionId } = req.params;
  const history = sessionStore.getHistory(sessionId);
  
  if (!history || history.messages.length === 0) {
    return res.status(404).json({ error: "Session not found" });
  }
  
  res.json({
    sessionId,
    messageCount: history.messages.length,
    created: history.created,
    lastAccessed: history.lastAccessed,
    topics: history.metadata.topics,
    entities: history.metadata.entities,
    expertsUsed: history.metadata.expertsUsed
  });
});

// Clear session endpoint
app.delete("/api/session/:sessionId", (req, res) => {
  const { sessionId } = req.params;
  const deleted = sessionStore.delete(sessionId);
  
  res.json({
    success: deleted,
    message: deleted ? "Session cleared" : "Session not found"
  });
});

app.listen(PORT, () => {
  console.log(`Ask Any Expert API running on http://localhost:${PORT}`);
  console.log(`Health check: http://localhost:${PORT}/health`);
});