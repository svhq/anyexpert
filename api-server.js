const express = require("express");
const cors = require("cors");
const workflowEngine = require("./src/workflow-engine");
const logger = require("./src/utils/logger");

const app = express();
const PORT = 3000;

app.use(cors());
app.use(express.json());

// Health check
app.get("/health", (req, res) => {
  res.json({ status: "healthy", service: "Ask Any Expert API" });
});

// Main ask endpoint
app.post("/api/ask", async (req, res) => {
  const { question } = req.body;
  
  if (!question) {
    return res.status(400).json({ error: "Question is required" });
  }
  
  logger.info({ question }, "Processing API query");
  
  try {
    const result = await workflowEngine.answer(question);
    
    res.json({
      answer: result.content,
      sources: result.sources,
      searchPerformed: result.searchPerformed,
      agent: result.metadata?.agent || "unknown",
      confidence: result.metadata?.finalConfidence,
      metadata: result.metadata
    });
    
  } catch (error) {
    logger.error({ error: error.message }, "Query processing failed");
    res.status(500).json({ 
      error: "Failed to process question",
      message: error.message 
    });
  }
});

app.listen(PORT, () => {
  console.log(`Ask Any Expert API running on http://localhost:${PORT}`);
  console.log(`Health check: http://localhost:${PORT}/health`);
});
