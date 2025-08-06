const express = require("express");
const cors = require("cors");
const workflowEngine = require("./src/workflow-engine");
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

// E2B diagnostic endpoint
app.get("/api/debug/e2b", async (req, res) => {
  const e2bManager = require("./src/e2b-manager-v3");
  const config = require("./config");
  const agent = require("./src/unified-agent-modular");
  
  const diagnostics = {
    environment: {
      E2B_API_KEY: process.env.E2B_API_KEY ? "Set (hidden)" : "NOT SET",
      E2B_TEMPLATE_ID: process.env.E2B_TEMPLATE_ID || "Not set (using default)",
      API_MODE: process.env.API_MODE || "full",
      NODE_ENV: process.env.NODE_ENV || "development",
      OPENROUTER_API_KEY: process.env.OPENROUTER_API_KEY ? "Set (hidden)" : "NOT SET",
      OPENROUTER_MODEL: process.env.OPENROUTER_MODEL || "Not set",
      MAX_ROUNDS: process.env.MAX_ROUNDS || "Not set",
      CONF_THRESHOLD: process.env.CONF_THRESHOLD || "Not set",
      actualMaxSteps: agent.maxSteps,
      actualConfThreshold: agent.confidenceThreshold
    },
    config: {
      templateId: e2bManager.templateId,
      maxRetries: e2bManager.maxRetries,
      metrics: e2bManager.getMetrics()
    },
    test: {
      simpleCalculation: null,
      error: null
    }
  };
  
  // Test E2B if API key is set
  if (process.env.E2B_API_KEY) {
    try {
      const result = await e2bManager.executeCode("print(2 + 2)", {
        timeout: 5000,
        requestId: "diagnostic-test"
      });
      diagnostics.test.simpleCalculation = result;
    } catch (error) {
      diagnostics.test.error = error.message;
    }
  } else {
    diagnostics.test.error = "E2B_API_KEY not configured";
  }
  
  res.json(diagnostics);
});

// Test endpoint for debugging tool calls
app.post("/api/test-tools", async (req, res) => {
  const { question } = req.body;
  const config = require("./config");
  
  try {
    // Direct OpenRouter call with tools
    const response = await fetch(`${config.openrouter.baseUrl}/chat/completions`, {
      method: 'POST',
      headers: {
        'Authorization': `Bearer ${config.openrouter.apiKey}`,
        'Content-Type': 'application/json',
        'HTTP-Referer': 'https://askanyexpert.ai',
        'X-Title': 'Ask Any Expert'
      },
      body: JSON.stringify({
        model: config.openrouter.model,
        messages: [
          {
            role: 'system',
            content: 'You are a helpful assistant. Use the run_code tool when calculations are needed.'
          },
          {
            role: 'user',
            content: question || 'What is 2+2?'
          }
        ],
        tools: [{
          type: 'function',
          function: {
            name: 'run_code',
            description: 'Execute Python code',
            parameters: {
              type: 'object',
              properties: {
                code: { type: 'string', description: 'Python code to execute' }
              },
              required: ['code']
            }
          }
        }],
        tool_choice: 'auto',
        max_tokens: 1000
      })
    });
    
    const data = await response.json();
    
    res.json({
      model: config.openrouter.model,
      hasToolCalls: !!(data.choices?.[0]?.message?.tool_calls),
      toolCallsCount: data.choices?.[0]?.message?.tool_calls?.length || 0,
      content: data.choices?.[0]?.message?.content,
      toolCalls: data.choices?.[0]?.message?.tool_calls,
      raw: data
    });
    
  } catch (error) {
    res.status(500).json({ error: error.message });
  }
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
