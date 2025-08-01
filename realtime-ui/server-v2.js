const express = require('express');
const http = require('http');
const WebSocket = require('ws');
const path = require('path');
const fs = require('fs');

// Import the workflow engine
const workflowEngine = require('../src/workflow-engine');

const app = express();
const server = http.createServer(app);
const wss = new WebSocket.Server({ server });

// Serve static files
app.use(express.static(path.join(__dirname, 'public')));
app.use(express.json());

// Store active connections
const clients = new Set();

// WebSocket connection handler
wss.on('connection', (ws) => {
  clients.add(ws);
  console.log('Client connected. Total clients:', clients.size);
  
  ws.on('close', () => {
    clients.delete(ws);
    console.log('Client disconnected. Total clients:', clients.size);
  });
});

// Broadcast to all connected clients
function broadcast(data) {
  const message = JSON.stringify(data);
  clients.forEach(client => {
    if (client.readyState === WebSocket.OPEN) {
      client.send(message);
    }
  });
}

// Enhanced logging interceptors
const originalLog = console.log;
const originalError = console.error;
const originalInfo = console.info;

console.log = function(...args) {
  originalLog.apply(console, args);
  const message = args.join(' ');
  
  // Detect special patterns in logs
  if (message.includes('Math agent')) {
    broadcast({
      type: 'agent_start',
      agent: 'math',
      reason: 'Complex mathematical calculation detected',
      timestamp: new Date().toISOString()
    });
  } else if (message.includes('Code agent')) {
    broadcast({
      type: 'agent_start',
      agent: 'code',
      reason: 'Code analysis or execution required',
      timestamp: new Date().toISOString()
    });
  }
  
  broadcast({
    type: 'log',
    level: 'info',
    timestamp: new Date().toISOString(),
    message: message
  });
};

console.error = function(...args) {
  originalError.apply(console, args);
  broadcast({
    type: 'log',
    level: 'error',
    timestamp: new Date().toISOString(),
    message: args.join(' ')
  });
};

// Intercept fetch to track API calls
const originalFetch = global.fetch;
let apiCallCount = 0;

global.fetch = async function(...args) {
  const url = args[0];
  const options = args[1] || {};
  
  // Track OpenRouter API calls
  if (url && url.includes('openrouter.ai')) {
    apiCallCount++;
    broadcast({
      type: 'api_call',
      service: 'OpenRouter',
      count: apiCallCount,
      timestamp: new Date().toISOString()
    });
  }
  
  // Track code execution calls
  if (url && url.includes('localhost:3001/run_code')) {
    const body = JSON.parse(options.body || '{}');
    broadcast({
      type: 'code_execution',
      language: body.language,
      code: body.source,
      timestamp: new Date().toISOString()
    });
  }
  
  // Track Serper API calls
  if (url && url.includes('serper.dev')) {
    const body = JSON.parse(options.body || '{}');
    broadcast({
      type: 'search_query',
      query: body.q,
      timestamp: new Date().toISOString()
    });
  }
  
  try {
    const response = await originalFetch.apply(this, args);
    
    // Clone response to read it
    const cloned = response.clone();
    
    // Track code execution results
    if (url && url.includes('localhost:3001/run_code')) {
      try {
        const result = await cloned.json();
        broadcast({
          type: 'code_execution',
          output: result.stdout,
          error: result.stderr,
          exitCode: result.exitCode,
          timestamp: new Date().toISOString()
        });
      } catch (e) {
        // Ignore JSON parsing errors
      }
    }
    
    // Track search results
    if (url && url.includes('serper.dev')) {
      try {
        const result = await cloned.json();
        broadcast({
          type: 'search_query',
          results: result.organic ? result.organic.slice(0, 5) : [],
          timestamp: new Date().toISOString()
        });
      } catch (e) {
        // Ignore JSON parsing errors
      }
    }
    
    return response;
  } catch (error) {
    broadcast({
      type: 'api_error',
      service: url,
      error: error.message,
      timestamp: new Date().toISOString()
    });
    throw error;
  }
};

// API endpoint to process questions
app.post('/api/ask', async (req, res) => {
  const { question } = req.body;
  apiCallCount = 0; // Reset counter
  
  broadcast({
    type: 'status',
    message: 'Starting to process question...',
    timestamp: new Date().toISOString()
  });
  
  broadcast({
    type: 'question',
    content: question,
    timestamp: new Date().toISOString()
  });
  
  try {
    const startTime = Date.now();
    
    // Intercept logger for more detailed workflow tracking
    const logger = require('../src/utils/logger');
    const originalLoggerInfo = logger.info.bind(logger);
    
    logger.info = function(message, metadata) {
      originalLoggerInfo(message, metadata);
      
      if (metadata) {
        // Enhanced step tracking
        if (metadata.step) {
          broadcast({
            type: 'step',
            step: metadata.step,
            metadata: metadata,
            timestamp: new Date().toISOString()
          });
        }
        
        // Track search decisions
        if (metadata.type === 'search_decision') {
          broadcast({
            type: 'decision',
            needsSearch: metadata.needsSearch,
            rationale: metadata.rationale,
            analysisMethod: metadata.analysisMethod,
            timestamp: new Date().toISOString()
          });
        }
        
        // Track query types
        if (metadata.queryType) {
          broadcast({
            type: 'query_type',
            queryType: metadata.queryType,
            complexity: metadata.complexity,
            codingType: metadata.codingType,
            timestamp: new Date().toISOString()
          });
        }
        
        // Track search planning
        if (metadata.type === 'search_plan') {
          broadcast({
            type: 'search_plan',
            queries: metadata.queries,
            timestamp: new Date().toISOString()
          });
        }
        
        // Track synthesis
        if (metadata.type === 'synthesis_start') {
          broadcast({
            type: 'synthesis',
            sourceCount: metadata.sourceCount,
            timestamp: new Date().toISOString()
          });
        }
      }
    };
    
    // Process the question with progress callback
    const response = await workflowEngine.answer(question, { 
      timeout: 300000,
      onProgress: (progress) => {
        broadcast({
          type: 'progress',
          ...progress,
          timestamp: new Date().toISOString()
        });
      }
    });
    
    const duration = Date.now() - startTime;
    
    // Send final response
    broadcast({
      type: 'response',
      content: response.content,
      duration: duration,
      searchUsed: response.searchPerformed,
      codeExecuted: response.codeExecuted,
      apiCalls: apiCallCount,
      timestamp: new Date().toISOString()
    });
    
    res.json({
      success: true,
      response: response.content,
      duration: duration,
      metadata: {
        searchUsed: response.searchPerformed,
        codeExecuted: response.codeExecuted,
        apiCalls: apiCallCount
      }
    });
    
  } catch (error) {
    broadcast({
      type: 'error',
      message: error.message,
      stack: error.stack,
      timestamp: new Date().toISOString()
    });
    
    res.status(500).json({
      success: false,
      error: error.message
    });
  }
});

// Default route serves the advanced UI
app.get('/', (req, res) => {
  res.sendFile(path.join(__dirname, 'public', 'index-v2.html'));
});

// Start server
const PORT = process.env.UI_PORT || 3003;
server.listen(PORT, () => {
  console.log(`Advanced Realtime UI server running on http://localhost:${PORT}`);
  console.log(`WebSocket server ready for connections`);
  console.log(`E2B service should be running on port 3001`);
});