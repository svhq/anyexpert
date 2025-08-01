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

// Override console.log to capture logs
const originalLog = console.log;
const originalError = console.error;
const originalInfo = console.info;

console.log = function(...args) {
  originalLog.apply(console, args);
  broadcast({
    type: 'log',
    level: 'info',
    timestamp: new Date().toISOString(),
    message: args.join(' ')
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

// API endpoint to process questions
app.post('/api/ask', async (req, res) => {
  const { question } = req.body;
  
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
    // Track steps
    let steps = [];
    const startTime = Date.now();
    
    // Intercept logger to capture workflow steps
    const logger = require('../src/utils/logger');
    const originalLoggerInfo = logger.info.bind(logger);
    logger.info = function(message, metadata) {
      originalLoggerInfo(message, metadata);
      
      // Broadcast workflow steps
      if (metadata) {
        if (metadata.step) {
          broadcast({
            type: 'step',
            step: metadata.step,
            metadata: metadata,
            timestamp: new Date().toISOString()
          });
        }
        
        if (metadata.type === 'search_decision') {
          broadcast({
            type: 'decision',
            needsSearch: metadata.needsSearch,
            rationale: metadata.rationale,
            timestamp: new Date().toISOString()
          });
        }
        
        if (metadata.queryType) {
          broadcast({
            type: 'query_type',
            queryType: metadata.queryType,
            complexity: metadata.complexity,
            timestamp: new Date().toISOString()
          });
        }
      }
    };
    
    // Process the question
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
      timestamp: new Date().toISOString()
    });
    
    res.json({
      success: true,
      response: response.content,
      duration: duration,
      metadata: {
        searchUsed: response.searchPerformed,
        codeExecuted: response.codeExecuted
      }
    });
    
  } catch (error) {
    broadcast({
      type: 'error',
      message: error.message,
      timestamp: new Date().toISOString()
    });
    
    res.status(500).json({
      success: false,
      error: error.message
    });
  }
});

// Start server
const PORT = process.env.UI_PORT || 3002;
server.listen(PORT, () => {
  console.log(`Realtime UI server running on http://localhost:${PORT}`);
  console.log(`WebSocket server ready for connections`);
});