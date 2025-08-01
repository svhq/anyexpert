const express = require('express');
const http = require('http');
const WebSocket = require('ws');
const path = require('path');
const fs = require('fs');

// Load environment variables from parent directory
require('dotenv').config({ path: path.join(__dirname, '..', '.env') });

const app = express();
const server = http.createServer(app);
const wss = new WebSocket.Server({ server });

// Create logs directory
const logsDir = path.join(__dirname, 'session-logs');
if (!fs.existsSync(logsDir)) {
  fs.mkdirSync(logsDir, { recursive: true });
}

// Current log file
let logFile = null;
let logStream = null;

// Store WebSocket clients
const clients = new Set();

// Middleware
app.use(express.json());
app.use(express.static(path.join(__dirname, 'public')));

// WebSocket connection
wss.on('connection', (ws) => {
  clients.add(ws);
  console.log('Client connected. Total:', clients.size);
  
  ws.on('close', () => {
    clients.delete(ws);
    console.log('Client disconnected. Total:', clients.size);
  });
  
  ws.on('error', (err) => {
    console.error('WebSocket error:', err);
  });
});

// Broadcast and log
function broadcast(data) {
  const entry = {
    ...data,
    timestamp: new Date().toISOString()
  };
  
  // Send to clients
  const message = JSON.stringify(entry);
  clients.forEach(client => {
    if (client.readyState === WebSocket.OPEN) {
      client.send(message);
    }
  });
  
  // Write to file
  if (logStream) {
    logStream.write(JSON.stringify(entry) + '\n');
  }
}

// Start logging session
function startLogging(questionId) {
  if (logStream) {
    logStream.end();
  }
  
  const timestamp = new Date().toISOString().replace(/[:.]/g, '-');
  logFile = path.join(logsDir, `session-${timestamp}.jsonl`);
  logStream = fs.createWriteStream(logFile);
  
  // Write header
  logStream.write(JSON.stringify({
    type: 'session_start',
    questionId: questionId,
    timestamp: new Date().toISOString()
  }) + '\n');
  
  // Update latest pointer
  fs.writeFileSync(
    path.join(logsDir, 'latest-session.jsonl'),
    JSON.stringify({ redirect: logFile }) + '\n'
  );
  
  console.log('Started logging to:', logFile);
  return logFile;
}

// Simple proxy to workflow engine
app.post('/api/ask', async (req, res) => {
  const { question } = req.body;
  const questionId = `q-${Date.now()}`;
  
  try {
    // Start logging
    const logPath = startLogging(questionId);
    
    broadcast({
      type: 'question',
      questionId: questionId,
      question: question
    });
    
    // Import workflow engine here to avoid startup issues
    broadcast({ type: 'log', level: 'info', message: 'Loading workflow engine...' });
    const workflowEngine = require('../src/workflow-engine');
    broadcast({ type: 'log', level: 'info', message: 'Workflow engine loaded' });
    
    // Intercept logger BEFORE using workflow engine
    const logger = require('../src/utils/logger');
    const originalLoggerInfo = logger.info.bind(logger);
    logger.info = function(message, metadata) {
      originalLoggerInfo(message, metadata);
      broadcast({
        type: 'log',
        level: 'info',
        message: `[LOGGER] ${message}`,
        metadata: metadata
      });
    };
    
    // Intercept console logs
    const originalLog = console.log;
    const originalError = console.error;
    
    console.log = (...args) => {
      originalLog(...args);
      broadcast({
        type: 'log',
        level: 'info',
        message: args.join(' ')
      });
    };
    
    console.error = (...args) => {
      originalError(...args);
      broadcast({
        type: 'log',
        level: 'error',
        message: args.join(' ')
      });
    };
    
    // Process question
    const startTime = Date.now();
    broadcast({ type: 'log', level: 'info', message: 'Calling workflow engine...' });
    const response = await workflowEngine.answer(question, { timeout: 300000 });
    const duration = Date.now() - startTime;
    broadcast({ type: 'log', level: 'info', message: `Response received in ${duration}ms` });
    
    // Restore console
    console.log = originalLog;
    console.error = originalError;
    
    // Log completion
    broadcast({
      type: 'response',
      content: response.content,
      duration: duration
    });
    
    broadcast({
      type: 'session_end',
      success: true
    });
    
    if (logStream) {
      logStream.end();
    }
    
    res.json({
      success: true,
      response: response.content,
      duration: duration,
      logFile: logPath
    });
    
  } catch (error) {
    console.error('Error:', error);
    
    broadcast({
      type: 'error',
      error: error.message
    });
    
    if (logStream) {
      logStream.end();
    }
    
    res.status(500).json({
      success: false,
      error: error.message
    });
  }
});

// Health check
app.get('/api/health', (req, res) => {
  res.json({ 
    status: 'ok', 
    clients: clients.size,
    logsDir: logsDir
  });
});

// Start server
const PORT = 3005;
server.listen(PORT, () => {
  console.log(`Simple logging server running on http://localhost:${PORT}`);
  console.log(`Logs saved to: ${logsDir}`);
});