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

// Create logs directory
const logsDir = path.join(__dirname, 'session-logs');
if (!fs.existsSync(logsDir)) {
  fs.mkdirSync(logsDir, { recursive: true });
}

// Current session log file
let currentLogFile = null;
let logStream = null;

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

// Broadcast to all connected clients AND save to file
function broadcast(data) {
  const message = JSON.stringify(data);
  
  // Send to all WebSocket clients
  clients.forEach(client => {
    if (client.readyState === WebSocket.OPEN) {
      client.send(message);
    }
  });
  
  // Write to log file
  if (logStream) {
    logStream.write(JSON.stringify({
      ...data,
      timestamp: new Date().toISOString()
    }) + '\n');
  }
}

// Start new log session
function startNewLogSession(questionId) {
  // Close previous log if exists
  if (logStream) {
    logStream.end();
  }
  
  // Create new log file
  const timestamp = new Date().toISOString().replace(/[:.]/g, '-');
  currentLogFile = path.join(logsDir, `session-${questionId}-${timestamp}.jsonl`);
  logStream = fs.createWriteStream(currentLogFile, { flags: 'a' });
  
  // Write session header
  logStream.write(JSON.stringify({
    type: 'session_start',
    questionId: questionId,
    timestamp: new Date().toISOString(),
    file: currentLogFile
  }) + '\n');
  
  // Copy to latest for easy access (Windows doesn't always support symlinks)
  const latestPath = path.join(logsDir, 'latest-session.jsonl');
  fs.writeFileSync(latestPath, `{"type":"redirect","actualFile":"${currentLogFile.replace(/\\/g, '/')}"}\n`);
  
  console.log(`📝 Log session started: ${currentLogFile}`);
  
  return currentLogFile;
}

// Enhanced logging interceptors
const originalLog = console.log;
const originalError = console.error;

console.log = function(...args) {
  originalLog.apply(console, args);
  const message = args.join(' ');
  
  // Detect special patterns
  if (message.includes('Math agent')) {
    broadcast({
      type: 'agent_start',
      agent: 'math',
      reason: 'Complex mathematical calculation detected',
      message: message
    });
  } else if (message.includes('Code agent')) {
    broadcast({
      type: 'agent_start',
      agent: 'code',
      reason: 'Code analysis or execution required',
      message: message
    });
  }
  
  broadcast({
    type: 'log',
    level: 'info',
    message: message
  });
};

console.error = function(...args) {
  originalError.apply(console, args);
  broadcast({
    type: 'log',
    level: 'error',
    message: args.join(' ')
  });
};

// Intercept fetch to track API calls
const originalFetch = global.fetch;
let apiCallCount = 0;

global.fetch = async function(...args) {
  const url = args[0];
  const options = args[1] || {};
  
  // Log all API calls
  broadcast({
    type: 'api_call_start',
    url: url,
    method: options.method || 'GET',
    hasBody: !!options.body
  });
  
  // Track specific services
  if (url && url.includes('openrouter.ai')) {
    apiCallCount++;
    const body = options.body ? JSON.parse(options.body) : {};
    broadcast({
      type: 'openrouter_call',
      model: body.model,
      messageCount: body.messages ? body.messages.length : 0,
      temperature: body.temperature,
      maxTokens: body.max_tokens
    });
  }
  
  // Track code execution
  if (url && url.includes('localhost:3001/run_code')) {
    const body = JSON.parse(options.body || '{}');
    broadcast({
      type: 'code_execution_request',
      language: body.language,
      codeLength: body.source ? body.source.length : 0,
      code: body.source,
      timeout: body.timeout
    });
  }
  
  // Track Serper
  if (url && url.includes('serper.dev')) {
    const body = JSON.parse(options.body || '{}');
    broadcast({
      type: 'search_request',
      query: body.q,
      num: body.num
    });
  }
  
  const startTime = Date.now();
  
  try {
    const response = await originalFetch.apply(this, args);
    const duration = Date.now() - startTime;
    
    // Clone to read response
    const cloned = response.clone();
    
    // Log response
    broadcast({
      type: 'api_call_complete',
      url: url,
      status: response.status,
      duration: duration,
      ok: response.ok
    });
    
    // Track specific responses
    if (url && url.includes('localhost:3001/run_code')) {
      try {
        const result = await cloned.json();
        broadcast({
          type: 'code_execution_result',
          stdout: result.stdout,
          stderr: result.stderr,
          exitCode: result.exitCode,
          executionTime: duration
        });
      } catch (e) {
        broadcast({
          type: 'code_execution_error',
          error: 'Failed to parse response',
          responseStatus: response.status
        });
      }
    }
    
    if (url && url.includes('serper.dev')) {
      try {
        const result = await cloned.json();
        broadcast({
          type: 'search_results',
          resultsCount: result.organic ? result.organic.length : 0,
          results: result.organic ? result.organic.slice(0, 3).map(r => ({
            title: r.title,
            link: r.link,
            snippet: r.snippet
          })) : []
        });
      } catch (e) {
        broadcast({
          type: 'search_error',
          error: 'Failed to parse search results'
        });
      }
    }
    
    if (url && url.includes('openrouter.ai')) {
      try {
        const result = await cloned.json();
        if (result.choices && result.choices[0]) {
          broadcast({
            type: 'openrouter_response',
            responseLength: result.choices[0].message.content.length,
            finishReason: result.choices[0].finish_reason,
            usage: result.usage
          });
        }
      } catch (e) {
        broadcast({
          type: 'openrouter_error',
          error: 'Failed to parse OpenRouter response',
          responseText: await cloned.text()
        });
      }
    }
    
    return response;
  } catch (error) {
    const duration = Date.now() - startTime;
    broadcast({
      type: 'api_call_error',
      url: url,
      error: error.message,
      duration: duration,
      stack: error.stack
    });
    throw error;
  }
};

// API endpoint to process questions
app.post('/api/ask', async (req, res) => {
  const { question, questionId = 'unknown' } = req.body;
  apiCallCount = 0;
  
  // Start new log session
  const logFile = startNewLogSession(questionId);
  
  broadcast({
    type: 'session_info',
    logFile: logFile,
    questionId: questionId
  });
  
  broadcast({
    type: 'question_start',
    questionId: questionId,
    questionLength: question.length,
    question: question
  });
  
  try {
    const startTime = Date.now();
    
    // Enhanced logger interception
    const logger = require('../src/utils/logger');
    const originalLoggerInfo = logger.info.bind(logger);
    
    logger.info = function(message, metadata) {
      originalLoggerInfo(message, metadata);
      
      // Log ALL metadata
      broadcast({
        type: 'workflow_log',
        message: message,
        metadata: metadata
      });
      
      if (metadata) {
        // Specific handling for different types
        if (metadata.step) {
          broadcast({
            type: 'workflow_step',
            step: metadata.step,
            requestId: metadata.requestId,
            fullMetadata: metadata
          });
        }
        
        if (metadata.type === 'search_decision') {
          broadcast({
            type: 'search_decision',
            needsSearch: metadata.needsSearch,
            rationale: metadata.rationale,
            analysisMethod: metadata.analysisMethod,
            query: metadata.query
          });
        }
        
        if (metadata.queryType) {
          broadcast({
            type: 'query_classification',
            queryType: metadata.queryType,
            complexity: metadata.complexity,
            codingType: metadata.codingType,
            requestId: metadata.requestId
          });
        }
        
        if (metadata.error) {
          broadcast({
            type: 'workflow_error',
            error: metadata.error,
            stack: metadata.stack,
            duration: metadata.duration,
            requestId: metadata.requestId
          });
        }
      }
    };
    
    // Process question
    const response = await workflowEngine.answer(question, { 
      timeout: 300000,
      onProgress: (progress) => {
        broadcast({
          type: 'progress_update',
          ...progress
        });
      }
    });
    
    const duration = Date.now() - startTime;
    
    // Final response
    broadcast({
      type: 'question_complete',
      questionId: questionId,
      duration: duration,
      responseLength: response.content.length,
      searchUsed: response.searchPerformed,
      codeExecuted: response.codeExecuted,
      apiCalls: apiCallCount,
      response: response.content
    });
    
    // Close log session
    broadcast({
      type: 'session_end',
      questionId: questionId,
      success: true,
      totalDuration: duration
    });
    
    if (logStream) {
      logStream.end();
      console.log(`📝 Log session saved: ${currentLogFile}`);
    }
    
    res.json({
      success: true,
      response: response.content,
      duration: duration,
      logFile: currentLogFile,
      metadata: {
        searchUsed: response.searchPerformed,
        codeExecuted: response.codeExecuted,
        apiCalls: apiCallCount
      }
    });
    
  } catch (error) {
    broadcast({
      type: 'question_error',
      questionId: questionId,
      error: error.message,
      stack: error.stack
    });
    
    broadcast({
      type: 'session_end',
      questionId: questionId,
      success: false,
      error: error.message
    });
    
    if (logStream) {
      logStream.end();
    }
    
    res.status(500).json({
      success: false,
      error: error.message,
      logFile: currentLogFile
    });
  }
});

// Endpoint to get latest log
app.get('/api/latest-log', (req, res) => {
  const latestPath = path.join(logsDir, 'latest-session.jsonl');
  if (fs.existsSync(latestPath)) {
    res.sendFile(latestPath);
  } else {
    res.status(404).json({ error: 'No log file found' });
  }
});

// Endpoint to list all logs
app.get('/api/logs', (req, res) => {
  const files = fs.readdirSync(logsDir)
    .filter(f => f.endsWith('.jsonl') && f !== 'latest-session.jsonl')
    .map(f => ({
      name: f,
      path: path.join(logsDir, f),
      size: fs.statSync(path.join(logsDir, f)).size,
      created: fs.statSync(path.join(logsDir, f)).birthtime
    }))
    .sort((a, b) => b.created - a.created);
  
  res.json({ logs: files });
});

// Serve the advanced UI
app.use(express.static(path.join(__dirname, 'public')));

// Start server
const PORT = process.env.UI_PORT || 3004;
server.listen(PORT, () => {
  console.log(`📊 Advanced UI with logging running on http://localhost:${PORT}`);
  console.log(`📝 Session logs will be saved to: ${logsDir}`);
  console.log(`🔗 Latest log always available at: ${path.join(logsDir, 'latest-session.jsonl')}`);
  console.log(`✅ E2B service should be running on port 3001`);
});