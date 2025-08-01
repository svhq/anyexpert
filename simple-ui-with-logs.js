const express = require('express');
const fs = require('fs');
const path = require('path');
const workflowEngine = require('./src/workflow-engine');

const app = express();
app.use(express.json());
app.use(express.static('public')); // Serve static files

// Create logs directory
const logsDir = path.join(__dirname, 'simple-logs');
if (!fs.existsSync(logsDir)) {
  fs.mkdirSync(logsDir, { recursive: true });
}

let currentLogFile = null;
let logs = [];

// Function to log messages
function log(type, message, metadata = {}) {
  const entry = {
    timestamp: new Date().toISOString(),
    type,
    message,
    ...metadata
  };
  
  logs.push(entry);
  
  // Write to file
  if (currentLogFile) {
    fs.appendFileSync(currentLogFile, JSON.stringify(entry) + '\n');
  }
  
  // Also console log
  console.log(`[${type}] ${message}`);
}

// Simple HTML UI with GPQA questions
app.get('/', (req, res) => {
  res.send(`
    <!DOCTYPE html>
    <html>
    <head>
      <title>Simple Ask UI with Logs</title>
      <style>
        body { font-family: monospace; padding: 20px; background: #1a1a1a; color: #e0e0e0; }
        .container { display: flex; gap: 20px; }
        .left { flex: 1; }
        .right { flex: 1; }
        select, button { width: 100%; padding: 10px; margin: 5px 0; }
        button { background: #00ff88; color: #000; font-weight: bold; cursor: pointer; }
        button:disabled { background: #666; cursor: not-allowed; }
        #output { white-space: pre-wrap; background: #2a2a2a; padding: 10px; margin: 10px 0; min-height: 300px; }
        #logs { background: #0a0a0a; padding: 10px; height: 400px; overflow-y: auto; font-size: 12px; }
        .log-entry { margin: 2px 0; padding: 2px; }
        .log-info { color: #4ecdc4; }
        .log-error { color: #ff6b6b; }
        .log-success { color: #00ff88; }
        h3 { color: #00ff88; }
      </style>
    </head>
    <body>
      <h1>Simple GPQA Tester with Logs</h1>
      
      <div class="container">
        <div class="left">
          <h3>Question</h3>
          <select id="questionSelect">
            <option value="">Select a GPQA question...</option>
          </select>
          <button id="askBtn" onclick="ask()">Ask Question</button>
          
          <h3>Response</h3>
          <div id="output">No response yet</div>
          
          <p style="margin-top: 20px; font-size: 12px;">
            Log file: <span id="logFile">None</span>
          </p>
        </div>
        
        <div class="right">
          <h3>Live Logs</h3>
          <button onclick="clearLogs()" style="width: auto; padding: 5px 10px; font-size: 12px;">Clear</button>
          <div id="logs"></div>
        </div>
      </div>
      
      <script>
        let questions = [];
        let polling = false;
        
        // Load GPQA questions
        fetch('/gpqa-questions.json')
          .then(r => r.json())
          .then(data => {
            questions = data.questions;
            const select = document.getElementById('questionSelect');
            questions.forEach((q, i) => {
              const opt = document.createElement('option');
              opt.value = i;
              opt.text = q.id + ': ' + q.question.substring(0, 60) + '...';
              select.add(opt);
            });
          });
        
        // Poll for logs
        async function pollLogs() {
          if (!polling) return;
          
          try {
            const response = await fetch('/logs');
            const data = await response.json();
            
            const logsDiv = document.getElementById('logs');
            logsDiv.innerHTML = '';
            
            data.logs.forEach(log => {
              const div = document.createElement('div');
              div.className = 'log-entry log-' + log.type;
              const time = new Date(log.timestamp).toLocaleTimeString();
              div.textContent = time + ' [' + log.type + '] ' + log.message;
              logsDiv.appendChild(div);
            });
            
            logsDiv.scrollTop = logsDiv.scrollHeight;
          } catch (e) {
            console.error('Poll error:', e);
          }
          
          if (polling) {
            setTimeout(pollLogs, 500);
          }
        }
        
        async function ask() {
          const select = document.getElementById('questionSelect');
          const btn = document.getElementById('askBtn');
          const output = document.getElementById('output');
          
          const index = select.value;
          if (!index) {
            alert('Please select a question');
            return;
          }
          
          const q = questions[index];
          const questionText = 'Question: ' + q.question + '\\n\\nOptions:\\n' +
            Object.entries(q.options).map(([k, v]) => k + ') ' + v).join('\\n');
          
          btn.disabled = true;
          output.textContent = 'Processing...';
          polling = true;
          pollLogs();
          
          try {
            const response = await fetch('/ask', {
              method: 'POST',
              headers: { 'Content-Type': 'application/json' },
              body: JSON.stringify({ 
                question: questionText,
                questionId: q.id,
                correctAnswer: q.correct_answer
              })
            });
            
            const result = await response.json();
            
            if (result.success) {
              output.textContent = result.response;
              document.getElementById('logFile').textContent = result.logFile || 'None';
            } else {
              output.textContent = 'Error: ' + result.error;
            }
          } catch (error) {
            output.textContent = 'Error: ' + error.message;
          } finally {
            btn.disabled = false;
            setTimeout(() => { polling = false; }, 2000);
          }
        }
        
        function clearLogs() {
          document.getElementById('logs').innerHTML = '';
        }
        
        // Start polling when page loads
        polling = true;
        pollLogs();
      </script>
    </body>
    </html>
  `);
});

// API endpoint to get logs
app.get('/logs', (req, res) => {
  res.json({ logs: logs.slice(-100) }); // Last 100 logs
});

// API endpoint to ask questions
app.post('/ask', async (req, res) => {
  const { question, questionId, correctAnswer } = req.body;
  
  // Reset logs
  logs = [];
  
  // Create new log file
  const timestamp = new Date().toISOString().replace(/[:.]/g, '-');
  currentLogFile = path.join(logsDir, `session-${questionId}-${timestamp}.jsonl`);
  
  log('info', `Starting new session for question ${questionId}`);
  log('info', `Log file: ${currentLogFile}`);
  log('info', 'Question received, length: ' + question.length);
  
  // Intercept logger to capture workflow logs
  const logger = require('./src/utils/logger');
  const originalInfo = logger.info.bind(logger);
  
  logger.info = function(message, metadata) {
    originalInfo(message, metadata);
    if (metadata) {
      log('workflow', message, metadata);
    }
  };
  
  // Intercept console
  const originalLog = console.log;
  const originalError = console.error;
  
  console.log = (...args) => {
    originalLog(...args);
    log('console', args.join(' '));
  };
  
  console.error = (...args) => {
    originalError(...args);
    log('error', args.join(' '));
  };
  
  try {
    log('info', 'Calling workflow engine...');
    const startTime = Date.now();
    
    const response = await workflowEngine.answer(question, { timeout: 300000 });
    
    const duration = Date.now() - startTime;
    log('success', `Response received in ${duration}ms`);
    
    // Check if correct
    if (correctAnswer) {
      const responseUpper = response.content.toUpperCase();
      const isCorrect = responseUpper.includes(correctAnswer + ')') || 
                       responseUpper.includes('ANSWER IS ' + correctAnswer) ||
                       responseUpper.includes('ANSWER: ' + correctAnswer);
      log('info', `Correct answer: ${correctAnswer}, Model correct: ${isCorrect}`);
    }
    
    // Restore console
    console.log = originalLog;
    console.error = originalError;
    logger.info = originalInfo;
    
    res.json({
      success: true,
      response: response.content,
      duration: duration,
      logFile: path.basename(currentLogFile)
    });
    
  } catch (error) {
    log('error', 'Error: ' + error.message);
    
    // Restore console
    console.log = originalLog;
    console.error = originalError;
    logger.info = originalInfo;
    
    res.json({
      success: false,
      error: error.message
    });
  }
});

// Serve GPQA questions
app.get('/gpqa-questions.json', (req, res) => {
  res.sendFile(path.join(__dirname, 'gpqa-10-questions-clean.json'));
});

const PORT = 3007;
app.listen(PORT, () => {
  console.log(`Simple UI with logs running on http://localhost:${PORT}`);
  console.log(`Logs will be saved to: ${logsDir}`);
});