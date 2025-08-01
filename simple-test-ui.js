const express = require('express');
const fs = require('fs');
const path = require('path');
const workflowEngine = require('./src/workflow-engine');

const app = express();
app.use(express.json());

// Load GPQA questions
const gpqaQuestions = JSON.parse(fs.readFileSync(path.join(__dirname, 'gpqa-10-questions-clean.json'), 'utf8'));

// Store logs in memory
let logs = [];
let isProcessing = false;

// Capture all console output
const originalLog = console.log;
const originalError = console.error;
const originalWarn = console.warn;

function captureLog(type, args) {
  const message = args.map(arg => 
    typeof arg === 'object' ? JSON.stringify(arg, null, 2) : String(arg)
  ).join(' ');
  
  logs.push({
    timestamp: new Date().toISOString(),
    type,
    message
  });
  
  // Keep only last 1000 logs
  if (logs.length > 1000) logs.shift();
}

console.log = (...args) => {
  originalLog(...args);
  captureLog('log', args);
};

console.error = (...args) => {
  originalError(...args);
  captureLog('error', args);
};

console.warn = (...args) => {
  originalWarn(...args);
  captureLog('warn', args);
};

// Intercept logger
setTimeout(() => {
  try {
    const logger = require('./src/utils/logger');
    ['info', 'error', 'warn', 'debug'].forEach(method => {
      if (logger[method]) {
        const original = logger[method].bind(logger);
        logger[method] = (...args) => {
          original(...args);
          captureLog(`logger.${method}`, args);
        };
      }
    });
  } catch (e) {
    console.error('Could not intercept logger:', e.message);
  }
}, 100);

// Simple UI
app.get('/', (req, res) => {
  res.send(`
    <!DOCTYPE html>
    <html>
    <head>
      <title>Ask Any Expert - Test UI</title>
      <style>
        body { font-family: monospace; background: #1a1a1a; color: #e0e0e0; padding: 20px; margin: 0; }
        .container { max-width: 1400px; margin: 0 auto; }
        textarea { width: 100%; background: #0a0a0a; color: #e0e0e0; border: 1px solid #333; padding: 10px; }
        button { padding: 10px 20px; background: #00ff88; color: #000; border: none; font-weight: bold; cursor: pointer; }
        button:disabled { background: #555; color: #999; }
        .output { background: #0a0a0a; border: 1px solid #333; padding: 15px; margin: 10px 0; min-height: 200px; white-space: pre-wrap; word-wrap: break-word; }
        .logs { background: #0a0a0a; border: 1px solid #333; padding: 10px; max-height: 400px; overflow-y: auto; font-size: 12px; }
        .log-entry { margin: 2px 0; padding: 2px; }
        .log-log { color: #e0e0e0; }
        .log-error { color: #ff6b6b; }
        .log-warn { color: #ffeb3b; }
        .log-logger-info { color: #4ecdc4; }
        h2 { color: #00ff88; }
        .status { padding: 10px; background: #0a0a0a; border: 1px solid #333; margin: 10px 0; }
        .examples { margin: 10px 0; }
        .example { display: inline-block; margin: 5px; padding: 5px 10px; background: #222; border: 1px solid #444; cursor: pointer; }
        .example:hover { background: #333; }
      </style>
    </head>
    <body>
      <div class="container">
        <h1>Ask Any Expert - Test UI</h1>
        
        <h2>Ask a Question</h2>
        <textarea id="question" rows="4" placeholder="Enter your question...">Write a Python function to calculate the factorial of 10 and tell me the result.</textarea>
        
        <div class="examples">
          Example questions:
          <span class="example" onclick="setQuestion('What is 2 + 2?')">Simple Math</span>
          <span class="example" onclick="setQuestion('Write a Python function to calculate the factorial of 10 and tell me the result.')">Code Execution</span>
          <span class="example" onclick="setQuestion('What is the capital of France?')">General Knowledge</span>
          <span class="example" onclick="setQuestion('What are the latest developments in quantum computing?')">Search Required</span>
        </div>
        
        <h2>GPQA Questions</h2>
        <select id="gpqaSelect" onchange="selectGPQA()" style="width: 100%; padding: 5px; background: #0a0a0a; color: #e0e0e0; border: 1px solid #333;">
          <option value="">Select a GPQA question...</option>
          ${gpqaQuestions.questions.map((q, i) => 
            `<option value="${i}">${q.id}: ${q.question.substring(0, 80)}...</option>`
          ).join('')}
        </select>
        
        <button id="askBtn" onclick="askQuestion()">Ask Question</button>
        
        <div class="status" id="status">Status: Ready</div>
        
        <h2>Response</h2>
        <div class="output" id="response">No response yet</div>
        
        <h2>Live Logs (<span id="logCount">0</span> entries) <button onclick="clearLogs()">Clear</button></h2>
        <div class="logs" id="logs"></div>
      </div>
      
      <script>
        let pollingInterval = null;
        const gpqaQuestions = ${JSON.stringify(gpqaQuestions.questions)};
        
        function setQuestion(q) {
          document.getElementById('question').value = q;
        }
        
        function selectGPQA() {
          const select = document.getElementById('gpqaSelect');
          if (select.value === '') return;
          
          const q = gpqaQuestions[parseInt(select.value)];
          const questionText = 'Question: ' + q.question + '\\n\\nOptions:\\n' +
            Object.entries(q.options).map(([k, v]) => k + ') ' + v).join('\\n');
          
          document.getElementById('question').value = questionText;
        }
        
        async function askQuestion() {
          const question = document.getElementById('question').value.trim();
          if (!question) {
            alert('Please enter a question');
            return;
          }
          
          const btn = document.getElementById('askBtn');
          const status = document.getElementById('status');
          const response = document.getElementById('response');
          
          btn.disabled = true;
          status.textContent = 'Status: Processing...';
          response.textContent = 'Processing...';
          
          // Start polling logs
          if (pollingInterval) clearInterval(pollingInterval);
          pollingInterval = setInterval(updateLogs, 500);
          
          try {
            const startTime = Date.now();
            const res = await fetch('/ask', {
              method: 'POST',
              headers: { 'Content-Type': 'application/json' },
              body: JSON.stringify({ question })
            });
            
            const result = await res.json();
            const duration = Date.now() - startTime;
            
            if (result.success) {
              response.textContent = result.response;
              status.textContent = 'Status: Complete in ' + duration + 'ms';
            } else {
              response.textContent = 'Error: ' + result.error;
              status.textContent = 'Status: Error';
            }
          } catch (error) {
            response.textContent = 'Error: ' + error.message;
            status.textContent = 'Status: Error';
          } finally {
            btn.disabled = false;
            setTimeout(() => {
              if (pollingInterval) {
                clearInterval(pollingInterval);
                pollingInterval = null;
              }
            }, 2000);
          }
        }
        
        async function updateLogs() {
          try {
            const res = await fetch('/logs');
            const data = await res.json();
            
            document.getElementById('logCount').textContent = data.logs.length;
            
            const logsDiv = document.getElementById('logs');
            logsDiv.innerHTML = data.logs.map(log => {
              const time = new Date(log.timestamp).toLocaleTimeString();
              const className = 'log-entry log-' + log.type.replace('.', '-');
              return '<div class="' + className + '">' + time + ' [' + log.type + '] ' + 
                     log.message.replace(/</g, '&lt;').replace(/>/g, '&gt;') + '</div>';
            }).join('');
            
            logsDiv.scrollTop = logsDiv.scrollHeight;
          } catch (e) {
            console.error('Failed to update logs:', e);
          }
        }
        
        function clearLogs() {
          fetch('/clear-logs', { method: 'POST' });
          document.getElementById('logs').innerHTML = '';
          document.getElementById('logCount').textContent = '0';
        }
        
        // Initial log update
        updateLogs();
        
        // Poll logs every second when not processing
        setInterval(() => {
          if (!pollingInterval) updateLogs();
        }, 1000);
      </script>
    </body>
    </html>
  `);
});

app.get('/logs', (req, res) => {
  res.json({ logs });
});

app.post('/clear-logs', (req, res) => {
  logs = [];
  res.json({ success: true });
});

app.post('/ask', async (req, res) => {
  const { question } = req.body;
  
  if (isProcessing) {
    return res.json({ success: false, error: 'Already processing a question' });
  }
  
  isProcessing = true;
  console.log('Processing question:', question);
  
  try {
    const startTime = Date.now();
    const response = await workflowEngine.answer(question, { timeout: 300000 });
    const duration = Date.now() - startTime;
    
    console.log(`Response received in ${duration}ms`);
    console.log(`Response length: ${response.content.length} characters`);
    
    res.json({
      success: true,
      response: response.content,
      duration
    });
  } catch (error) {
    console.error('Error:', error.message);
    res.json({
      success: false,
      error: error.message
    });
  } finally {
    isProcessing = false;
  }
});

const PORT = 3011;
app.listen(PORT, () => {
  console.log(`Simple Test UI running on http://localhost:${PORT}`);
});