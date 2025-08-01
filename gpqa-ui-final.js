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

// UI with GPQA questions
app.get('/', (req, res) => {
  res.send(`
    <!DOCTYPE html>
    <html>
    <head>
      <title>GPQA Test with Full Logs</title>
      <style>
        body { font-family: monospace; background: #1a1a1a; color: #e0e0e0; padding: 20px; margin: 0; }
        .container { max-width: 1400px; margin: 0 auto; }
        textarea, select { width: 100%; background: #0a0a0a; color: #e0e0e0; border: 1px solid #333; padding: 10px; }
        button { padding: 10px 20px; background: #00ff88; color: #000; border: none; font-weight: bold; cursor: pointer; margin: 10px 0; }
        button:disabled { background: #555; color: #999; }
        .output { background: #0a0a0a; border: 1px solid #333; padding: 15px; margin: 10px 0; min-height: 200px; white-space: pre-wrap; word-wrap: break-word; }
        .logs { background: #0a0a0a; border: 1px solid #333; padding: 10px; max-height: 400px; overflow-y: auto; font-size: 12px; }
        .log-entry { margin: 2px 0; padding: 2px; }
        .log-log { color: #e0e0e0; }
        .log-error { color: #ff6b6b; }
        h2 { color: #00ff88; }
        .status { padding: 10px; background: #0a0a0a; border: 1px solid #333; margin: 10px 0; }
        .gpqa-info { background: #222; padding: 10px; margin: 10px 0; border: 1px solid #444; }
      </style>
    </head>
    <body>
      <div class="container">
        <h1>GPQA Test with Full Logs</h1>
        
        <h2>Select GPQA Question</h2>
        <select id="gpqaSelect" onchange="selectGPQA()" style="padding: 10px; font-size: 14px;">
          <option value="">Select a GPQA question...</option>
          ${gpqaQuestions.questions.map((q, i) => 
            `<option value="${i}">${q.id}: ${q.question.substring(0, 100)}...</option>`
          ).join('')}
        </select>
        
        <div class="gpqa-info" id="gpqaInfo" style="display: none;">
          <strong>Correct Answer:</strong> <span id="correctAnswer"></span>
        </div>
        
        <h2>Question</h2>
        <textarea id="question" rows="8" placeholder="Select a GPQA question or enter your own..."></textarea>
        
        <button id="askBtn" onclick="askQuestion()">Ask Question</button>
        
        <div class="status" id="status">Status: Ready</div>
        
        <h2>Response</h2>
        <div class="output" id="response">No response yet</div>
        
        <h2>Live Logs (<span id="logCount">0</span> entries) <button onclick="clearLogs()">Clear</button></h2>
        <div class="logs" id="logs"></div>
      </div>
      
      <script>
        let pollingInterval = null;
        const gpqaData = ${JSON.stringify(gpqaQuestions.questions)};
        let selectedGPQA = null;
        
        function selectGPQA() {
          const select = document.getElementById('gpqaSelect');
          if (select.value === '') {
            document.getElementById('gpqaInfo').style.display = 'none';
            selectedGPQA = null;
            return;
          }
          
          selectedGPQA = gpqaData[parseInt(select.value)];
          const questionText = 'Question: ' + selectedGPQA.question + '\\n\\nOptions:\\n' +
            Object.entries(selectedGPQA.options).map(([k, v]) => k + ') ' + v).join('\\n');
          
          document.getElementById('question').value = questionText;
          document.getElementById('correctAnswer').textContent = selectedGPQA.correct_answer;
          document.getElementById('gpqaInfo').style.display = 'block';
        }
        
        async function askQuestion() {
          const question = document.getElementById('question').value.trim();
          if (!question) {
            alert('Please enter a question or select a GPQA question');
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
              
              // Check if correct for GPQA
              if (selectedGPQA) {
                const responseUpper = result.response.toUpperCase();
                const correct = selectedGPQA.correct_answer;
                const isCorrect = responseUpper.includes(correct + ')') || 
                                 responseUpper.includes('ANSWER IS ' + correct) ||
                                 responseUpper.includes('ANSWER: ' + correct);
                status.textContent = 'Status: Complete in ' + duration + 'ms | Answer: ' + 
                                   (isCorrect ? 'CORRECT ✓' : 'INCORRECT ✗');
              } else {
                status.textContent = 'Status: Complete in ' + duration + 'ms';
              }
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
  console.log('Processing question:', question.substring(0, 100) + '...');
  
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

const PORT = 3012;
app.listen(PORT, () => {
  console.log(`GPQA UI Final running on http://localhost:${PORT}`);
  console.log(`Loaded ${gpqaQuestions.questions.length} GPQA questions`);
});