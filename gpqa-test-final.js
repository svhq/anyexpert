const express = require('express');
const fs = require('fs');
const path = require('path');
const workflowEngine = require('./src/workflow-engine');

const app = express();
app.use(express.json({ limit: '50mb' }));

// Load GPQA questions
const gpqaQuestions = JSON.parse(fs.readFileSync(path.join(__dirname, 'gpqa-10-questions-clean.json'), 'utf8'));

// Simple logging to file
const logsDir = path.join(__dirname, 'gpqa-test-logs');
if (!fs.existsSync(logsDir)) {
  fs.mkdirSync(logsDir, { recursive: true });
}

let currentSession = null;

function logToFile(type, message) {
  if (currentSession) {
    const entry = {
      timestamp: new Date().toISOString(),
      type,
      message
    };
    fs.appendFileSync(currentSession.logFile, JSON.stringify(entry) + '\n');
  }
}

// Intercept console
const originalLog = console.log;
const originalError = console.error;

console.log = (...args) => {
  originalLog(...args);
  logToFile('log', args.join(' '));
};

console.error = (...args) => {
  originalError(...args);
  logToFile('error', args.join(' '));
};

// HTML UI
app.get('/', (req, res) => {
  res.send(`
    <!DOCTYPE html>
    <html>
    <head>
      <title>GPQA Test - Final Working Version</title>
      <style>
        * { box-sizing: border-box; }
        body { 
          font-family: 'Consolas', 'Monaco', monospace; 
          background: #0a0a0a; 
          color: #e0e0e0; 
          margin: 0;
          padding: 20px;
        }
        .container { max-width: 1200px; margin: 0 auto; }
        .header { text-align: center; margin-bottom: 30px; }
        h1 { color: #00ff88; margin: 0; }
        .subtitle { color: #888; margin-top: 10px; }
        
        .questions-grid {
          display: grid;
          grid-template-columns: repeat(auto-fill, minmax(300px, 1fr));
          gap: 10px;
          margin-bottom: 20px;
        }
        
        .question-card {
          background: #1a1a1a;
          border: 1px solid #333;
          padding: 15px;
          cursor: pointer;
          transition: all 0.2s;
        }
        
        .question-card:hover {
          border-color: #00ff88;
          background: #222;
        }
        
        .question-card.selected {
          border-color: #00ff88;
          background: #1a3a1a;
        }
        
        .question-id { color: #00ff88; font-weight: bold; }
        .question-preview { 
          color: #ccc; 
          font-size: 12px; 
          margin-top: 5px;
          overflow: hidden;
          text-overflow: ellipsis;
          white-space: nowrap;
        }
        
        .control-panel {
          background: #1a1a1a;
          border: 1px solid #333;
          padding: 20px;
          margin-bottom: 20px;
        }
        
        .question-display {
          background: #0a0a0a;
          border: 1px solid #444;
          padding: 15px;
          margin: 15px 0;
          white-space: pre-wrap;
          font-size: 13px;
          max-height: 300px;
          overflow-y: auto;
        }
        
        button {
          background: #00ff88;
          color: #000;
          border: none;
          padding: 12px 30px;
          font-weight: bold;
          cursor: pointer;
          font-size: 16px;
          transition: all 0.2s;
        }
        
        button:hover { background: #00cc66; }
        button:disabled { 
          background: #444; 
          color: #888; 
          cursor: not-allowed; 
        }
        
        .status {
          margin: 15px 0;
          padding: 10px;
          background: #222;
          border-left: 3px solid #00ff88;
        }
        
        .response-section {
          background: #1a1a1a;
          border: 1px solid #333;
          padding: 20px;
          margin-bottom: 20px;
        }
        
        .response-content {
          background: #0a0a0a;
          border: 1px solid #444;
          padding: 15px;
          white-space: pre-wrap;
          font-size: 13px;
          max-height: 500px;
          overflow-y: auto;
        }
        
        .logs-section {
          background: #1a1a1a;
          border: 1px solid #333;
          padding: 20px;
        }
        
        .logs-content {
          background: #0a0a0a;
          border: 1px solid #444;
          padding: 10px;
          font-size: 11px;
          height: 300px;
          overflow-y: auto;
          font-family: monospace;
        }
        
        .log-entry { margin: 2px 0; }
        .log-entry.error { color: #ff6b6b; }
        .log-entry.log { color: #e0e0e0; }
        
        .correct { color: #00ff88; }
        .incorrect { color: #ff6b6b; }
        
        .warning {
          background: #332200;
          border: 1px solid #ffaa00;
          padding: 10px;
          margin: 10px 0;
          color: #ffaa00;
        }
      </style>
    </head>
    <body>
      <div class="container">
        <div class="header">
          <h1>GPQA Test Interface</h1>
          <div class="subtitle">Working Version - E2B Service Connected</div>
        </div>
        
        <div class="warning">
          ⚠️ Note: numpy is not available in E2B environment. Questions requiring numpy may fail code execution.
        </div>
        
        <h2>Select a GPQA Question:</h2>
        <div class="questions-grid">
          ${gpqaQuestions.questions.map((q, i) => `
            <div class="question-card" onclick="selectQuestion(${i})" id="q${i}">
              <div class="question-id">${q.id}</div>
              <div class="question-preview">${q.question}</div>
            </div>
          `).join('')}
        </div>
        
        <div class="control-panel">
          <h3>Selected Question:</h3>
          <div class="question-display" id="questionDisplay">No question selected</div>
          <button id="askBtn" onclick="askQuestion()" disabled>Ask Question</button>
          <div class="status" id="status">Status: Select a question to begin</div>
        </div>
        
        <div class="response-section">
          <h3>Response:</h3>
          <div class="response-content" id="response">No response yet</div>
        </div>
        
        <div class="logs-section">
          <h3>Session Logs: <button onclick="clearLogs()" style="float: right; padding: 5px 10px; font-size: 12px;">Clear</button></h3>
          <div class="logs-content" id="logs"></div>
        </div>
      </div>
      
      <script>
        let selectedQuestion = null;
        let selectedIndex = -1;
        const questions = ${JSON.stringify(gpqaQuestions.questions)};
        let logs = [];
        
        function selectQuestion(index) {
          // Update UI
          document.querySelectorAll('.question-card').forEach(card => card.classList.remove('selected'));
          document.getElementById('q' + index).classList.add('selected');
          
          selectedIndex = index;
          selectedQuestion = questions[index];
          
          // Display full question
          const display = document.getElementById('questionDisplay');
          let text = 'Question: ' + selectedQuestion.question + '\\n\\nOptions:\\n';
          Object.entries(selectedQuestion.options).forEach(([k, v]) => {
            text += k + ') ' + v + '\\n';
          });
          text += '\\nCorrect Answer: ' + selectedQuestion.correct_answer;
          
          display.textContent = text;
          document.getElementById('askBtn').disabled = false;
          document.getElementById('status').textContent = 'Status: Ready to ask question';
        }
        
        async function askQuestion() {
          if (!selectedQuestion) return;
          
          const btn = document.getElementById('askBtn');
          const status = document.getElementById('status');
          const responseDiv = document.getElementById('response');
          
          btn.disabled = true;
          status.textContent = 'Status: Processing... (this may take 30-60 seconds)';
          responseDiv.textContent = 'Processing...';
          
          // Clear logs
          logs = [];
          updateLogs();
          
          // Start polling logs
          const pollInterval = setInterval(async () => {
            try {
              const res = await fetch('/session-logs');
              if (res.ok) {
                const data = await res.json();
                logs = data.logs;
                updateLogs();
              }
            } catch (e) {}
          }, 1000);
          
          const questionText = 'Question: ' + selectedQuestion.question + '\\n\\nOptions:\\n' +
            Object.entries(selectedQuestion.options).map(([k, v]) => k + ') ' + v).join('\\n');
          
          try {
            const startTime = Date.now();
            
            const response = await fetch('/ask', {
              method: 'POST',
              headers: { 'Content-Type': 'application/json' },
              body: JSON.stringify({ 
                question: questionText,
                questionId: selectedQuestion.id,
                correctAnswer: selectedQuestion.correct_answer
              })
            });
            
            const result = await response.json();
            const duration = Date.now() - startTime;
            
            clearInterval(pollInterval);
            
            if (result.success) {
              responseDiv.textContent = result.response;
              
              // Check if correct
              const responseUpper = result.response.toUpperCase();
              const correct = selectedQuestion.correct_answer;
              const isCorrect = responseUpper.includes(correct + ')') || 
                               responseUpper.includes('ANSWER IS ' + correct) ||
                               responseUpper.includes('ANSWER: ' + correct) ||
                               responseUpper.includes('**' + correct + ')**');
              
              status.innerHTML = 'Status: Complete in ' + (duration/1000).toFixed(1) + 's | Answer: ' + 
                                '<span class="' + (isCorrect ? 'correct' : 'incorrect') + '">' +
                                (isCorrect ? 'CORRECT ✓' : 'INCORRECT ✗') + '</span>';
              
              addLog('info', 'Test complete. Correct answer: ' + correct + ', Model correct: ' + isCorrect);
            } else {
              responseDiv.textContent = 'Error: ' + result.error;
              status.textContent = 'Status: Error occurred';
              addLog('error', 'Error: ' + result.error);
            }
          } catch (error) {
            clearInterval(pollInterval);
            responseDiv.textContent = 'Error: ' + error.message;
            status.textContent = 'Status: Request failed';
            addLog('error', 'Request error: ' + error.message);
          } finally {
            btn.disabled = false;
            
            // Final log fetch
            setTimeout(async () => {
              try {
                const res = await fetch('/session-logs');
                if (res.ok) {
                  const data = await res.json();
                  logs = data.logs;
                  updateLogs();
                }
              } catch (e) {}
            }, 1000);
          }
        }
        
        function addLog(type, message) {
          logs.push({
            timestamp: new Date().toISOString(),
            type,
            message
          });
          updateLogs();
        }
        
        function updateLogs() {
          const logsDiv = document.getElementById('logs');
          logsDiv.innerHTML = logs.map(log => {
            const time = new Date(log.timestamp).toLocaleTimeString();
            return '<div class="log-entry ' + log.type + '">' + 
                   time + ' [' + log.type + '] ' + 
                   log.message.replace(/</g, '&lt;').replace(/>/g, '&gt;') + 
                   '</div>';
          }).join('');
          logsDiv.scrollTop = logsDiv.scrollHeight;
        }
        
        function clearLogs() {
          logs = [];
          updateLogs();
        }
      </script>
    </body>
    </html>
  `);
});

// Session logs endpoint
app.get('/session-logs', (req, res) => {
  if (!currentSession) {
    return res.json({ logs: [] });
  }
  
  try {
    const logs = fs.readFileSync(currentSession.logFile, 'utf8')
      .split('\n')
      .filter(line => line.trim())
      .map(line => JSON.parse(line));
    res.json({ logs });
  } catch (e) {
    res.json({ logs: [] });
  }
});

// Ask endpoint
app.post('/ask', async (req, res) => {
  const { question, questionId, correctAnswer } = req.body;
  
  // Create new session
  const timestamp = new Date().toISOString().replace(/[:.]/g, '-');
  currentSession = {
    logFile: path.join(logsDir, `${questionId}-${timestamp}.jsonl`)
  };
  
  logToFile('info', `Starting test for question ${questionId}`);
  logToFile('info', `Correct answer: ${correctAnswer}`);
  
  try {
    console.log(`\nProcessing ${questionId}...`);
    const startTime = Date.now();
    
    const response = await workflowEngine.answer(question, { timeout: 120000 });
    
    const duration = Date.now() - startTime;
    logToFile('info', `Response received in ${duration}ms`);
    
    res.json({
      success: true,
      response: response.content,
      duration
    });
    
  } catch (error) {
    console.error('Error:', error);
    logToFile('error', error.stack || error.message);
    
    res.json({
      success: false,
      error: error.message
    });
  } finally {
    currentSession = null;
  }
});

const PORT = 3013;
app.listen(PORT, () => {
  console.log(`\nGPQA Test Final running on http://localhost:${PORT}`);
  console.log(`E2B service status: Connected at ${require('./config').codeExecution.serviceUrl}`);
  console.log(`Loaded ${gpqaQuestions.questions.length} GPQA questions`);
  console.log('\nNOTE: numpy is not available in E2B environment');
  console.log('Questions requiring numpy may fail code execution\n');
});