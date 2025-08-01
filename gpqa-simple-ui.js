const express = require('express');
const fs = require('fs');
const path = require('path');
const workflowEngine = require('./src/workflow-engine');

const app = express();
app.use(express.json());

// Create logs directory
const logsDir = path.join(__dirname, 'gpqa-logs');
if (!fs.existsSync(logsDir)) {
  fs.mkdirSync(logsDir, { recursive: true });
}

let currentLogFile = null;
let logs = [];

// Intercept all console outputs
const originalLog = console.log;
const originalError = console.error;
const originalWarn = console.warn;
const originalInfo = console.info;

function captureLog(type, args) {
  const message = args.map(arg => 
    typeof arg === 'object' ? JSON.stringify(arg, null, 2) : String(arg)
  ).join(' ');
  
  const entry = {
    timestamp: new Date().toISOString(),
    type,
    message
  };
  
  logs.push(entry);
  
  if (currentLogFile) {
    fs.appendFileSync(currentLogFile, JSON.stringify(entry) + '\n');
  }
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

console.info = (...args) => {
  originalInfo(...args);
  captureLog('info', args);
};

// Load GPQA questions
const gpqaQuestions = JSON.parse(fs.readFileSync(path.join(__dirname, 'gpqa-10-questions-clean.json'), 'utf8'));

// Simple HTML UI
app.get('/', (req, res) => {
  res.send(`
    <!DOCTYPE html>
    <html>
    <head>
      <title>GPQA Test UI</title>
      <style>
        body { 
          font-family: monospace; 
          padding: 20px; 
          background: #000; 
          color: #0f0; 
          margin: 0;
        }
        .container { 
          display: flex; 
          gap: 20px; 
          height: calc(100vh - 40px);
        }
        .left { 
          flex: 1; 
          display: flex;
          flex-direction: column;
        }
        .right { 
          flex: 1; 
          display: flex;
          flex-direction: column;
        }
        #response {
          background: #111;
          border: 1px solid #0f0;
          padding: 10px;
          margin: 10px 0;
          max-height: 300px;
          overflow-y: auto;
          white-space: pre-wrap;
        }
        .questions {
          background: #111;
          border: 1px solid #0f0;
          padding: 10px;
          margin-bottom: 10px;
          max-height: 200px;
          overflow-y: auto;
        }
        .question-btn {
          display: block;
          width: 100%;
          padding: 5px;
          margin: 2px 0;
          background: #222;
          color: #0f0;
          border: 1px solid #0f0;
          cursor: pointer;
          text-align: left;
          font-family: monospace;
        }
        .question-btn:hover {
          background: #333;
        }
        .question-display {
          background: #111;
          border: 1px solid #0f0;
          padding: 10px;
          flex: 1;
          overflow-y: auto;
          white-space: pre-wrap;
        }
        button#askBtn {
          padding: 10px;
          background: #0f0;
          color: #000;
          border: none;
          font-weight: bold;
          cursor: pointer;
          margin: 10px 0;
        }
        button#askBtn:disabled {
          background: #555;
          cursor: not-allowed;
        }
        #logs {
          background: #111;
          border: 1px solid #0f0;
          padding: 10px;
          flex: 1;
          overflow-y: auto;
          font-size: 12px;
        }
        .log-entry {
          margin: 2px 0;
          white-space: pre-wrap;
          word-wrap: break-word;
        }
        .log-error { color: #f00; }
        .log-warn { color: #ff0; }
        h3 { color: #0f0; margin: 10px 0; }
        #status {
          padding: 10px;
          background: #111;
          border: 1px solid #0f0;
          margin: 10px 0;
        }
      </style>
    </head>
    <body>
      <div class="container">
        <div class="left">
          <h3>GPQA Questions</h3>
          <div class="questions">
            ${gpqaQuestions.questions.map((q, i) => 
              `<button class="question-btn" onclick="selectQuestion(${i})">${q.id}: ${q.question.substring(0, 50)}...</button>`
            ).join('')}
          </div>
          
          <h3>Selected Question</h3>
          <div class="question-display" id="questionDisplay">No question selected</div>
          
          <button id="askBtn" onclick="askQuestion()" disabled>Ask Question</button>
          
          <div id="status">Status: Ready</div>
          
          <h3>Full Response</h3>
          <div id="response">No response yet</div>
        </div>
        
        <div class="right">
          <h3>Full Logs</h3>
          <button onclick="clearLogs()" style="padding: 5px; margin-bottom: 10px;">Clear Logs</button>
          <div id="logs"></div>
        </div>
      </div>
      
      <script>
        let selectedQuestion = null;
        const questions = ${JSON.stringify(gpqaQuestions.questions)};
        let pollingInterval = null;
        
        function selectQuestion(index) {
          selectedQuestion = questions[index];
          const display = document.getElementById('questionDisplay');
          display.textContent = 'Question: ' + selectedQuestion.question + '\\n\\n';
          display.textContent += 'Options:\\n';
          Object.entries(selectedQuestion.options).forEach(([k, v]) => {
            display.textContent += k + ') ' + v + '\\n';
          });
          display.textContent += '\\nCorrect Answer: ' + selectedQuestion.correct_answer;
          
          document.getElementById('askBtn').disabled = false;
        }
        
        async function askQuestion() {
          if (!selectedQuestion) return;
          
          const btn = document.getElementById('askBtn');
          const status = document.getElementById('status');
          
          btn.disabled = true;
          status.textContent = 'Status: Processing...';
          
          // Start polling logs
          startPolling();
          
          const questionText = 'Question: ' + selectedQuestion.question + '\\n\\nOptions:\\n' +
            Object.entries(selectedQuestion.options).map(([k, v]) => k + ') ' + v).join('\\n');
          
          try {
            const response = await fetch('/ask', {
              method: 'POST',
              headers: { 'Content-Type': 'application/json' },
              body: JSON.stringify({ 
                question: questionText,
                questionId: selectedQuestion.id
              })
            });
            
            const result = await response.json();
            
            if (result.success) {
              // Show full response
              document.getElementById('response').textContent = result.response;
              
              // Update status
              status.textContent = 'Status: Complete in ' + result.duration + 'ms';
              
              // Check if correct
              const responseUpper = result.response.toUpperCase();
              const correct = selectedQuestion.correct_answer;
              const isCorrect = responseUpper.includes(correct + ')') || 
                               responseUpper.includes('ANSWER IS ' + correct) ||
                               responseUpper.includes('ANSWER: ' + correct);
              status.textContent += ' | ' + (isCorrect ? 'CORRECT' : 'INCORRECT');
            } else {
              status.textContent = 'Status: Error - ' + result.error;
              document.getElementById('response').textContent = 'Error: ' + result.error;
            }
          } catch (error) {
            status.textContent = 'Status: Error - ' + error.message;
          } finally {
            btn.disabled = false;
            setTimeout(() => stopPolling(), 2000);
          }
        }
        
        function startPolling() {
          pollingInterval = setInterval(async () => {
            try {
              const response = await fetch('/logs');
              const data = await response.json();
              updateLogs(data.logs);
            } catch (e) {
              console.error('Poll error:', e);
            }
          }, 500);
        }
        
        function stopPolling() {
          if (pollingInterval) {
            clearInterval(pollingInterval);
            pollingInterval = null;
          }
        }
        
        function updateLogs(logEntries) {
          const logsDiv = document.getElementById('logs');
          logsDiv.innerHTML = '';
          
          logEntries.forEach(log => {
            const div = document.createElement('div');
            div.className = 'log-entry log-' + log.type;
            const time = new Date(log.timestamp).toLocaleTimeString();
            div.textContent = time + ' [' + log.type + '] ' + log.message;
            logsDiv.appendChild(div);
          });
          
          logsDiv.scrollTop = logsDiv.scrollHeight;
        }
        
        function clearLogs() {
          document.getElementById('logs').innerHTML = '';
        }
        
        // Initial log poll
        startPolling();
      </script>
    </body>
    </html>
  `);
});

// API endpoint to get logs
app.get('/logs', (req, res) => {
  res.json({ logs: logs }); // Return ALL logs
});

// API endpoint to ask questions
app.post('/ask', async (req, res) => {
  const { question, questionId } = req.body;
  
  // Reset logs for new question
  logs = [];
  
  // Create new log file
  const timestamp = new Date().toISOString().replace(/[:.]/g, '-');
  currentLogFile = path.join(logsDir, `gpqa-${questionId}-${timestamp}.jsonl`);
  
  console.log(`Starting question ${questionId}`);
  console.log(`Question length: ${question.length} characters`);
  console.log(`Log file: ${currentLogFile}`);
  
  // Also intercept the logger module
  try {
    const logger = require('./src/utils/logger');
    const originalInfo = logger.info.bind(logger);
    const originalError = logger.error.bind(logger);
    const originalWarn = logger.warn.bind(logger);
    const originalDebug = logger.debug ? logger.debug.bind(logger) : null;
    
    logger.info = function(...args) {
      originalInfo(...args);
      captureLog('logger-info', args);
    };
    
    logger.error = function(...args) {
      originalError(...args);
      captureLog('logger-error', args);
    };
    
    logger.warn = function(...args) {
      originalWarn(...args);
      captureLog('logger-warn', args);
    };
    
    if (originalDebug) {
      logger.debug = function(...args) {
        originalDebug(...args);
        captureLog('logger-debug', args);
      };
    }
  } catch (e) {
    console.log('Could not intercept logger:', e.message);
  }
  
  // Intercept agent logs
  try {
    const mathAgent = require('./src/agents/math-agent');
    const codeAgent = require('./src/agents/code-agent');
    
    // Log math agent activity
    const originalMathExecute = mathAgent.execute;
    mathAgent.execute = async function(query, options) {
      console.log('[MATH AGENT] Starting execution for query:', query);
      const result = await originalMathExecute.call(this, query, options);
      console.log('[MATH AGENT] Result:', JSON.stringify(result, null, 2));
      return result;
    };
    
    // Log code agent activity
    const originalCodeExecute = codeAgent.execute;
    codeAgent.execute = async function(query, options) {
      console.log('[CODE AGENT] Starting execution for query:', query);
      const result = await originalCodeExecute.call(this, query, options);
      console.log('[CODE AGENT] Result:', JSON.stringify(result, null, 2));
      return result;
    };
  } catch (e) {
    console.log('Could not intercept agents:', e.message);
  }
  
  try {
    console.log('Calling workflow engine...');
    const startTime = Date.now();
    
    const response = await workflowEngine.answer(question, { timeout: 300000 });
    
    const duration = Date.now() - startTime;
    console.log(`Response received in ${duration}ms`);
    console.log(`Full response length: ${response.content.length} characters`);
    console.log('FULL RESPONSE:', response.content);
    
    res.json({
      success: true,
      response: response.content,
      duration: duration,
      logFile: path.basename(currentLogFile)
    });
    
  } catch (error) {
    console.error('Error processing question:', error.message);
    console.error('Stack trace:', error.stack);
    
    res.json({
      success: false,
      error: error.message
    });
  }
});

const PORT = 3008;
app.listen(PORT, () => {
  console.log(`GPQA Simple UI running on http://localhost:${PORT}`);
  console.log(`Logs will be saved to: ${logsDir}`);
  console.log(`Loaded ${gpqaQuestions.questions.length} GPQA questions`);
});