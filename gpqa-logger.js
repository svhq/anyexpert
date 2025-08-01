const express = require('express');
const fs = require('fs');
const path = require('path');

// First setup logging before importing anything else
const logsDir = path.join(__dirname, 'gpqa-logs-v2');
if (!fs.existsSync(logsDir)) {
  fs.mkdirSync(logsDir, { recursive: true });
}

let logStream = null;
const allLogs = [];

function writeLog(entry) {
  allLogs.push(entry);
  if (logStream) {
    logStream.write(JSON.stringify(entry) + '\n');
  }
  // Also write to console
  const time = new Date(entry.timestamp).toLocaleTimeString();
  console.error(`[${time}] [${entry.type}] ${entry.message}`);
}

// Override console methods BEFORE requiring any modules
const originalConsole = {
  log: console.log,
  error: console.error,
  warn: console.warn,
  info: console.info
};

console.log = (...args) => {
  originalConsole.log(...args);
  writeLog({
    timestamp: new Date().toISOString(),
    type: 'console.log',
    message: args.map(a => typeof a === 'object' ? JSON.stringify(a, null, 2) : String(a)).join(' ')
  });
};

console.error = (...args) => {
  originalConsole.error(...args);
  writeLog({
    timestamp: new Date().toISOString(),
    type: 'console.error',
    message: args.map(a => typeof a === 'object' ? JSON.stringify(a, null, 2) : String(a)).join(' ')
  });
};

console.warn = (...args) => {
  originalConsole.warn(...args);
  writeLog({
    timestamp: new Date().toISOString(),
    type: 'console.warn',
    message: args.map(a => typeof a === 'object' ? JSON.stringify(a, null, 2) : String(a)).join(' ')
  });
};

console.info = (...args) => {
  originalConsole.info(...args);
  writeLog({
    timestamp: new Date().toISOString(),
    type: 'console.info',
    message: args.map(a => typeof a === 'object' ? JSON.stringify(a, null, 2) : String(a)).join(' ')
  });
};

// Now require modules
const workflowEngine = require('./src/workflow-engine');

// Override logger
try {
  const logger = require('./src/utils/logger');
  const methods = ['info', 'error', 'warn', 'debug'];
  methods.forEach(method => {
    if (logger[method]) {
      const original = logger[method].bind(logger);
      logger[method] = (...args) => {
        original(...args);
        writeLog({
          timestamp: new Date().toISOString(),
          type: `logger.${method}`,
          message: args.map(a => typeof a === 'object' ? JSON.stringify(a, null, 2) : String(a)).join(' ')
        });
      };
    }
  });
} catch (e) {
  writeLog({
    timestamp: new Date().toISOString(),
    type: 'error',
    message: 'Failed to override logger: ' + e.message
  });
}

// Load GPQA questions
const gpqaQuestions = JSON.parse(fs.readFileSync(path.join(__dirname, 'gpqa-10-questions-clean.json'), 'utf8'));

// Express app
const app = express();
app.use(express.json());

app.get('/', (req, res) => {
  res.send(`
    <!DOCTYPE html>
    <html>
    <head>
      <title>GPQA Logger Test</title>
      <style>
        body { font-family: monospace; background: #1a1a1a; color: #e0e0e0; padding: 20px; }
        .container { display: flex; gap: 20px; height: calc(100vh - 40px); }
        .left, .right { flex: 1; display: flex; flex-direction: column; }
        .box { background: #0a0a0a; border: 1px solid #333; padding: 10px; margin: 10px 0; overflow-y: auto; }
        .questions { max-height: 200px; }
        .question-btn { display: block; width: 100%; padding: 5px; margin: 2px 0; background: #222; color: #e0e0e0; border: 1px solid #444; cursor: pointer; text-align: left; }
        .question-btn:hover { background: #333; }
        #questionDisplay, #response { flex: 1; white-space: pre-wrap; font-size: 12px; }
        #logs { flex: 1; font-size: 11px; line-height: 1.3; }
        button { padding: 10px; background: #00ff88; color: #000; border: none; font-weight: bold; cursor: pointer; }
        button:disabled { background: #555; color: #999; }
        .log-entry { margin: 1px 0; }
        .log-console-log { color: #e0e0e0; }
        .log-console-error { color: #ff6b6b; }
        .log-logger-info { color: #4ecdc4; }
        .log-logger-error { color: #ff6b6b; }
        h3 { margin: 5px 0; color: #00ff88; }
      </style>
    </head>
    <body>
      <div class="container">
        <div class="left">
          <h3>GPQA Questions</h3>
          <div class="box questions">
            ${gpqaQuestions.questions.map((q, i) => 
              `<button class="question-btn" onclick="selectQuestion(${i})">${q.id}</button>`
            ).join('')}
          </div>
          
          <h3>Question</h3>
          <div class="box" id="questionDisplay">Select a question</div>
          
          <button id="askBtn" onclick="askQuestion()" disabled>Ask Question</button>
          
          <h3>Response</h3>
          <div class="box" id="response">No response yet</div>
        </div>
        
        <div class="right">
          <h3>Live Logs (${allLogs.length} entries)</h3>
          <div class="box" id="logs"></div>
        </div>
      </div>
      
      <script>
        let selectedQuestion = null;
        const questions = ${JSON.stringify(gpqaQuestions.questions)};
        
        function selectQuestion(index) {
          selectedQuestion = questions[index];
          const display = document.getElementById('questionDisplay');
          display.textContent = 'Question: ' + selectedQuestion.question + '\\n\\nOptions:\\n';
          Object.entries(selectedQuestion.options).forEach(([k, v]) => {
            display.textContent += k + ') ' + v + '\\n';
          });
          document.getElementById('askBtn').disabled = false;
        }
        
        async function askQuestion() {
          if (!selectedQuestion) return;
          
          document.getElementById('askBtn').disabled = true;
          document.getElementById('response').textContent = 'Processing...';
          
          const questionText = 'Question: ' + selectedQuestion.question + '\\n\\nOptions:\\n' +
            Object.entries(selectedQuestion.options).map(([k, v]) => k + ') ' + v).join('\\n');
          
          try {
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
            document.getElementById('response').textContent = result.success ? result.response : 'Error: ' + result.error;
            
          } catch (error) {
            document.getElementById('response').textContent = 'Error: ' + error.message;
          } finally {
            document.getElementById('askBtn').disabled = false;
          }
        }
        
        // Poll logs
        setInterval(async () => {
          try {
            const response = await fetch('/logs');
            const data = await response.json();
            const logsDiv = document.getElementById('logs');
            
            // Update log count
            document.querySelector('h3').textContent = 'Live Logs (' + data.total + ' entries)';
            
            logsDiv.innerHTML = data.logs.map(log => {
              const time = new Date(log.timestamp).toLocaleTimeString();
              const className = 'log-entry log-' + log.type.replace('.', '-');
              return '<div class="' + className + '">' + time + ' [' + log.type + '] ' + 
                     log.message.replace(/</g, '&lt;').replace(/>/g, '&gt;') + '</div>';
            }).join('');
            
            logsDiv.scrollTop = logsDiv.scrollHeight;
          } catch (e) {}
        }, 500);
      </script>
    </body>
    </html>
  `);
});

app.get('/logs', (req, res) => {
  res.json({ 
    logs: allLogs.slice(-500), // Last 500 for display
    total: allLogs.length 
  });
});

app.post('/ask', async (req, res) => {
  const { question, questionId, correctAnswer } = req.body;
  
  // Create new log file
  const timestamp = new Date().toISOString().replace(/[:.]/g, '-');
  const logFile = path.join(logsDir, `gpqa-${questionId}-${timestamp}.jsonl`);
  logStream = fs.createWriteStream(logFile);
  
  writeLog({
    timestamp: new Date().toISOString(),
    type: 'session',
    message: `Starting question ${questionId}, log file: ${logFile}`
  });
  
  try {
    const startTime = Date.now();
    const response = await workflowEngine.answer(question, { timeout: 300000 });
    const duration = Date.now() - startTime;
    
    writeLog({
      timestamp: new Date().toISOString(),
      type: 'result',
      message: `Complete in ${duration}ms. Response length: ${response.content.length} chars`
    });
    
    // Check correctness
    if (correctAnswer) {
      const responseUpper = response.content.toUpperCase();
      const isCorrect = responseUpper.includes(correctAnswer + ')') || 
                       responseUpper.includes('ANSWER IS ' + correctAnswer) ||
                       responseUpper.includes('ANSWER: ' + correctAnswer);
      
      writeLog({
        timestamp: new Date().toISOString(),
        type: 'validation',
        message: `Correct answer: ${correctAnswer}, Model answered correctly: ${isCorrect}`
      });
    }
    
    res.json({
      success: true,
      response: response.content,
      duration: duration
    });
    
  } catch (error) {
    writeLog({
      timestamp: new Date().toISOString(),
      type: 'error',
      message: error.stack || error.message
    });
    
    res.json({
      success: false,
      error: error.message
    });
  } finally {
    if (logStream) {
      logStream.end();
      logStream = null;
    }
  }
});

const PORT = 3009;
app.listen(PORT, () => {
  writeLog({
    timestamp: new Date().toISOString(),
    type: 'server',
    message: `GPQA Logger running on http://localhost:${PORT}`
  });
});