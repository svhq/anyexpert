const express = require('express');
const fs = require('fs');
const path = require('path');
const workflowEngine = require('./src/workflow-engine');

const app = express();
app.use(express.json());

// Load GPQA questions
const gpqaQuestions = JSON.parse(fs.readFileSync(path.join(__dirname, 'gpqa-10-questions-clean.json'), 'utf8'));

// Simple UI
app.get('/', (req, res) => {
  res.send(`
    <!DOCTYPE html>
    <html>
    <head>
      <title>GPQA Test</title>
      <style>
        body { font-family: monospace; background: #000; color: #0f0; padding: 20px; }
        .container { max-width: 1200px; margin: 0 auto; }
        .questions { background: #111; border: 1px solid #0f0; padding: 10px; margin-bottom: 20px; }
        .question-btn { display: block; width: 100%; padding: 5px; margin: 2px 0; background: #222; color: #0f0; border: 1px solid #0f0; cursor: pointer; text-align: left; }
        .question-btn:hover { background: #333; }
        .selected { background: #040; }
        pre { background: #111; border: 1px solid #0f0; padding: 10px; white-space: pre-wrap; word-wrap: break-word; }
        button { padding: 10px 20px; background: #0f0; color: #000; border: none; font-weight: bold; cursor: pointer; }
        button:disabled { background: #555; color: #999; }
        #status { margin: 10px 0; padding: 10px; background: #111; border: 1px solid #0f0; }
      </style>
    </head>
    <body>
      <div class="container">
        <h1>GPQA Test</h1>
        
        <div class="questions">
          <h3>Select a question:</h3>
          ${gpqaQuestions.questions.map((q, i) => 
            `<button class="question-btn" onclick="selectQuestion(${i})" id="q${i}">${q.id}: ${q.question.substring(0, 80)}...</button>`
          ).join('')}
        </div>
        
        <div id="status">Status: Ready</div>
        
        <button id="askBtn" onclick="askQuestion()" disabled>Ask Selected Question</button>
        
        <h3>Selected Question:</h3>
        <pre id="questionDisplay">No question selected</pre>
        
        <h3>Response:</h3>
        <pre id="response">No response yet</pre>
        
        <h3>Debug Log:</h3>
        <pre id="debugLog" style="max-height: 300px; overflow-y: auto;">Waiting for question...</pre>
      </div>
      
      <script>
        let selectedQuestion = null;
        let selectedIndex = -1;
        const questions = ${JSON.stringify(gpqaQuestions.questions)};
        
        function log(msg) {
          const debugLog = document.getElementById('debugLog');
          debugLog.textContent += new Date().toLocaleTimeString() + ' - ' + msg + '\\n';
          debugLog.scrollTop = debugLog.scrollHeight;
        }
        
        function selectQuestion(index) {
          // Clear previous selection
          document.querySelectorAll('.question-btn').forEach(btn => btn.classList.remove('selected'));
          document.getElementById('q' + index).classList.add('selected');
          
          selectedIndex = index;
          selectedQuestion = questions[index];
          
          const display = document.getElementById('questionDisplay');
          display.textContent = 'Question: ' + selectedQuestion.question + '\\n\\nOptions:\\n';
          Object.entries(selectedQuestion.options).forEach(([k, v]) => {
            display.textContent += k + ') ' + v + '\\n';
          });
          display.textContent += '\\nCorrect Answer: ' + selectedQuestion.correct_answer;
          
          document.getElementById('askBtn').disabled = false;
          log('Selected question ' + selectedQuestion.id);
        }
        
        async function askQuestion() {
          if (!selectedQuestion) return;
          
          const btn = document.getElementById('askBtn');
          const status = document.getElementById('status');
          const responseDiv = document.getElementById('response');
          
          btn.disabled = true;
          status.textContent = 'Status: Processing...';
          responseDiv.textContent = 'Processing...';
          
          log('Sending question to server...');
          
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
            
            log('Response received in ' + duration + 'ms');
            
            if (result.success) {
              responseDiv.textContent = result.response;
              
              // Check if correct
              const responseUpper = result.response.toUpperCase();
              const correct = selectedQuestion.correct_answer;
              const isCorrect = responseUpper.includes(correct + ')') || 
                               responseUpper.includes('ANSWER IS ' + correct) ||
                               responseUpper.includes('ANSWER: ' + correct);
              
              status.textContent = 'Status: Complete in ' + result.duration + 'ms | ' + 
                                  (isCorrect ? 'CORRECT ✓' : 'INCORRECT ✗');
              
              log('Answer validation: ' + (isCorrect ? 'CORRECT' : 'INCORRECT'));
              log('Response length: ' + result.response.length + ' characters');
              
              if (result.debugInfo) {
                log('Debug info: ' + JSON.stringify(result.debugInfo));
              }
            } else {
              status.textContent = 'Status: Error';
              responseDiv.textContent = 'Error: ' + result.error;
              log('Error: ' + result.error);
            }
          } catch (error) {
            status.textContent = 'Status: Error';
            responseDiv.textContent = 'Error: ' + error.message;
            log('Request error: ' + error.message);
          } finally {
            btn.disabled = false;
          }
        }
      </script>
    </body>
    </html>
  `);
});

// Simple logging to file
const logsDir = path.join(__dirname, 'gpqa-simple-logs');
if (!fs.existsSync(logsDir)) {
  fs.mkdirSync(logsDir, { recursive: true });
}

let currentLogFile = null;

function logToFile(type, message) {
  if (currentLogFile) {
    const entry = {
      timestamp: new Date().toISOString(),
      type,
      message
    };
    fs.appendFileSync(currentLogFile, JSON.stringify(entry) + '\n');
  }
}

// Intercept console methods
const originalLog = console.log;
const originalError = console.error;

console.log = (...args) => {
  originalLog(...args);
  logToFile('console.log', args.join(' '));
};

console.error = (...args) => {
  originalError(...args);
  logToFile('console.error', args.join(' '));
};

// API endpoint
app.post('/ask', async (req, res) => {
  const { question, questionId, correctAnswer } = req.body;
  
  // Create log file for this question
  const timestamp = new Date().toISOString().replace(/[:.]/g, '-');
  currentLogFile = path.join(logsDir, `${questionId}-${timestamp}.jsonl`);
  
  logToFile('start', `Processing question ${questionId}`);
  console.log(`Processing question ${questionId}`);
  console.log(`Log file: ${currentLogFile}`);
  
  try {
    console.log('Calling workflow engine...');
    const startTime = Date.now();
    
    const response = await workflowEngine.answer(question, { timeout: 300000 });
    
    const duration = Date.now() - startTime;
    console.log(`Response received in ${duration}ms`);
    console.log(`Response length: ${response.content.length} characters`);
    
    // Log full response
    logToFile('response', response.content);
    
    // Check correctness
    if (correctAnswer) {
      const responseUpper = response.content.toUpperCase();
      const isCorrect = responseUpper.includes(correctAnswer + ')') || 
                       responseUpper.includes('ANSWER IS ' + correctAnswer) ||
                       responseUpper.includes('ANSWER: ' + correctAnswer);
      
      logToFile('validation', `Correct answer: ${correctAnswer}, Model correct: ${isCorrect}`);
    }
    
    res.json({
      success: true,
      response: response.content,
      duration: duration,
      logFile: path.basename(currentLogFile)
    });
    
  } catch (error) {
    console.error('Error:', error);
    logToFile('error', error.stack || error.message);
    
    res.json({
      success: false,
      error: error.message
    });
  }
});

const PORT = 3010;
app.listen(PORT, () => {
  console.log(`GPQA Test Simple running on http://localhost:${PORT}`);
  console.log(`Logs will be saved to: ${logsDir}`);
});