const express = require('express');
const workflowEngine = require('./src/workflow-engine');

const app = express();
app.use(express.json());

// Simple HTML UI
app.get('/', (req, res) => {
  res.send(`
    <!DOCTYPE html>
    <html>
    <head>
      <title>Simple Ask UI</title>
      <style>
        body { font-family: monospace; padding: 20px; }
        textarea { width: 100%; height: 100px; }
        button { padding: 10px 20px; margin: 10px 0; }
        #output { white-space: pre-wrap; background: #f0f0f0; padding: 10px; }
        .log { color: #666; font-size: 0.9em; }
      </style>
    </head>
    <body>
      <h1>Simple Ask UI</h1>
      <textarea id="question" placeholder="Enter your question...">What is 2 + 2?</textarea>
      <br>
      <button onclick="ask()">Ask</button>
      <div id="logs"></div>
      <h3>Response:</h3>
      <div id="output">No response yet</div>
      
      <script>
        async function ask() {
          const question = document.getElementById('question').value;
          const logs = document.getElementById('logs');
          const output = document.getElementById('output');
          
          logs.innerHTML = '<div class="log">Sending question...</div>';
          output.textContent = 'Processing...';
          
          try {
            const response = await fetch('/ask', {
              method: 'POST',
              headers: { 'Content-Type': 'application/json' },
              body: JSON.stringify({ question })
            });
            
            const result = await response.json();
            
            if (result.success) {
              output.textContent = result.response;
              logs.innerHTML += '<div class="log">Completed in ' + result.duration + 'ms</div>';
            } else {
              output.textContent = 'Error: ' + result.error;
            }
          } catch (error) {
            output.textContent = 'Error: ' + error.message;
          }
        }
      </script>
    </body>
    </html>
  `);
});

// Simple API endpoint
app.post('/ask', async (req, res) => {
  const { question } = req.body;
  
  console.log('Received question:', question);
  
  try {
    const startTime = Date.now();
    const response = await workflowEngine.answer(question, { timeout: 300000 });
    const duration = Date.now() - startTime;
    
    console.log('Response received in', duration, 'ms');
    
    res.json({
      success: true,
      response: response.content,
      duration: duration
    });
  } catch (error) {
    console.error('Error:', error);
    res.json({
      success: false,
      error: error.message
    });
  }
});

const PORT = 3006;
app.listen(PORT, () => {
  console.log(`Simple UI running on http://localhost:${PORT}`);
});