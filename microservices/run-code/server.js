const express = require('express');
const cors = require('cors');
const { Sandbox } = require('e2b');

const app = express();
const PORT = process.env.PORT || 3001;

// Middleware
app.use(cors());
app.use(express.json({ limit: '1mb' }));

// Health check endpoint
app.get('/health', (req, res) => {
  res.json({ status: 'healthy', timestamp: new Date().toISOString() });
});

// Main code execution endpoint
app.post('/run_code', async (req, res) => {
  const { language, source, timeout = 6000 } = req.body;
  
  // Validate input
  if (!language || !source) {
    return res.status(400).json({
      error: 'Missing required fields: language and source'
    });
  }
  
  // Validate supported languages
  const supportedLanguages = ['python', 'javascript', 'bash'];
  if (!supportedLanguages.includes(language.toLowerCase())) {
    return res.status(400).json({
      error: `Unsupported language: ${language}. Supported: ${supportedLanguages.join(', ')}`
    });
  }
  
  // Validate timeout (max 10 seconds for safety)
  const maxTimeout = 10000;
  const actualTimeout = Math.min(timeout, maxTimeout);
  
  let sandbox;
  
  try {
    // Create E2B sandbox with Python environment
    sandbox = await Sandbox.create({
      apiKey: process.env.E2B_API_KEY,
      template: 'base'  // Base template supports multiple languages
    });
    
    let result;
    
    if (language.toLowerCase() === 'python') {
      // Execute Python code directly
      result = await sandbox.runCode({
        code: source,
        language: 'python',
        timeout: actualTimeout / 1000  // E2B uses seconds
      });
    } else if (language.toLowerCase() === 'javascript') {
      // Execute JavaScript via Node.js
      result = await sandbox.runCode({
        code: source,
        language: 'javascript',
        timeout: actualTimeout / 1000
      });
    } else if (language.toLowerCase() === 'bash') {
      // Execute bash commands
      result = await sandbox.runBash({
        commands: source.split('\n'),
        timeout: actualTimeout / 1000
      });
    }
    
    // Format response
    res.json({
      stdout: result.logs?.stdout || result.stdout || '',
      stderr: result.logs?.stderr || result.stderr || result.error || '',
      exitCode: result.error ? 1 : 0
    });
    
  } catch (error) {
    console.error('Code execution error:', error);
    
    res.json({
      stdout: '',
      stderr: error.message || 'Execution failed',
      exitCode: 1
    });
    
  } finally {
    // Always clean up the sandbox
    if (sandbox) {
      try {
        await sandbox.close();
      } catch (closeError) {
        console.error('Error closing sandbox:', closeError);
      }
    }
  }
});

// Error handling middleware
app.use((error, req, res, next) => {
  console.error('Server error:', error);
  res.status(500).json({
    error: 'Internal server error',
    stdout: '',
    stderr: error.message,
    exitCode: 1
  });
});

// Start server
app.listen(PORT, () => {
  console.log(`Code execution service running on port ${PORT}`);
  console.log(`Health check: http://localhost:${PORT}/health`);
  console.log(`E2B API Key: ${process.env.E2B_API_KEY ? 'Configured' : 'Missing!'}`);
});

module.exports = app;