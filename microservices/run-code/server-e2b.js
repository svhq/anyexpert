const express = require('express');
const cors = require('cors');
const { Sandbox } = require('e2b');
const path = require('path');

// Load environment variables from parent directory
require('dotenv').config({ path: path.join(__dirname, '../../.env') });

const app = express();
const PORT = process.env.PORT || 3001;

// Middleware
app.use(cors());
app.use(express.json({ limit: '1mb' }));

// Health check endpoint
app.get('/health', (req, res) => {
  res.json({ 
    status: 'healthy', 
    timestamp: new Date().toISOString(),
    e2bKey: process.env.E2B_API_KEY ? 'configured' : 'missing'
  });
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
  
  let sandbox;
  
  try {
    console.log(`Creating E2B sandbox for ${language}...`);
    
    // Create sandbox with the E2B API key
    sandbox = await Sandbox.create({
      apiKey: process.env.E2B_API_KEY
    });
    
    console.log('Sandbox created, executing code...');
    
    let stdout = '';
    let stderr = '';
    
    if (language === 'python') {
      // Execute Python code
      const result = await sandbox.runPython(source);
      stdout = result.stdout || '';
      stderr = result.stderr || '';
    } else if (language === 'javascript') {
      // Execute JavaScript via Node.js in sandbox
      const jsCode = `node -e "${source.replace(/"/g, '\\"')}"`;
      const result = await sandbox.runBash(jsCode);
      stdout = result.stdout || '';
      stderr = result.stderr || '';
    } else if (language === 'bash') {
      // Execute bash commands
      const result = await sandbox.runBash(source);
      stdout = result.stdout || '';
      stderr = result.stderr || '';
    }
    
    console.log('Execution complete');
    
    res.json({
      stdout,
      stderr,
      exitCode: stderr ? 1 : 0
    });
    
  } catch (error) {
    console.error('E2B execution error:', error);
    
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
        console.log('Sandbox closed');
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
  console.log(`E2B Code execution service running on port ${PORT}`);
  console.log(`Health check: http://localhost:${PORT}/health`);
  console.log(`E2B API Key: ${process.env.E2B_API_KEY ? 'Configured ✓' : 'Missing ✗'}`);
});

module.exports = app;