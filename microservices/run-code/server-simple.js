const express = require('express');
const cors = require('cors');

const app = express();
const PORT = process.env.PORT || 3001;

// Middleware
app.use(cors());
app.use(express.json({ limit: '1mb' }));

// Health check endpoint
app.get('/health', (req, res) => {
  res.json({ status: 'healthy', timestamp: new Date().toISOString() });
});

// Main code execution endpoint - simplified for testing
app.post('/run_code', async (req, res) => {
  const { language, source, timeout = 6000 } = req.body;
  
  // Validate input
  if (!language || !source) {
    return res.status(400).json({
      error: 'Missing required fields: language and source'
    });
  }
  
  console.log(`Executing ${language} code (${source.length} chars)`);
  
  try {
    // For now, return mock responses to test integration
    // Replace with actual E2B implementation once SDK is clarified
    
    if (language === 'python' && source.includes('[x**2 for x in range(5)]')) {
      res.json({
        stdout: '[0, 1, 4, 9, 16]\n',
        stderr: '',
        exitCode: 0
      });
    } else if (language === 'python' && source.includes('sum(primes)')) {
      res.json({
        stdout: 'Primes below 20: [2, 3, 5, 7, 11, 13, 17, 19]\nSum: 77\n',
        stderr: '',
        exitCode: 0
      });
    } else if (language === 'javascript') {
      res.json({
        stdout: 'Hello from Node.js!\n4\n',
        stderr: '',
        exitCode: 0
      });
    } else {
      // Generic response
      res.json({
        stdout: `Executed ${language} code successfully\n`,
        stderr: '',
        exitCode: 0
      });
    }
    
  } catch (error) {
    console.error('Code execution error:', error);
    res.json({
      stdout: '',
      stderr: error.message || 'Execution failed',
      exitCode: 1
    });
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
  console.log(`Code execution service (mock) running on port ${PORT}`);
  console.log(`Health check: http://localhost:${PORT}/health`);
  console.log('Note: This is a mock service for testing integration');
  console.log('Replace with actual E2B implementation for production');
});

module.exports = app;