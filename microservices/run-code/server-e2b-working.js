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
  const { language, source, timeout = 30000 } = req.body;
  
  // Validate input
  if (!language || !source) {
    return res.status(400).json({
      error: 'Missing required fields: language and source'
    });
  }
  
  let sandbox;
  
  try {
    console.log(`Creating E2B sandbox for ${language}...`);
    
    // Create sandbox with API key and custom template (if available)
    const sandboxOptions = { apiKey: process.env.E2B_API_KEY };
    if (process.env.E2B_TEMPLATE_ID) {
      sandboxOptions.template = process.env.E2B_TEMPLATE_ID;
    }
    sandbox = await Sandbox.create(sandboxOptions);
    
    console.log('Sandbox created, executing code...');
    
    let result;
    
    console.log(`Running ${language} code (${source.length} chars)`);
    console.log('Code preview:', source.substring(0, 100));
    
    if (language === 'python') {
      console.log('Executing Python code...');
      try {
        // First try to write the file using E2B's file API
        await sandbox.files.write('/tmp/script.py', source);
        console.log('File written successfully, executing...');
        result = await sandbox.commands.run('python3 /tmp/script.py');
      } catch (writeError) {
        console.log('File write failed, trying alternative method:', writeError.message);
        // Fallback to python -c for simple scripts
        const pythonCmd = `python3 -c "${source.replace(/"/g, '\\"').replace(/\n/g, '\\n')}"`;
        result = await sandbox.commands.run(pythonCmd);
      }
    } else if (language === 'javascript') {
      // Execute JavaScript using node -e
      const jsCmd = `node -e "${source.replace(/"/g, '\\"').replace(/\n/g, '\\n')}"`;
      console.log('Executing JavaScript...');
      result = await sandbox.commands.run(jsCmd);
    } else if (language === 'bash') {
      // Execute bash commands directly
      console.log('Executing Bash...');
      result = await sandbox.commands.run(source);
    } else {
      throw new Error(`Unsupported language: ${language}`);
    }
    
    console.log('Execution result:', {
      stdout: result.stdout ? `${result.stdout.length} chars: ${result.stdout.substring(0, 100)}...` : '(empty)',
      stderr: result.stderr ? `${result.stderr.length} chars: ${result.stderr.substring(0, 100)}...` : '(empty)',
      exitCode: result.exitCode
    });
    
    res.json({
      stdout: result.stdout || '',
      stderr: result.stderr || '',
      exitCode: result.exitCode || 0
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
    if (sandbox && sandbox.kill) {
      try {
        await sandbox.kill();
        console.log('Sandbox killed');
      } catch (killError) {
        console.error('Error killing sandbox:', killError);
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
  console.log(`E2B Code execution service (working) running on port ${PORT}`);
  console.log(`Health check: http://localhost:3001/health`);
  console.log(`E2B API Key: ${process.env.E2B_API_KEY ? 'Configured ✓' : 'Missing ✗'}`);
});

module.exports = app;