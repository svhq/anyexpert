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
  
  let sandbox;
  
  try {
    console.log(`Creating E2B sandbox for ${language}...`);
    
    // Create sandbox with the E2B API key
    sandbox = await Sandbox.create();
    
    console.log('Sandbox created, executing code...');
    
    let result;
    
    if (language === 'python') {
      // Write Python code to a file and execute it
      await sandbox.filesystem.write('/tmp/script.py', source);
      result = await sandbox.process.start({
        cmd: 'python3 /tmp/script.py',
        timeout: timeout / 1000
      });
      await result.wait();
    } else if (language === 'javascript') {
      // Write JS code to a file and execute with node
      await sandbox.filesystem.write('/tmp/script.js', source);
      result = await sandbox.process.start({
        cmd: 'node /tmp/script.js',
        timeout: timeout / 1000
      });
      await result.wait();
    } else if (language === 'bash') {
      // Execute bash directly
      result = await sandbox.process.start({
        cmd: source,
        timeout: timeout / 1000
      });
      await result.wait();
    } else {
      throw new Error(`Unsupported language: ${language}`);
    }
    
    console.log('Execution complete');
    
    res.json({
      stdout: result.stdout || '',
      stderr: result.stderr || '',
      exitCode: result.exitCode || 0
    });
    
  } catch (error) {
    console.error('E2B execution error:', error);
    
    // Try alternative approach using code-interpreter template
    try {
      // If basic execution fails, try with code-interpreter template
      if (sandbox) await sandbox.close();
      
      // Use code-interpreter specific template
      sandbox = await Sandbox.create({ template: 'base' });
      
      let execResult;
      if (language === 'python') {
        // Execute Python using direct command
        const pythonCmd = `python3 -c "${source.replace(/"/g, '\\"').replace(/\n/g, '\\n')}"`;
        execResult = await sandbox.process.start({ cmd: pythonCmd });
        await execResult.wait();
      } else if (language === 'javascript') {
        // Execute JS using direct command
        const jsCmd = `node -e "${source.replace(/"/g, '\\"').replace(/\n/g, '\\n')}"`;
        execResult = await sandbox.process.start({ cmd: jsCmd });
        await execResult.wait();
      } else {
        execResult = await sandbox.process.start({ cmd: source });
        await execResult.wait();
      }
      
      res.json({
        stdout: execResult.stdout || '',
        stderr: execResult.stderr || '',
        exitCode: execResult.exitCode || 0
      });
      
    } catch (fallbackError) {
      console.error('Fallback execution also failed:', fallbackError);
      
      res.json({
        stdout: '',
        stderr: error.message || 'Execution failed',
        exitCode: 1
      });
    }
    
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
  console.log(`E2B Code execution service (v1 API) running on port ${PORT}`);
  console.log(`Health check: http://localhost:3001/health`);
  console.log(`E2B API Key: ${process.env.E2B_API_KEY ? 'Configured ✓' : 'Missing ✗'}`);
});

module.exports = app;