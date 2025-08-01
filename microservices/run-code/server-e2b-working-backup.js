const express = require('express');
const cors = require('cors');
const { Sandbox } = require('e2b');
const path = require('path');

// Load environment variables from parent directory
require('dotenv').config({ path: path.join(__dirname, '../../.env') });

const app = express();
const PORT = process.env.PORT || 3002;

// Middleware
app.use(cors());
app.use(express.json({ limit: '1mb' }));

// Global sandbox instance
let sandbox = null;
let sandboxReady = false;

// Initialize sandbox on startup
async function initializeSandbox() {
  try {
    console.log('Initializing persistent E2B sandbox...');
    
    const sandboxOptions = { 
      apiKey: process.env.E2B_API_KEY,
      onStdout: (output) => console.log('Sandbox output:', output),
      onStderr: (error) => console.error('Sandbox error:', error)
    };
    
    if (process.env.E2B_TEMPLATE_ID) {
      sandboxOptions.template = process.env.E2B_TEMPLATE_ID;
    }
    
    sandbox = await Sandbox.create(sandboxOptions);
    console.log('Sandbox created, installing libraries...');
    
    // Run the bootstrap command to install libraries
    const bootstrapCmd = `apt-get update -y && DEBIAN_FRONTEND=noninteractive apt-get install -y build-essential git curl ca-certificates libgl1 libglib2.0-0 tzdata && python -m pip install -q --upgrade pip && pip install -q numpy pandas sympy requests beautifulsoup4 lxml matplotlib seaborn networkx scikit-learn pillow scipy statsmodels plotly mpmath`;
    
    console.log('Running bootstrap command...');
    const bootstrapResult = await sandbox.commands.run(bootstrapCmd);
    console.log('Bootstrap complete:', bootstrapResult.exitCode === 0 ? 'Success' : 'Failed');
    
    if (bootstrapResult.exitCode !== 0) {
      console.error('Bootstrap stderr:', bootstrapResult.stderr);
    }
    
    // Test that libraries are installed
    const testResult = await sandbox.commands.run('python3 -c "import numpy, pandas, sympy; print(\\"Libraries ready\\")"');
    if (testResult.exitCode === 0) {
      console.log('✅ Libraries verified:', testResult.stdout);
      sandboxReady = true;
    } else {
      console.error('❌ Library test failed:', testResult.stderr);
    }
    
  } catch (error) {
    console.error('Failed to initialize sandbox:', error);
    sandbox = null;
    sandboxReady = false;
  }
}

// Initialize sandbox on startup
initializeSandbox();

// Health check endpoint
app.get('/health', (req, res) => {
  res.json({ 
    status: sandboxReady ? 'healthy' : 'initializing', 
    timestamp: new Date().toISOString(),
    e2bKey: process.env.E2B_API_KEY ? 'configured' : 'missing',
    sandboxReady
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
  
  // Check if sandbox is ready
  if (!sandbox || !sandboxReady) {
    return res.json({
      stdout: '',
      stderr: 'Sandbox not ready. Please try again in a few seconds.',
      exitCode: 1
    });
  }
  
  try {
    console.log(`Running ${language} code (${source.length} chars)`);
    console.log('Code preview:', source.substring(0, 100));
    
    let result;
    
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
  }
});

// Restart sandbox endpoint
app.post('/restart_sandbox', async (req, res) => {
  try {
    if (sandbox && sandbox.kill) {
      await sandbox.kill();
    }
    sandboxReady = false;
    await initializeSandbox();
    res.json({ status: 'restarted', sandboxReady });
  } catch (error) {
    res.status(500).json({ error: error.message });
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

// Graceful shutdown
process.on('SIGINT', async () => {
  console.log('\nShutting down gracefully...');
  if (sandbox && sandbox.kill) {
    try {
      await sandbox.kill();
      console.log('Sandbox killed');
    } catch (error) {
      console.error('Error killing sandbox:', error);
    }
  }
  process.exit(0);
});

// Start server
app.listen(PORT, () => {
  console.log(`E2B Code execution service (persistent) running on port ${PORT}`);
  console.log(`Health check: http://localhost:${PORT}/health`);
  console.log(`E2B API Key: ${process.env.E2B_API_KEY ? 'Configured ✓' : 'Missing ✗'}`);
});

module.exports = app;