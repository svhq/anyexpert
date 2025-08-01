const express = require('express');
const cors = require('cors');
const { Sandbox } = require('e2b');
const path = require('path');

// Load environment variables from parent directory
require('dotenv').config({ path: path.join(__dirname, '../../.env') });

const app = express();
const PORT = process.env.PORT || 3001;

// Global sandbox instance
let globalSandbox = null;
let sandboxLastUsed = Date.now();

// Middleware
app.use(cors());
app.use(express.json({ limit: '1mb' }));

// Get or create sandbox
async function getSandbox() {
  // Check if we need a new sandbox (none exists, closed, or idle too long)
  if (!globalSandbox || globalSandbox.closed || (Date.now() - sandboxLastUsed > 3500000)) {
    console.log('Creating new sandbox...');
    if (globalSandbox && !globalSandbox.closed) {
      try {
        await globalSandbox.kill();
      } catch (e) {
        console.error('Error killing old sandbox:', e);
      }
    }
    
    const sandboxOptions = { 
      apiKey: process.env.E2B_API_KEY,
      timeout: 3600000 // 1 hour
    };
    
    if (process.env.E2B_TEMPLATE_ID) {
      sandboxOptions.template = process.env.E2B_TEMPLATE_ID;
      console.log('Using template:', process.env.E2B_TEMPLATE_ID);
    }
    
    globalSandbox = await Sandbox.create(sandboxOptions);
    console.log('Sandbox created successfully');
  } else {
    // Reset timeout
    console.log('Reusing existing sandbox');
    await globalSandbox.setTimeout(3600000);
  }
  
  sandboxLastUsed = Date.now();
  return globalSandbox;
}

// Health check endpoint with sandbox verification
app.get('/health', async (req, res) => {
  try {
    const sandbox = await getSandbox();
    // Quick test to verify sandbox is working
    const result = await sandbox.commands.run('echo "healthy"');
    
    res.json({ 
      status: result.exitCode === 0 ? 'healthy' : 'unhealthy',
      timestamp: new Date().toISOString(),
      e2bKey: process.env.E2B_API_KEY ? 'configured' : 'missing',
      templateId: process.env.E2B_TEMPLATE_ID || 'none',
      sandboxStatus: 'active',
      sandboxAge: Date.now() - sandboxLastUsed
    });
  } catch (error) {
    res.status(500).json({ 
      status: 'unhealthy',
      error: 'sandbox-error',
      message: error.message
    });
  }
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
  
  try {
    console.log(`Executing ${language} code...`);
    
    // Get or reuse sandbox
    const sandbox = await getSandbox();
    
    let result;
    
    console.log(`Running ${language} code (${source.length} chars)`);
    console.log('Code preview:', source.substring(0, 100));
    console.log('Template ID:', process.env.E2B_TEMPLATE_ID || 'default');
    
    if (language === 'python') {
      console.log('Executing Python code...');
      console.log('Template:', process.env.E2B_TEMPLATE_ID || 'default');
      
      // Install libraries dynamically if needed (even with custom template)
      // Check if libraries are already installed in the template
      console.log('Checking for numpy...');
      let numpyInstalled = false;
      try {
        const checkCmd = await sandbox.commands.run('python -c "import numpy"');
        numpyInstalled = checkCmd.exitCode === 0;
      } catch (checkError) {
        // Numpy not installed, which is fine
        numpyInstalled = false;
      }
      console.log('Numpy installed:', numpyInstalled);
      
      if (!numpyInstalled) {
        // Check if code needs scientific libraries and install if needed
        const needsScientific = source.includes('numpy') || source.includes('np.') || 
                               source.includes('scipy') || source.includes('sympy') ||
                               source.includes('pandas') || source.includes('pd.');
        
        const needsChemistry = source.includes('rdkit') || source.includes('Chem') || 
                              source.includes('rdchiral') || source.includes('rdChemReactions') ||
                              source.includes('SMILES') || source.includes('Diels-Alder');
        
        if (needsScientific) {
          console.log('Installing scientific libraries...');
          const installCmd = 'pip install -q numpy scipy sympy pandas matplotlib';
          const installResult = await sandbox.commands.run(installCmd);
          console.log('Install result:', installResult.exitCode);
          if (installResult.exitCode !== 0) {
            console.log('Install stderr:', installResult.stderr);
          }
        }
        
        if (needsChemistry) {
          console.log('Installing chemistry libraries (RDKit, RDChiral)...');
          // RDKit needs numpy<2 for now due to compatibility
          const chemCmd = 'pip install -q "numpy<2" rdkit-pypi rdchiral';
          const chemResult = await sandbox.commands.run(chemCmd);
          console.log('Chemistry install result:', chemResult.exitCode);
          if (chemResult.exitCode !== 0) {
            console.log('Chemistry install stderr:', chemResult.stderr);
          }
        }
      } else {
        console.log('Libraries already installed in template');
      }
      
      // Write Python code to a file and execute it
      const escaped = source.replace(/'/g, "'\\''");
      const writeCmd = `cat > /tmp/script.py << 'EOF'
${source}
EOF`;
      
      console.log('Writing Python script...');
      const writeResult = await sandbox.commands.run(writeCmd);
      if (writeResult.exitCode !== 0) {
        console.log('Write failed:', writeResult.stderr);
      }
      
      console.log('Executing Python script...');
      result = await sandbox.commands.run('python /tmp/script.py');
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
    console.error('Error stack:', error.stack);
    
    // If sandbox error, reset it
    if (error.message && error.message.includes('sandbox')) {
      globalSandbox = null;
    }
    
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
  console.log(`E2B Code execution service (custom template) running on port ${PORT}`);
  console.log(`Health check: http://localhost:3001/health`);
  console.log(`E2B API Key: ${process.env.E2B_API_KEY ? 'Configured ✓' : 'Missing ✗'}`);
  console.log(`E2B Template: ${process.env.E2B_TEMPLATE_ID || 'Default'}`);
  console.log(`Sandbox mode: Warm (reusing sandbox across requests)`);
});

// Graceful shutdown
process.on('SIGINT', async () => {
  console.log('\nShutting down E2B service...');
  if (globalSandbox && !globalSandbox.closed) {
    try {
      await globalSandbox.kill();
      console.log('Sandbox cleaned up');
    } catch (e) {
      console.error('Error cleaning up sandbox:', e);
    }
  }
  process.exit(0);
});

module.exports = app;