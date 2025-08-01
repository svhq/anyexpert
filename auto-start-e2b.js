const { spawn } = require('child_process');
const path = require('path');

// Load environment variables
require('dotenv').config();

/**
 * Check if the E2B service is running
 * @returns {Promise<boolean>} True if service is healthy
 */
async function checkService() {
  try {
    const response = await fetch('http://localhost:3001/health');
    const data = await response.json();
    console.log('✅ E2B service already running:', data.status);
    return true;
  } catch (err) {
    console.log('⚠️  E2B service not running');
    return false;
  }
}

/**
 * Start the E2B service in the background
 */
async function startService() {
  console.log('🚀 Starting E2B service...');
  
  const servicePath = path.join(__dirname, 'microservices', 'run-code');
  
  // Start the service
  const service = spawn('npm', ['start'], {
    cwd: servicePath,
    detached: true,
    stdio: ['ignore', 'pipe', 'pipe'],
    shell: true
  });
  
  // Optional: Log output for debugging
  service.stdout.on('data', (data) => {
    console.log(`[E2B] ${data.toString().trim()}`);
  });
  
  service.stderr.on('data', (data) => {
    console.error(`[E2B Error] ${data.toString().trim()}`);
  });
  
  // Detach the process so it continues running
  service.unref();
  
  // Wait for service to fully start
  console.log('⏳ Waiting for service to initialize...');
  await new Promise(resolve => setTimeout(resolve, 5000));
  
  // Verify it started successfully
  const maxRetries = 3;
  for (let i = 0; i < maxRetries; i++) {
    const running = await checkService();
    if (running) {
      console.log('🎉 E2B service started successfully!');
      return;
    }
    
    if (i < maxRetries - 1) {
      console.log(`⏳ Retry ${i + 1}/${maxRetries}...`);
      await new Promise(resolve => setTimeout(resolve, 2000));
    }
  }
  
  throw new Error('Failed to start E2B service after multiple attempts');
}

/**
 * Ensure the E2B service is running
 * @param {boolean} silent - If true, suppress console output
 */
async function ensureServiceRunning(silent = false) {
  try {
    if (!silent) console.log('🔍 Checking E2B service status...');
    
    const isRunning = await checkService();
    if (!isRunning) {
      await startService();
    } else if (!silent) {
      console.log('✅ E2B service is ready!');
    }
    
    return true;
  } catch (error) {
    console.error('❌ Error with E2B service:', error.message);
    console.log('💡 You can manually start it with: node start-e2b-service.js');
    return false;
  }
}

/**
 * Stop the E2B service (if you need to)
 */
async function stopService() {
  try {
    // This would need to track the process ID or use PM2
    console.log('⚠️  Stopping service requires process management (use PM2 or manual Ctrl+C)');
  } catch (error) {
    console.error('Error stopping service:', error.message);
  }
}

// Export functions
module.exports = {
  checkService,
  startService,
  ensureServiceRunning,
  stopService
};

// If run directly, ensure service is running
if (require.main === module) {
  ensureServiceRunning()
    .then(success => {
      if (success) {
        console.log('\n📝 You can now run queries that require code execution!');
        console.log('📊 Test with: node test-code-execution-manual.js');
      }
      process.exit(success ? 0 : 1);
    })
    .catch(err => {
      console.error('Unexpected error:', err);
      process.exit(1);
    });
}