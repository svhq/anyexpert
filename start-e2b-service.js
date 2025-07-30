const { spawn } = require('child_process');
const path = require('path');

// Load environment variables
require('dotenv').config();

console.log('Starting E2B code execution service...');
console.log('E2B API Key:', process.env.E2B_API_KEY ? 'Found ✓' : 'Missing ✗');

const service = spawn('npm', ['start'], {
  cwd: path.join(__dirname, 'microservices', 'run-code'),
  shell: true,
  env: { ...process.env }
});

service.stdout.on('data', (data) => {
  console.log(`[E2B Service] ${data.toString().trim()}`);
});

service.stderr.on('data', (data) => {
  console.error(`[E2B Error] ${data.toString().trim()}`);
});

service.on('error', (error) => {
  console.error('Failed to start service:', error);
});

// Wait for service to start
setTimeout(() => {
  console.log('\nTesting E2B service...');
  
  fetch('http://localhost:3001/health')
    .then(res => res.json())
    .then(data => {
      console.log('✅ Service health:', data);
      console.log('\nReady to execute code with real E2B sandboxes!');
      console.log('Test with: node test-code-integration.js');
    })
    .catch(err => {
      console.error('❌ Health check failed:', err.message);
    });
}, 3000);

// Keep process running
process.stdin.resume();