const { spawn } = require('child_process');
const path = require('path');

console.log('Starting code execution service...');

const codeService = spawn('npm', ['start'], {
  cwd: path.join(__dirname, 'microservices', 'run-code'),
  shell: true
});

codeService.stdout.on('data', (data) => {
  console.log(`[Code Service] ${data}`);
});

codeService.stderr.on('data', (data) => {
  console.error(`[Code Service Error] ${data}`);
});

codeService.on('error', (error) => {
  console.error('Failed to start code service:', error);
});

// Wait a bit for service to start
setTimeout(() => {
  console.log('\nCode service should be running. Testing health check...');
  
  fetch('http://localhost:3001/health')
    .then(res => res.json())
    .then(data => {
      console.log('✅ Service is healthy:', data);
      console.log('\nYou can now run tests in another terminal:');
      console.log('  node test-code-integration.js');
    })
    .catch(err => {
      console.error('❌ Health check failed:', err.message);
    });
}, 3000);

// Keep process running
process.stdin.resume();