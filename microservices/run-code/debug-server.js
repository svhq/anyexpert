// Debug wrapper for E2B server
const { spawn } = require('child_process');
const path = require('path');

console.log('Starting E2B server with debug logging...');

const server = spawn('node', ['server-e2b-working.js'], {
  cwd: __dirname,
  env: process.env,
  stdio: 'inherit'
});

server.on('error', (err) => {
  console.error('Failed to start server:', err);
});

server.on('exit', (code) => {
  console.log(`Server exited with code ${code}`);
});

// Keep the wrapper running
process.stdin.resume();