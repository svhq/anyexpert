# E2B Code Execution Service Guide

## How It Works

The Ask Any Expert system uses **E2B.dev** for secure code execution in sandboxed environments. Here's how it works:

1. **Service Architecture**:
   - The E2B service runs as a separate microservice on port 3001
   - It provides a REST API for code execution requests
   - The math-agent and code-agent communicate with this service
   - E2B creates isolated sandboxes for each code execution

2. **When It's Needed**:
   - Mathematical calculations requiring actual computation
   - Code generation with execution requests
   - Data analysis and visualization tasks
   - Any query explicitly asking to "run" or "execute" code

3. **Service Lifecycle**:
   - The service does NOT auto-start with the main application
   - It remains running until you stop it (Ctrl+C) or close the terminal
   - It does NOT auto-shutdown after queries

## Starting the Service

### Method 1: Manual Start (Recommended for Development)
```bash
# In one terminal, start the E2B service
node start-e2b-service.js

# In another terminal, run your main application
node index.js
```

### Method 2: Using Batch Files (Windows)
```bash
# Double-click or run:
start-code-service.bat
```

### Method 3: Combined Start (All Services)
```bash
node start-services.js
```

## Checking Service Status

### Quick Health Check
```bash
curl http://localhost:3001/health
```

Expected response:
```json
{
  "status": "healthy",
  "timestamp": "2025-07-31T02:14:51.023Z",
  "e2bKey": "configured"
}
```

### From Node.js
```javascript
fetch('http://localhost:3001/health')
  .then(res => res.json())
  .then(data => console.log('Service status:', data))
  .catch(err => console.error('Service is down:', err.message));
```

## Ensuring It's Always Working

### Option 1: Process Manager (PM2) - RECOMMENDED
Install PM2 globally:
```bash
npm install -g pm2
```

Start the service with PM2:
```bash
pm2 start microservices/run-code/server-e2b-working.js --name e2b-service
pm2 save
pm2 startup  # Makes it start on system boot
```

Check status:
```bash
pm2 status e2b-service
pm2 logs e2b-service
```

### Option 2: Windows Service (Windows Only)
Use `node-windows` to install as a Windows service:
```bash
npm install -g node-windows
# Then create a service installer script
```

### Option 3: Docker Container
```dockerfile
FROM node:18
WORKDIR /app
COPY microservices/run-code .
RUN npm install
EXPOSE 3001
CMD ["node", "server-e2b-working.js"]
```

### Option 4: Auto-Start Script
Create a launcher that checks and starts the service if needed:
```javascript
// auto-start-e2b.js
const { spawn } = require('child_process');

async function checkService() {
  try {
    const response = await fetch('http://localhost:3001/health');
    const data = await response.json();
    console.log('✅ E2B service already running:', data.status);
    return true;
  } catch (err) {
    console.log('⚠️  E2B service not running, starting it...');
    return false;
  }
}

async function startService() {
  const service = spawn('node', ['start-e2b-service.js'], {
    detached: true,
    stdio: 'ignore'
  });
  service.unref();
  
  // Wait for service to start
  await new Promise(resolve => setTimeout(resolve, 5000));
  
  // Verify it started
  const running = await checkService();
  if (!running) {
    throw new Error('Failed to start E2B service');
  }
}

async function ensureServiceRunning() {
  const isRunning = await checkService();
  if (!isRunning) {
    await startService();
  }
}

// Auto-check on require
ensureServiceRunning();

module.exports = { ensureServiceRunning };
```

## Integration with Main App

### Automatic Service Check
Add to your main `index.js`:
```javascript
// At the top of index.js
const { ensureServiceRunning } = require('./auto-start-e2b');

// Before starting the main app
async function start() {
  await ensureServiceRunning();
  // ... rest of your app startup
}
```

### Graceful Fallback
The system already handles missing E2B service gracefully:
- If code execution fails, it falls back to showing code without execution
- Math questions still get answered with formulas and explanations
- No queries fail completely due to missing E2B service

## Troubleshooting

### Service Won't Start
1. Check if port 3001 is already in use:
   ```bash
   netstat -an | findstr :3001
   ```

2. Verify E2B API key is set:
   ```bash
   echo %E2B_API_KEY%
   ```

3. Check logs:
   ```bash
   cd microservices/run-code
   npm start
   ```

### Connection Errors
- Error: `ECONNREFUSED` - Service is not running
- Error: `ETIMEDOUT` - Service is overloaded or network issue
- Error: `401 Unauthorized` - E2B API key is invalid

### Performance Issues
- First execution may be slow (sandbox creation)
- Subsequent executions are faster (warm sandboxes)
- Complex computations may take longer (normal)

## Best Practices

1. **For Development**:
   - Start service manually when needed
   - Use two terminal windows (one for E2B, one for main app)

2. **For Production**:
   - Use PM2 or system service manager
   - Implement health checks and auto-restart
   - Monitor service logs

3. **For Testing**:
   - Always check service health before running tests
   - Use the test scripts to verify functionality

## Quick Reference

```bash
# Start service
node start-e2b-service.js

# Check health
curl http://localhost:3001/health

# Test execution
node test-code-execution-manual.js

# View logs (if using PM2)
pm2 logs e2b-service

# Stop service (if using PM2)
pm2 stop e2b-service
```

## Environment Variables

Required in `.env`:
```
E2B_API_KEY=your_api_key_here
CODE_EXECUTION_URL=http://localhost:3001  # Optional, this is default
```

## Summary

- **No**, you don't need to start E2B every time before a query
- **No**, it doesn't switch off by itself
- **Yes**, you need to start it once per session
- **Best practice**: Use PM2 to keep it running permanently