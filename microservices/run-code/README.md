# Run Code Microservice

A lightweight microservice for executing code in a secure E2B sandbox environment.

## Setup

1. Install dependencies:
```bash
npm install
```

2. Set E2B API key:
```bash
export E2B_API_KEY="your-e2b-api-key"
```

3. Start the service:
```bash
npm start
```

## API

### POST /run_code

Execute code in a sandboxed environment.

**Request:**
```json
{
  "language": "python" | "javascript" | "bash",
  "source": "print('Hello World')",
  "timeout": 6000
}
```

**Response:**
```json
{
  "stdout": "Hello World\n",
  "stderr": "",
  "exitCode": 0
}
```

### GET /health

Health check endpoint.

**Response:**
```json
{
  "status": "healthy",
  "timestamp": "2025-07-30T12:00:00.000Z"
}
```

## Supported Languages

- **Python**: Direct execution in E2B Python environment
- **JavaScript**: Executed via Node.js in sandbox
- **Bash**: Shell commands executed safely

## Safety Features

- **Timeout**: Maximum 10 seconds execution time
- **Resource limits**: 512MB RAM, 0.5 vCPU (E2B defaults)
- **No internet access**: Sandbox has no outbound connectivity
- **Automatic cleanup**: Sandbox destroyed after each execution

## Deployment

### Fly.io (Recommended)

1. Install Fly CLI and login
2. Create app: `fly apps create run-code-service`
3. Set environment: `fly secrets set E2B_API_KEY="your-key"`
4. Deploy: `fly deploy`

Cost: ~$1/month for minimal usage

### Local Development

```bash
# Start service
npm run dev

# Test endpoint
curl -X POST http://localhost:3001/run_code \
  -H "Content-Type: application/json" \
  -d '{"language":"python","source":"print(2+2)"}'
```