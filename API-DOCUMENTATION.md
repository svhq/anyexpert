# AskAnyExpert API Documentation

## Base URL
```
http://localhost:3000
```

## Endpoints

### 1. Health Check
Check if the API server is running.

**Endpoint:** `GET /health`

**Response:**
```json
{
  "status": "healthy",
  "service": "Ask Any Expert API"
}
```

### 2. Ask Question (Main Endpoint)
Submit a question to the AI expert system.

**Endpoint:** `POST /api/ask`

**Headers:**
```
Content-Type: application/json
```

**Request Body:**
```json
{
  "question": "Your question here"
}
```

**Response:**
```json
{
  "answer": "The detailed answer to your question...",
  "sources": [
    {
      "title": "Source Title",
      "url": "https://example.com",
      "snippet": "Relevant snippet from the source"
    }
  ],
  "searchPerformed": true,
  "agent": "expert",
  "confidence": 0.95,
  "metadata": {
    "agent": "expert",
    "finalConfidence": 0.95,
    "steps": 2,
    "searchQueries": ["query1", "query2"],
    "toolsUsed": ["search", "reason"]
  }
}
```

**Response Fields:**
- `answer` (string): The main answer to the question
- `sources` (array): List of sources used (if web search was performed)
- `searchPerformed` (boolean): Whether web search was used
- `agent` (string): Which expert agent was used
- `confidence` (number): Confidence score (0-1)
- `metadata` (object): Additional processing information

**Error Response:**
```json
{
  "error": "Failed to process question",
  "message": "Specific error details"
}
```

## Example Usage

### JavaScript/React Example:
```javascript
async function askQuestion(question) {
  try {
    const response = await fetch('http://localhost:3000/api/ask', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ question })
    });
    
    const data = await response.json();
    return data;
  } catch (error) {
    console.error('Error:', error);
    throw error;
  }
}

// Usage
const result = await askQuestion("What are the latest AI developments?");
console.log(result.answer);
```

### cURL Example:
```bash
curl -X POST http://localhost:3000/api/ask \
  -H "Content-Type: application/json" \
  -d '{"question": "What is quantum computing?"}'
```

## Running the API Server

1. Install dependencies:
```bash
cd askanyexpert
npm install
```

2. Set up environment variables in `.env`:
```env
OPENROUTER_API_KEY=your_key
SERPER_API_KEY=your_key
E2B_API_KEY=your_key
USE_MODULAR_SYSTEM=true
API_MODE=full
```

3. Start the server:
```bash
node api-server.js
```

The server will start on `http://localhost:3000`

## CORS Configuration
The API has CORS enabled, allowing requests from any origin. This makes it easy to connect from your Lovable React frontend.

## Rate Limiting
Currently no rate limiting is implemented. Consider adding for production use.

## Authentication
Currently no authentication required. Consider adding API keys for production.