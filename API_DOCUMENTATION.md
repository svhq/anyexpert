# Ask Any Expert - API Documentation

## Base URL
```
http://localhost:8080
```

## Endpoints

### 1. Ask Question
**Endpoint:** `POST /api/ask`

**Description:** Submit a question to the AI expert system. Maintains conversation context through sessions.

**Request Body:**
```json
{
  "question": "string (required) - The question to ask",
  "sessionId": "string (optional) - Session ID for conversation continuity"
}
```

**Response:** `200 OK`
```json
{
  "answer": "string - The AI-generated response in markdown format",
  "sessionId": "string - Session ID to use for follow-up questions",
  "metadata": {
    "requestId": "string - Unique request identifier",
    "processingTime": "number - Time taken in milliseconds",
    "confidence": "number - Confidence score (0-1)",
    "searchesPerformed": "number - Number of web searches performed",
    "expertsConsulted": ["array - Expert personas used"],
    "contextUsed": "boolean - Whether conversation context was used",
    "isFollowUp": "boolean - Whether detected as follow-up question",
    "sources": [
      {
        "number": "number - Citation number",
        "title": "string - Source title",
        "url": "string - Source URL",
        "snippet": "string - Relevant excerpt"
      }
    ],
    "modernData": {
      "bayesianPosterior": "number - Bayesian confidence (0-1)",
      "contextType": "string - Type of context detected",
      "informationScores": ["array - Information overlap scores"],
      "attentionActivation": "object - Attention pattern data"
    }
  }
}
```

**Error Responses:**
- `400 Bad Request` - Missing or invalid question
- `500 Internal Server Error` - Processing error
- `503 Service Unavailable` - System overloaded

**Example Request:**
```bash
curl -X POST http://localhost:8080/api/ask \
  -H "Content-Type: application/json" \
  -d '{
    "question": "What is quantum computing?",
    "sessionId": "session-123"
  }'
```

### 2. Health Check
**Endpoint:** `GET /api/health`

**Description:** Check if the API server is running and healthy.

**Response:** `200 OK`
```json
{
  "status": "healthy",
  "timestamp": "2025-01-08T12:00:00.000Z",
  "version": "1.0.0"
}
```

**Example Request:**
```bash
curl http://localhost:8080/api/health
```

## Session Management

### How Sessions Work
- Sessions are automatically created on first request
- Session ID is returned in response
- Sessions expire after 1 hour (configurable via SESSION_TTL)
- Use the same sessionId for follow-up questions to maintain context

### Session Best Practices
1. Store sessionId in your application state
2. Send sessionId with each request in a conversation
3. Clear sessionId to start a new conversation
4. Handle expired sessions gracefully

## Response Times

Expected response times vary by query complexity:
- **Simple queries:** 5-10 seconds
- **Search-required queries:** 15-25 seconds
- **Complex multi-step queries:** 20-40 seconds
- **Maximum timeout:** 60 seconds (recommended)

## Rate Limiting

Currently no rate limiting implemented. Please be reasonable with request frequency.

## WebSocket Support (Future)

WebSocket support for streaming responses is planned but not yet implemented.

## Example Integration

### JavaScript/Fetch
```javascript
async function askQuestion(question, sessionId = null) {
  const response = await fetch('http://localhost:8080/api/ask', {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json'
    },
    body: JSON.stringify({
      question,
      sessionId
    })
  });
  
  if (!response.ok) {
    throw new Error(`API error: ${response.status}`);
  }
  
  return await response.json();
}

// Usage
const result = await askQuestion("What is AI?");
console.log(result.answer);
console.log("Session ID:", result.sessionId);

// Follow-up question
const followUp = await askQuestion(
  "How does it work?", 
  result.sessionId
);
```

### Python/Requests
```python
import requests

def ask_question(question, session_id=None):
    url = "http://localhost:8080/api/ask"
    payload = {"question": question}
    
    if session_id:
        payload["sessionId"] = session_id
    
    response = requests.post(url, json=payload)
    response.raise_for_status()
    return response.json()

# Usage
result = ask_question("What is machine learning?")
print(result["answer"])
print(f"Session ID: {result['sessionId']}")

# Follow-up
follow_up = ask_question(
    "What are the main types?", 
    result["sessionId"]
)
```

## Response Format Details

### Answer Format
- Responses are in **Markdown format**
- May include headers, lists, code blocks, and emphasis
- Citations appear as `[1]`, `[2]` inline with corresponding sources in metadata

### Expert Personas
The system automatically selects expert personas like:
- Dr. Evelyn Reed (various specializations)
- Dr. Marcus Chen (technical expert)
- Dr. Anya Sharma (business expert)

### Context Detection
The `contextUsed` flag indicates when the system uses previous conversation context. The `modernData` object provides detailed context analysis metrics.

## Error Handling

### Common Error Scenarios
1. **Missing Question**
   - Status: 400
   - Response: `{ "error": "Question is required" }`

2. **Server Error**
   - Status: 500
   - Response: `{ "error": "Internal server error", "details": "..." }`

3. **Timeout**
   - Implement client-side timeout of 60 seconds
   - Retry once if timeout occurs

## CORS

CORS is enabled for all origins. The server sends appropriate CORS headers:
```
Access-Control-Allow-Origin: *
Access-Control-Allow-Methods: GET, POST, OPTIONS
Access-Control-Allow-Headers: Content-Type
```

## Environment Variables

The backend uses these environment variables (frontend doesn't need to set these):
```
OPENROUTER_API_KEY=<key>
SERPER_API_KEY=<key>
USE_MODERN_CONTEXT=true
MAX_CONTEXT_TOKENS=2000
SESSION_TTL=3600000
```

## Testing the API

### Quick Test Script
```bash
# Health check
curl http://localhost:8080/api/health

# Ask a question
curl -X POST http://localhost:8080/api/ask \
  -H "Content-Type: application/json" \
  -d '{"question": "What is the capital of France?"}'

# Follow-up with session
curl -X POST http://localhost:8080/api/ask \
  -H "Content-Type: application/json" \
  -d '{
    "question": "What is its population?",
    "sessionId": "YOUR_SESSION_ID_HERE"
  }'
```

## Upcoming Features

- WebSocket support for streaming
- Batch question processing
- User authentication
- Rate limiting
- Response caching
- Custom expert selection

## Support

For API issues or feature requests, check the repository issues or README.md for more information.