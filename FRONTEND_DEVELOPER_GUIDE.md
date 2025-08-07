# Frontend Developer Integration Guide

## 🎯 Quick Start for Frontend Development

Welcome! This guide will help you understand the Ask Any Expert system to build a frontend interface. The system is a web-augmented AI expert that maintains conversation context and provides sourced, intelligent responses.

## 📋 Essential Files to Read (In Order)

### 1. **API Server - Your Main Integration Point**
**File:** `api-server.js`
- **Purpose:** REST API server that handles all client requests
- **Key Endpoints:**
  - `POST /api/ask` - Main endpoint for questions
  - `GET /api/health` - Health check endpoint
- **Important:** This is your primary integration point

### 2. **Configuration**
**File:** `config.js`
- **Purpose:** System configuration and API keys
- **What to know:** Understand available models and capabilities

**File:** `.env.example` (create from `.env`)
- **Purpose:** Environment variables needed
- **Important:** Frontend doesn't need API keys, backend handles them

### 3. **Core Workflow Engine**
**File:** `src/workflow-engine-modular.js`
- **Purpose:** Main orchestrator that processes queries
- **What to know:** Understand the multi-step processing flow

### 4. **Session Management**
**File:** `src/session-store.js`
- **Purpose:** Manages conversation sessions
- **Key concept:** Sessions maintain context across multiple queries
- **TTL:** Default 1 hour (3600000ms)

### 5. **Context System**
**File:** `src/modern-context-manager.js`
- **Purpose:** Advanced context understanding using AI techniques
- **What to know:** How context affects responses

## 🔌 API Integration Details

### Main Endpoint: `/api/ask`

**Request:**
```javascript
POST http://localhost:8080/api/ask
Content-Type: application/json

{
  "question": "What is machine learning?",
  "sessionId": "optional-session-id"  // Omit for new conversation
}
```

**Response:**
```javascript
{
  "answer": "Machine learning is...",
  "sessionId": "generated-session-id",  // Save this for follow-ups!
  "metadata": {
    "requestId": "abc123",
    "processingTime": 15234,
    "confidence": 0.92,
    "searchesPerformed": 2,
    "expertsConsulted": ["Dr. Evelyn Reed"],
    "contextUsed": true,
    "isFollowUp": false,
    "sources": [
      {
        "number": 1,
        "title": "Article Title",
        "url": "https://...",
        "snippet": "..."
      }
    ],
    "modernData": {
      "bayesianPosterior": 0.75,
      "contextType": "continuation",
      "informationScores": [0.8, 0.6, 0.4]
    }
  }
}
```

### Health Check Endpoint

```javascript
GET http://localhost:8080/api/health

Response:
{
  "status": "healthy",
  "timestamp": "2025-01-08T..."
}
```

## 💡 Frontend Implementation Recommendations

### 1. Session Management
```javascript
// Store sessionId in your frontend state/storage
let currentSessionId = null;

async function askQuestion(question) {
  const payload = { question };
  if (currentSessionId) {
    payload.sessionId = currentSessionId;
  }
  
  const response = await fetch('/api/ask', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(payload)
  });
  
  const data = await response.json();
  currentSessionId = data.sessionId; // Save for next query
  return data;
}
```

### 2. UI Features to Implement

**Essential:**
- Text input for questions
- Response display with markdown support
- Loading indicator (responses take 10-30 seconds)
- Session indicator (new conversation vs continuing)

**Recommended:**
- Source citations display with clickable links
- Expert persona indicator ("Speaking as Dr. X...")
- Confidence meter
- Context indicator (shows when using previous conversation)
- "New Conversation" button to clear session

**Advanced:**
- Response streaming (if you implement SSE)
- Search progress indicator
- Token usage tracking
- Response time display
- Export conversation feature

### 3. Handling Response Times

Responses typically take **10-30 seconds** due to:
- Web searching (when needed)
- Multi-step reasoning
- Context analysis

**Implement:**
- Clear loading states
- Optional: Progress messages ("Searching web...", "Analyzing context...")
- Timeout handling (60 seconds recommended)

### 4. Context Awareness Display

The system maintains context. Show users:
- When context is being used (metadata.contextUsed)
- Whether it's a follow-up (metadata.isFollowUp)
- Confidence level (metadata.confidence)
- Which expert is responding (metadata.expertsConsulted)

## 🏗️ System Architecture Overview

```
Frontend (Your Part)
    ↓
API Server (api-server.js)
    ↓
Workflow Engine (workflow-engine-modular.js)
    ↓
┌─────────────────────────────────────┐
│ • Context Manager (context analysis) │
│ • Query Analyzer (needs search?)     │
│ • Search Planner (generate queries)  │
│ • Web Search (Serper API)           │
│ • Info Synthesizer (final response)  │
└─────────────────────────────────────┘
```

## 📊 Important System Behaviors

1. **Automatic Expert Selection**: The system automatically chooses appropriate expert personas
2. **Smart Search**: Only searches when needed (not all queries require web search)
3. **Context Persistence**: Sessions maintain context for 1 hour
4. **Source Citations**: Web-sourced information includes numbered citations
5. **Confidence Scoring**: Each response has a confidence score (0-1)

## 🎨 UI/UX Best Practices

1. **Conversation Flow:**
   - Show clear conversation boundaries
   - Indicate when starting new conversation
   - Display context usage visually

2. **Response Display:**
   - Support markdown formatting
   - Highlight expert personas
   - Make citations interactive
   - Show confidence visually

3. **Error Handling:**
   - Network errors
   - Timeout (60s recommended)
   - Invalid session handling
   - Server unavailable

## 🚀 Quick Development Setup

1. **Start the backend:**
```bash
cd askanyexpert
npm install
npm start
# Server runs on http://localhost:8080
```

2. **Test the API:**
```bash
# Test health
curl http://localhost:8080/api/health

# Test a question
curl -X POST http://localhost:8080/api/ask \
  -H "Content-Type: application/json" \
  -d '{"question": "What is the speed of light?"}'
```

3. **CORS is enabled**, so you can develop frontend on any port

## 📝 Example Frontend Features

### Basic Chat Interface
```javascript
// Minimal implementation
const [messages, setMessages] = useState([]);
const [sessionId, setSessionId] = useState(null);
const [loading, setLoading] = useState(false);

const sendMessage = async (question) => {
  setLoading(true);
  setMessages(prev => [...prev, { role: 'user', content: question }]);
  
  try {
    const response = await fetch('http://localhost:8080/api/ask', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ question, sessionId })
    });
    
    const data = await response.json();
    setSessionId(data.sessionId);
    setMessages(prev => [...prev, { 
      role: 'assistant', 
      content: data.answer,
      metadata: data.metadata 
    }]);
  } catch (error) {
    console.error('Error:', error);
  } finally {
    setLoading(false);
  }
};
```

### Context Indicator
```javascript
// Show when context is being used
{message.metadata?.contextUsed && (
  <div className="context-indicator">
    🔗 Using conversation context
  </div>
)}
```

### Source Citations
```javascript
// Display sources
{message.metadata?.sources?.map(source => (
  <div key={source.number} className="citation">
    [{source.number}] 
    <a href={source.url} target="_blank">
      {source.title}
    </a>
  </div>
))}
```

## ⚠️ Important Notes

1. **API Keys**: Backend handles all API keys - frontend doesn't need them
2. **Rate Limiting**: Currently no rate limits, but be reasonable
3. **Session Cleanup**: Sessions auto-expire after 1 hour
4. **Response Format**: Responses are markdown - render accordingly
5. **Error Messages**: API returns appropriate HTTP status codes and error messages

## 🔍 Testing Your Frontend

Test these scenarios:
1. Single question (no context)
2. Follow-up question (with context)
3. Topic switch in conversation
4. New conversation after context
5. Long response handling
6. Source citation display
7. Network error recovery
8. Timeout handling (>60s)

## 📚 Additional Resources

- **README.md** - General project overview
- **test-realistic-conversations.js** - See example conversation flows
- **Development Phase Status** - Check README for current capabilities

## 💬 Response Examples

**Simple Query:**
- Question: "What is the speed of light?"
- Response Time: ~5-10 seconds
- Sources: Likely none (common knowledge)

**Complex Query:**
- Question: "Latest AI breakthroughs in 2025"
- Response Time: ~20-30 seconds
- Sources: Multiple web sources
- Contains citations [1], [2], etc.

**Follow-up Query:**
- Question: "How does it compare to last year?"
- Uses context from previous response
- metadata.contextUsed = true

## 🆘 Common Issues & Solutions

**Issue:** Responses take too long
**Solution:** Implement proper loading states, consider 60s timeout

**Issue:** Session not maintaining context
**Solution:** Ensure sessionId is saved and sent with each request

**Issue:** CORS errors
**Solution:** Backend has CORS enabled, check your frontend port

**Issue:** Markdown rendering issues
**Solution:** Use a markdown parser library

---

## Need More Info?

The backend developer can help with:
- API modifications
- Custom endpoints
- Response format changes
- Performance optimization
- Additional metadata

Good luck with the frontend! The system is designed to be frontend-agnostic, so you can use React, Vue, Angular, or vanilla JS.