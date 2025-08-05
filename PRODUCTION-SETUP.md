# Production Setup Guide

## E2B Rate Limit Issue

You're seeing "Rate limit exceeded" errors because you're likely on E2B's free tier, which has strict limits on API calls.

### Current Behavior
- System tries E2B first
- If rate limited, automatically falls back to local Python
- **Everything still works**, just with a slight delay

### Production Options

#### Option 1: Upgrade E2B (Recommended)
1. Visit https://e2b.dev/pricing
2. Choose a paid plan based on your expected usage
3. Update your E2B_API_KEY in .env

#### Option 2: Disable E2B, Use Local Python Only
If you want to avoid E2B entirely, modify `src/e2b-manager.js`:

```javascript
// In executeWithFallback method, skip E2B and go straight to local:
async executeWithFallback(code, options = {}) {
  // Skip E2B, use local Python directly
  return this.executeLocalPython(code, options);
}
```

#### Option 3: Keep Current Setup
- The fallback mechanism works well
- Users won't see errors (handled internally)
- Just a bit slower when rate limited

### Monitoring
Add environment variable to track:
```
E2B_FALLBACK_COUNT=0  # Track how often fallback is used
```

### Recommendation
For production with paying customers, Option 1 (upgrade E2B) is best for:
- Reliability
- Performance  
- Sandboxed security
- No local Python dependency