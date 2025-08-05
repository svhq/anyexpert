const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });

console.log('=== E2B API STATUS CHECK ===\n');

async function checkE2BStatus() {
  const apiKey = process.env.E2B_API_KEY;
  
  if (!apiKey) {
    console.error('❌ E2B_API_KEY not found in environment');
    return;
  }
  
  console.log('API Key:', apiKey.substring(0, 10) + '...' + apiKey.substring(apiKey.length - 4));
  console.log('\nTesting E2B API directly...\n');
  
  try {
    // Try to list sandboxes
    const response = await fetch('https://api.e2b.dev/sandboxes', {
      method: 'GET',
      headers: {
        'X-API-Key': apiKey,
        'Content-Type': 'application/json'
      }
    });
    
    console.log('API Response Status:', response.status, response.statusText);
    
    if (response.ok) {
      const data = await response.json();
      console.log('✅ API is accessible');
      console.log('Active sandboxes:', Array.isArray(data) ? data.length : 'Unknown');
    } else {
      const error = await response.text();
      console.log('❌ API Error:', error);
      
      if (response.status === 429) {
        console.log('\n⚠️  RATE LIMIT DETECTED');
        
        // Check rate limit headers
        const retryAfter = response.headers.get('Retry-After');
        const rateLimitRemaining = response.headers.get('X-RateLimit-Remaining');
        const rateLimitReset = response.headers.get('X-RateLimit-Reset');
        
        if (retryAfter) console.log('Retry after:', retryAfter, 'seconds');
        if (rateLimitRemaining) console.log('Rate limit remaining:', rateLimitRemaining);
        if (rateLimitReset) {
          const resetDate = new Date(parseInt(rateLimitReset) * 1000);
          console.log('Rate limit resets at:', resetDate.toISOString());
        }
      }
    }
    
  } catch (error) {
    console.error('❌ Network Error:', error.message);
  }
}

checkE2BStatus().catch(console.error);