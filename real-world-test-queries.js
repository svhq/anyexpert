const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const unifiedAgent = require('./src/unified-agent');
const e2bManager = require('./src/e2b-manager');
const fs = require('fs');

// 5 diverse real-world queries that users might ask experts
const queries = [
  {
    id: 1,
    category: "Health/Medical",
    query: "I've been experiencing persistent fatigue for 3 months despite sleeping 8 hours nightly. I'm 35, exercise regularly, and eat well. My doctor ran basic blood tests that came back normal. What other causes should I investigate, and what specific tests should I request?"
  },
  {
    id: 2,
    category: "Technology/Programming",
    query: "I'm building a React app that needs to handle real-time updates for 10,000+ concurrent users. Should I use WebSockets, Server-Sent Events, or long polling? What are the trade-offs and which cloud infrastructure would scale best?"
  },
  {
    id: 3,
    category: "Finance/Investment",
    query: "I'm 28 with $50k saved and want to start investing. My goal is retirement by 55. Should I prioritize maxing out 401k, opening a Roth IRA, or investing in index funds? How should I balance these with current inflation rates?"
  },
  {
    id: 4,
    category: "Science/Environment",
    query: "My garden soil pH is 7.8 and I want to grow blueberries which need acidic soil. What's the safest way to lower pH without harming existing plants nearby? How long will it take to see results?"
  },
  {
    id: 5,
    category: "Legal/Business",
    query: "I'm a freelance designer wanting to protect my work. Should I register copyrights for each project, use watermarks, or rely on automatic copyright? What's the most cost-effective approach for someone doing 20+ projects yearly?"
  }
];

// Monitor E2B
let globalE2bCount = 0;
const originalExecute = e2bManager.executeWithFallback.bind(e2bManager);
e2bManager.executeWithFallback = async function(code, options) {
  globalE2bCount++;
  console.log(`\\n[E2B Execution #${globalE2bCount}]`);
  const result = await originalExecute(code, options);
  return result;
};

async function testQuery(query, model) {
  console.log(`\\n${'='.repeat(80)}`);
  console.log(`Testing Query ${query.id} (${query.category}) with ${model}`);
  console.log('='.repeat(80));
  console.log(`\\nQuery: "${query.query.substring(0, 100)}..."`);
  
  const startTime = Date.now();
  const startE2b = globalE2bCount;
  
  try {
    const result = await unifiedAgent.process(query.query, {}, {
      requestId: `real-q${query.id}-${model}-${Date.now()}`
    });
    
    const elapsed = ((Date.now() - startTime) / 1000).toFixed(1);
    const e2bUsed = globalE2bCount - startE2b;
    
    // Save response
    const filename = `real-world-q${query.id}-${model.replace('/', '-')}-${Date.now()}.txt`;
    fs.writeFileSync(filename, result.content);
    
    console.log(`\\nResponse time: ${elapsed}s`);
    console.log(`E2B executions: ${e2bUsed}`);
    console.log(`Tools used: ${result.metadata?.toolsUsed?.join(', ') || 'None'}`);
    console.log(`Response length: ${result.content.length} chars`);
    console.log(`Saved to: ${filename}`);
    
    // Show preview
    console.log('\\n--- RESPONSE PREVIEW (first 500 chars) ---');
    console.log(result.content.substring(0, 500) + '...');
    
    return {
      query_id: query.id,
      category: query.category,
      model,
      response_time: elapsed,
      e2b_count: e2bUsed,
      tools_used: result.metadata?.toolsUsed || [],
      response_length: result.content.length,
      filename
    };
    
  } catch (error) {
    console.error('ERROR:', error.message);
    return {
      query_id: query.id,
      category: query.category,
      model,
      error: error.message
    };
  }
}

async function runComparison() {
  const results = {
    'x-ai/grok-3-mini': [],
    'google/gemini-2.0-flash-lite': []
  };
  
  // Test with grok-3-mini first
  console.log('\\n' + '='.repeat(80));
  console.log('PHASE 1: TESTING WITH GROK-3-MINI');
  console.log('='.repeat(80));
  
  process.env.OPENROUTER_MODEL = 'x-ai/grok-3-mini';
  for (const query of queries) {
    const result = await testQuery(query, 'x-ai/grok-3-mini');
    results['x-ai/grok-3-mini'].push(result);
    await new Promise(resolve => setTimeout(resolve, 2000));
  }
  
  // Switch to gemini and test
  console.log('\\n\\n' + '='.repeat(80));
  console.log('PHASE 2: TESTING WITH GEMINI-2.0-FLASH-LITE');
  console.log('='.repeat(80));
  
  process.env.OPENROUTER_MODEL = 'google/gemini-2.0-flash-lite';
  for (const query of queries) {
    const result = await testQuery(query, 'google/gemini-2.0-flash-lite');
    results['google/gemini-2.0-flash-lite'].push(result);
    await new Promise(resolve => setTimeout(resolve, 2000));
  }
  
  // Save results
  const summaryFile = `real-world-comparison-${Date.now()}.json`;
  fs.writeFileSync(summaryFile, JSON.stringify(results, null, 2));
  
  // Print summary
  console.log('\\n\\n' + '='.repeat(80));
  console.log('SUMMARY');
  console.log('='.repeat(80));
  
  console.log('\\nGrok-3-mini:');
  const grokAvgTime = results['x-ai/grok-3-mini'].reduce((sum, r) => sum + parseFloat(r.response_time || 0), 0) / 5;
  const grokAvgLength = results['x-ai/grok-3-mini'].reduce((sum, r) => sum + (r.response_length || 0), 0) / 5;
  console.log(`- Average response time: ${grokAvgTime.toFixed(1)}s`);
  console.log(`- Average response length: ${Math.round(grokAvgLength)} chars`);
  console.log(`- Total E2B executions: ${results['x-ai/grok-3-mini'].reduce((sum, r) => sum + (r.e2b_count || 0), 0)}`);
  
  console.log('\\nGemini-2.0-flash-lite:');
  const geminiAvgTime = results['google/gemini-2.0-flash-lite'].reduce((sum, r) => sum + parseFloat(r.response_time || 0), 0) / 5;
  const geminiAvgLength = results['google/gemini-2.0-flash-lite'].reduce((sum, r) => sum + (r.response_length || 0), 0) / 5;
  console.log(`- Average response time: ${geminiAvgTime.toFixed(1)}s`);
  console.log(`- Average response length: ${Math.round(geminiAvgLength)} chars`);
  console.log(`- Total E2B executions: ${results['google/gemini-2.0-flash-lite'].reduce((sum, r) => sum + (r.e2b_count || 0), 0)}`);
  
  console.log(`\\nResults saved to: ${summaryFile}`);
  console.log('\\nNOTE: Please manually review the full responses to assess quality!');
  
  await e2bManager.shutdown();
}

// Run the comparison
runComparison().catch(console.error);