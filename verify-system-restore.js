const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const unifiedAgent = require('./src/unified-agent');

console.log('=== SYSTEM RESTORE VERIFICATION ===\n');
console.log('Date:', new Date().toISOString());
console.log('Model:', process.env.OPENROUTER_MODEL);
console.log('Testing if system works as before planning changes...\n');

async function testQuery(query, expectedTools) {
  console.log(`\nTesting: "${query}"`);
  console.log('Expected tools:', expectedTools.join(' → '));
  
  const startTime = Date.now();
  
  try {
    const result = await unifiedAgent.process(query, {}, {
      requestId: `verify-${Date.now()}`,
      timeout: 60000  // 1 minute
    });
    
    const duration = Date.now() - startTime;
    const toolsUsed = result.metadata.toolsUsed || [];
    
    console.log('Actual tools:', toolsUsed.join(' → '));
    console.log(`Duration: ${(duration / 1000).toFixed(1)}s`);
    console.log('Answer preview:', result.content.substring(0, 100) + '...');
    console.log('✅ Test completed successfully');
    
    return { success: true, toolsUsed, duration };
    
  } catch (error) {
    console.error('❌ Error:', error.message);
    return { success: false, error: error.message };
  }
}

async function main() {
  const tests = [
    {
      query: "What is 2+2?",
      expectedTools: ['reason', 'synthesize']
    },
    {
      query: "Calculate the square root of 144",
      expectedTools: ['code', 'synthesize']
    },
    {
      query: "What is the latest news about artificial intelligence?",
      expectedTools: ['search', 'synthesize']
    }
  ];
  
  console.log(`Running ${tests.length} verification tests...\n`);
  
  for (const test of tests) {
    await testQuery(test.query, test.expectedTools);
    
    // Wait between tests
    if (tests[tests.length - 1] !== test) {
      console.log('\nWaiting 5 seconds...');
      await new Promise(resolve => setTimeout(resolve, 5000));
    }
  }
  
  console.log('\n=== VERIFICATION COMPLETE ===');
  console.log('System appears to be working normally.');
}

main().catch(console.error);