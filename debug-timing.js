const fetch = require('node-fetch');

async function debugTiming() {
  console.log('=== TIMING ANALYSIS ===\n');
  
  // Test 1: Direct LLM call speed
  console.log('Test 1: Direct LLM generation speed...');
  const startLLM = Date.now();
  
  const llmResponse = await fetch('https://openrouter.ai/api/v1/chat/completions', {
    method: 'POST',
    headers: {
      'Authorization': `Bearer ${process.env.OPENROUTER_API_KEY}`,
      'Content-Type': 'application/json',
    },
    body: JSON.stringify({
      model: 'openai/o1-mini',
      messages: [
        {
          role: 'user',
          content: 'Generate a very long Python script that analyzes functions. Make it at least 1000 lines with lots of comments.'
        }
      ],
      temperature: 0.1,
      max_tokens: 50000
    })
  });
  
  const llmTime = Date.now() - startLLM;
  const llmData = await llmResponse.json();
  const tokenCount = llmData.usage?.total_tokens || 0;
  const responseLength = llmData.choices?.[0]?.message?.content?.length || 0;
  
  console.log(`- Time: ${(llmTime/1000).toFixed(1)}s`);
  console.log(`- Tokens: ${tokenCount}`);
  console.log(`- Response length: ${responseLength} chars`);
  console.log(`- Speed: ${(tokenCount / (llmTime/1000)).toFixed(0)} tokens/sec\n`);
  
  // Test 2: E2B execution speed
  console.log('Test 2: E2B execution speed...');
  const simpleCode = `
import time
print("Start")
time.sleep(1)
print("End")
result = sum(range(1000000))
print(f"Sum: {result}")
`;
  
  const startE2B = Date.now();
  const e2bResponse = await fetch('http://localhost:3001/run_code', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      language: 'python',
      source: simpleCode,
      timeout: 10000
    })
  });
  
  const e2bTime = Date.now() - startE2B;
  const e2bResult = await e2bResponse.json();
  
  console.log(`- Time: ${(e2bTime/1000).toFixed(1)}s`);
  console.log(`- Exit code: ${e2bResult.exitCode}`);
  console.log(`- Output: ${e2bResult.stdout?.substring(0, 100)}...\n`);
  
  // Test 3: Large code execution
  console.log('Test 3: Large code execution...');
  const largeCode = `
# Large Python script
import numpy as np
import time

print("Starting large computation...")
start = time.time()

# Generate lots of output
for i in range(100):
    print(f"Iteration {i}: " + "x" * 100)
    
# Do some actual computation
matrix = np.random.rand(100, 100)
result = np.linalg.det(matrix)
print(f"\\nDeterminant: {result}")

end = time.time()
print(f"Total time: {end - start:.2f}s")
`;
  
  const startLargeE2B = Date.now();
  const largeE2bResponse = await fetch('http://localhost:3001/run_code', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      language: 'python',
      source: largeCode,
      timeout: 30000
    })
  });
  
  const largeE2bTime = Date.now() - startLargeE2B;
  const largeE2bResult = await largeE2bResponse.json();
  
  console.log(`- Time: ${(largeE2bTime/1000).toFixed(1)}s`);
  console.log(`- Exit code: ${largeE2bResult.exitCode}`);
  console.log(`- Output length: ${largeE2bResult.stdout?.length || 0} chars\n`);
  
  console.log('=== ANALYSIS ===');
  console.log('If 34K tokens took 75s total, and 140K tokens took 257s total:');
  console.log(`- 34K tokens generation time: ~${(34000 / (tokenCount / (llmTime/1000))).toFixed(0)}s`);
  console.log(`- 140K tokens generation time: ~${(140000 / (tokenCount / (llmTime/1000))).toFixed(0)}s`);
  console.log('\nThe extra 182s (257s - 75s) must be coming from somewhere else...');
}

require('dotenv').config();
debugTiming().catch(console.error);