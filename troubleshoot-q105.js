// Troubleshoot Q105 JSON parsing errors
const workflowEngine = require('./src/workflow-engine');
const fs = require('fs');

// Original Q105
const q105 = {
  question_id: 105,
  category: "chemistry",
  question: "Which of the following pairs of compounds are structural isomers?",
  options: [
    "A) CH₃CH₂OH and CH₃OCH₃",
    "B) CH₃CH₂CH₃ and CH₃CH₃",
    "C) CH₃COCH₃ and CH₃CH₂CHO",
    "D) CH₃CH₂Cl and CH₃CHCl₂",
    "E) CH₄ and C₂H₆",
    "F) H₂O and H₂O₂",
    "G) CO and CO₂",
    "H) NH₃ and N₂H₄",
    "I) HCl and HBr",
    "J) O₂ and O₃"
  ],
  correct_answer: "A"
};

async function testVariation(description, questionData, options = {}) {
  console.log('\n' + '='.repeat(80));
  console.log(`🧪 Testing: ${description}`);
  console.log('='.repeat(80));
  
  const query = `${questionData.question}\n\nOptions:\n${questionData.options.join('\n')}\n\nPlease select the correct answer and explain your reasoning. Make sure to clearly state which option (A, B, C, etc.) you are selecting as your final answer.`;
  
  console.log(`Query length: ${query.length} chars`);
  console.log(`Timeout: ${options.timeout || 600000}ms`);
  console.log(`Max tokens: ${options.max_tokens || 'default'}`);
  
  const startTime = Date.now();
  
  try {
    const response = await workflowEngine.answer(query, options);
    const duration = Date.now() - startTime;
    
    console.log(`✅ Success in ${Math.round(duration/1000)}s`);
    console.log(`Response length: ${response.content.length} chars`);
    
    // Check if response mentions the correct answer
    const mentions = ['A)', 'option A', 'answer is A', 'CH₃CH₂OH and CH₃OCH₃'].filter(term => 
      response.content.toLowerCase().includes(term.toLowerCase())
    );
    console.log(`Mentions of correct answer: ${mentions.join(', ') || 'none'}`);
    
    return { success: true, duration, responseLength: response.content.length };
    
  } catch (error) {
    const duration = Date.now() - startTime;
    console.log(`❌ Failed after ${Math.round(duration/1000)}s`);
    console.log(`Error: ${error.message}`);
    
    return { success: false, duration, error: error.message };
  }
}

async function runTroubleshooting() {
  console.log('🔍 TROUBLESHOOTING Q105 JSON PARSING ERRORS');
  console.log('Date:', new Date().toISOString());
  
  const results = [];
  
  // Test 1: Original format with default settings
  results.push(await testVariation('Original format (default settings)', q105));
  await new Promise(resolve => setTimeout(resolve, 10000));
  
  // Test 2: Shorter timeout
  results.push(await testVariation('Shorter timeout (2 minutes)', q105, { timeout: 120000 }));
  await new Promise(resolve => setTimeout(resolve, 10000));
  
  // Test 3: Explicit max_tokens
  results.push(await testVariation('With max_tokens=4000', q105, { max_tokens: 4000 }));
  await new Promise(resolve => setTimeout(resolve, 10000));
  
  // Test 4: Simplified question without subscripts
  const q105_simple = {
    ...q105,
    options: [
      "A) CH3CH2OH and CH3OCH3",
      "B) CH3CH2CH3 and CH3CH3",
      "C) CH3COCH3 and CH3CH2CHO",
      "D) CH3CH2Cl and CH3CHCl2",
      "E) CH4 and C2H6",
      "F) H2O and H2O2",
      "G) CO and CO2",
      "H) NH3 and N2H4",
      "I) HCl and HBr",
      "J) O2 and O3"
    ]
  };
  results.push(await testVariation('Without subscripts', q105_simple));
  await new Promise(resolve => setTimeout(resolve, 10000));
  
  // Test 5: Shorter question text
  const q105_short = {
    ...q105,
    question: "Which pairs are structural isomers?"
  };
  results.push(await testVariation('Shorter question text', q105_short));
  await new Promise(resolve => setTimeout(resolve, 10000));
  
  // Test 6: Fewer options
  const q105_fewer = {
    ...q105,
    options: q105.options.slice(0, 4) // Only A-D
  };
  results.push(await testVariation('Fewer options (A-D only)', q105_fewer));
  
  // Save results
  const troubleshootResults = {
    timestamp: new Date().toISOString(),
    question_id: 105,
    tests: results.map((r, i) => ({
      test: ['Original', 'Short timeout', 'Max tokens', 'No subscripts', 'Short question', 'Fewer options'][i],
      ...r
    }))
  };
  
  fs.writeFileSync('q105-troubleshoot-results.json', JSON.stringify(troubleshootResults, null, 2));
  
  // Summary
  console.log('\n' + '📊'.repeat(40));
  console.log('\n📈 TROUBLESHOOTING SUMMARY');
  console.log('─'.repeat(80));
  
  results.forEach((r, i) => {
    const testName = ['Original', 'Short timeout', 'Max tokens', 'No subscripts', 'Short question', 'Fewer options'][i];
    console.log(`${testName}: ${r.success ? '✅ Success' : '❌ Failed'} (${Math.round(r.duration/1000)}s)`);
  });
  
  const successCount = results.filter(r => r.success).length;
  console.log(`\nSuccess rate: ${successCount}/${results.length}`);
}

runTroubleshooting().catch(console.error);