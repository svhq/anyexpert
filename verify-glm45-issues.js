// Verify GLM-4.5 issues - check extraction and question quality
const fs = require('fs');
const workflowEngine = require('./src/workflow-engine');

// Load a few problematic questions to test
const testCases = [
  {
    id: 'orig-5',
    name: 'Psychology - Piaget',
    question: "According to Piaget's theory of cognitive development, a child who can think logically about concrete objects and events but struggles with abstract hypothetical thinking is most likely in which stage?",
    options: [
      'A) Sensorimotor stage',
      'B) Preoperational stage', 
      'C) Concrete operational stage',
      'D) Formal operational stage',
      'E) Postformal operational stage',
      'F) Symbolic stage',
      'G) Intuitive stage',
      'H) Schematic stage',
      'I) Abstract stage',
      'J) Logical stage'
    ],
    correct: 'C',
    glm45_selected: 'J'
  },
  {
    id: 'orig-8',
    name: 'CS - Merge Sort',
    question: 'In algorithm analysis, what is the time complexity of the merge sort algorithm in the worst case?',
    options: [
      'A) O(n)',
      'B) O(n log n)',
      'C) O(n²)',
      'D) O(log n)',
      'E) O(n³)',
      'F) O(2^n)',
      'G) O(n!)',
      'H) O(1)',
      'I) O(n^k)',
      'J) O(√n)'
    ],
    correct: 'B',
    glm45_selected: 'J'
  }
];

async function verifyQuestion(testCase) {
  console.log('\n' + '='.repeat(80));
  console.log(`\n🔍 Testing: ${testCase.name}`);
  console.log(`Previous result: Selected ${testCase.glm45_selected}, Correct: ${testCase.correct}`);
  
  // Construct query exactly as the test script does
  const query = `${testCase.question}\n\nOptions:\n${testCase.options.join('\n')}\n\nPlease select the correct answer and explain your reasoning.`;
  
  console.log('\n📤 Query being sent:');
  console.log(query);
  console.log('\n' + '-'.repeat(40));
  
  try {
    const response = await workflowEngine.answer(query, {});
    
    console.log('\n📤 Full Response:');
    console.log(response.content);
    console.log('\n' + '-'.repeat(40));
    
    // Try multiple extraction patterns
    const patterns = [
      // Standard patterns
      /(?:answer is|Answer:|correct answer is)[\s:]*([A-J])\)/i,
      /\*\*([A-J])\)\*\*/,
      /^([A-J])\)/m,
      /Option ([A-J])/i,
      // Additional patterns
      /select[\s:]*option[\s:]*([A-J])/i,
      /choose[\s:]*([A-J])\)/i,
      /\b([A-J])\)\s+is\s+(?:the\s+)?correct/i,
      /correct[\s:]+([A-J])\)/i
    ];
    
    let extractedAnswer = null;
    let matchedPattern = null;
    
    for (let i = 0; i < patterns.length; i++) {
      const match = response.content.match(patterns[i]);
      if (match) {
        extractedAnswer = match[1];
        matchedPattern = i;
        break;
      }
    }
    
    // Check if response mentions the correct answer at all
    const mentionsCorrect = response.content.includes(`${testCase.correct})`);
    
    // Check response quality
    const hasExplanation = response.content.length > 200;
    const hasOptions = testCase.options.some(opt => response.content.includes(opt.substring(3)));
    
    console.log('\n📊 ANALYSIS:');
    console.log(`- Response length: ${response.content.length} chars`);
    console.log(`- Has explanation: ${hasExplanation ? '✅' : '❌'}`);
    console.log(`- References options: ${hasOptions ? '✅' : '❌'}`);
    console.log(`- Mentions correct answer (${testCase.correct}): ${mentionsCorrect ? '✅' : '❌'}`);
    console.log(`- Extracted answer: ${extractedAnswer || 'NONE'}`);
    console.log(`- Matched pattern: ${matchedPattern !== null ? `Pattern ${matchedPattern}` : 'None'}`);
    console.log(`- Is correct: ${extractedAnswer === testCase.correct ? '✅' : '❌'}`);
    
    // Look for why it might select J
    if (extractedAnswer === 'J') {
      console.log('\n⚠️ Selected J - checking why:');
      const jIndex = response.content.lastIndexOf('J)');
      if (jIndex > -1) {
        const context = response.content.substring(Math.max(0, jIndex - 50), jIndex + 50);
        console.log('Context around J):', context);
      }
    }
    
    return {
      id: testCase.id,
      extracted: extractedAnswer,
      correct: testCase.correct,
      is_correct: extractedAnswer === testCase.correct,
      response_length: response.content.length,
      mentions_correct: mentionsCorrect
    };
    
  } catch (error) {
    console.error('\n❌ Error:', error.message);
    return { id: testCase.id, error: error.message };
  }
}

async function runVerification() {
  console.log('🧪 Verifying GLM-4.5 Issues');
  console.log('Model:', process.env.OPENROUTER_MODEL);
  
  const results = [];
  
  for (const testCase of testCases) {
    const result = await verifyQuestion(testCase);
    results.push(result);
    await new Promise(resolve => setTimeout(resolve, 3000));
  }
  
  console.log('\n' + '='.repeat(80));
  console.log('\n📊 SUMMARY:');
  results.forEach(r => {
    if (r.error) {
      console.log(`${r.id}: ERROR - ${r.error}`);
    } else {
      console.log(`${r.id}: Extracted ${r.extracted}, Correct ${r.correct} - ${r.is_correct ? '✅' : '❌'}`);
    }
  });
}

runVerification().catch(console.error);