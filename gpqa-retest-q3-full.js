const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const unifiedAgent = require('./src/unified-agent');
const e2bManager = require('./src/e2b-manager');
const fs = require('fs');

console.log('=== RETESTING GPQA Q3 WITH FULL RESPONSE ===\n');
console.log('Model:', process.env.OPENROUTER_MODEL);
console.log('Date:', new Date().toISOString());

// Verified Q3 data from CSV (note the typo "3 was treated" is in the original)
const q3 = {
  id: "gpqa-3",
  category: "Organic Chemistry",
  question: "trans-cinnamaldehyde was treated with methylmagnesium bromide, forming product 1. 1 was treated with pyridinium chlorochromate, forming product 2. 3 was treated with (dimethyl(oxo)-l6-sulfaneylidene)methane in DMSO at elevated temperature, forming product 3. how many carbon atoms are there in product 3?",
  options: {
    A: "11",
    B: "10", 
    C: "12",
    D: "14"
  },
  correct_answer: "A",
  explanation: "trans-cinnamaldehyde was treated with methylmagnesium bromide, forming (E)-4-phenylbut-3-en-2-ol. (E)-4-phenylbut-3-en-2-ol was treated with pyridinium chlorochromate, forming (E)-4-phenylbut-3-en-2-one. (E)-4-phenylbut-3-en-2-one was treated with (dimethyl(oxo)-l6-sulfaneylidene)methane in DMSO at elevated temperature, forming 1-(2-phenylcyclopropyl)ethan-1-one Chemical Formula: C11H12O."
};

console.log('VERIFIED DATA FROM SOURCE:');
console.log('- Question: ✓ Matches CSV exactly (including the "3" typo)');
console.log('- Options: ✓ All 4 options verified (11, 10, 12, 14)');
console.log('- Correct Answer: A (11) - confirmed from CSV');
console.log('- Note: CSV validators noted the typo "3 was treated" should be "2 was treated"\n');

// Format question
const formattedQuestion = `${q3.question}

Options:
A) ${q3.options.A}
B) ${q3.options.B}
C) ${q3.options.C}
D) ${q3.options.D}`;

console.log('Formatted Question:');
console.log(formattedQuestion);
console.log('\n' + '='.repeat(70) + '\n');

console.log('This is a hard graduate level organic chemistry question.');
console.log('Allowing extended timeout and full response...\n');

// Monitor E2B
let e2bExecutionCount = 0;
const originalExecute = e2bManager.executeWithFallback.bind(e2bManager);
e2bManager.executeWithFallback = async function(code, options) {
  e2bExecutionCount++;
  console.log(`\n=== E2B EXECUTION #${e2bExecutionCount} ===`);
  console.log('Code preview:', code.substring(0, 300) + '...');
  const result = await originalExecute(code, options);
  console.log('Execution completed');
  return result;
};

async function retestQuestion() {
  const startTime = Date.now();
  
  try {
    // Increase timeout for this complex question
    const result = await unifiedAgent.process(formattedQuestion, {}, {
      requestId: 'gpqa-q3-retest-' + Date.now(),
      timeout: 300000 // 5 minutes
    });
    
    const elapsedTime = ((Date.now() - startTime) / 1000).toFixed(1);
    
    // Save full response
    const filename = `gpqa-q3-retest-full-${Date.now()}.txt`;
    fs.writeFileSync(filename, result.content);
    
    console.log('\n=== RESPONSE COMPLETE ===');
    console.log('Response time:', elapsedTime, 'seconds');
    console.log('Response length:', result.content.length, 'characters');
    console.log('E2B executions:', e2bExecutionCount);
    console.log('Tools used:', result.metadata?.toolsUsed?.join(', ') || 'None');
    console.log('Full response saved to:', filename);
    
    // Show full response for manual checking
    console.log('\n' + '='.repeat(70));
    console.log('FULL RESPONSE:');
    console.log('='.repeat(70) + '\n');
    console.log(result.content);
    console.log('\n' + '='.repeat(70));
    
    // Manual extraction guidance
    console.log('\n=== MANUAL ANSWER EXTRACTION ===');
    console.log('The correct answer from CSV is: A (11 carbon atoms)');
    console.log('Please check the full response above for the model\'s final answer.');
    
    // Try automatic extraction
    const answerPatterns = [
      /final answer is[:\s]*\*?\*?\(?([A-D])\)?/i,
      /answer is[:\s]*\*?\*?\(?([A-D])\)?/i,
      /correct answer[:\s]*\*?\*?\(?([A-D])\)?/i,
      /Therefore[,\s]*\(?([A-D])\)?/i,
      /\*\*\(?([A-D])\)?\*\*/,
      /\boxed\{([A-D])\}/,
      /(\d+)\s*carbon atoms.*final answer/i,
      /final answer.*?(\d+)\s*carbon/i
    ];
    
    let extractedAnswer = null;
    for (const pattern of answerPatterns) {
      const match = result.content.match(pattern);
      if (match) {
        extractedAnswer = match[1];
        // If it's a number, map to letter
        if (extractedAnswer === "11") extractedAnswer = "A";
        else if (extractedAnswer === "10") extractedAnswer = "B";
        else if (extractedAnswer === "12") extractedAnswer = "C";
        else if (extractedAnswer === "14") extractedAnswer = "D";
        break;
      }
    }
    
    if (extractedAnswer) {
      console.log('\nAutomatically extracted answer:', extractedAnswer);
      console.log('Status:', extractedAnswer === q3.correct_answer ? 'CORRECT ✓' : 'INCORRECT ✗');
    }
    
    // Show explanation for reference
    console.log('\n=== CORRECT EXPLANATION (for reference) ===');
    console.log(q3.explanation);
    console.log('\nKey insight: This is a Johnson-Corey-Chaykovsky reaction that forms a cyclopropyl ketone.');
    
    return {
      question_id: q3.id,
      category: q3.category,
      correct_answer: q3.correct_answer,
      extracted_answer: extractedAnswer,
      response_time: elapsedTime,
      e2b_executions: e2bExecutionCount,
      response_file: filename,
      response_length: result.content.length,
      tools_used: result.metadata?.toolsUsed || []
    };
    
  } catch (error) {
    console.error('\n=== ERROR ===');
    console.error('Error:', error.message);
    if (error.message.includes('timeout')) {
      console.error('The question timed out again. This is a complex organic synthesis problem.');
    }
    return null;
  } finally {
    e2bManager.executeWithFallback = originalExecute;
    await e2bManager.shutdown();
  }
}

retestQuestion().then(result => {
  if (result) {
    console.log('\n=== FINAL SUMMARY ===');
    console.log(JSON.stringify(result, null, 2));
  }
}).catch(console.error);