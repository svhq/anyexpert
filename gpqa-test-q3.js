const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const unifiedAgent = require('./src/unified-agent');
const e2bManager = require('./src/e2b-manager');
const fs = require('fs');

// Load GPQA questions
const gpqaData = JSON.parse(fs.readFileSync('gpqa-10-questions-clean.json', 'utf-8'));
const q3 = gpqaData.questions[2];

console.log('=== TESTING GPQA Q3: ORGANIC CHEMISTRY ===\n');
console.log('Model:', process.env.OPENROUTER_MODEL);
console.log('Category:', q3.category);
console.log('Difficulty:', q3.metadata.difficulty);
console.log('Correct Answer:', q3.correct_answer, `(${q3.options[q3.correct_answer]})`);

// Format question with options
const formattedQuestion = `${q3.question}

Options:
A) ${q3.options.A}
B) ${q3.options.B}
C) ${q3.options.C}
D) ${q3.options.D}`;

console.log('\nFormatted Question:');
console.log(formattedQuestion);
console.log('\n' + '='.repeat(60) + '\n');
console.log('Processing...\n');

// Monitor E2B execution
let e2bExecutionCount = 0;
const originalExecuteWithFallback = e2bManager.executeWithFallback.bind(e2bManager);
e2bManager.executeWithFallback = async function(code, options) {
  e2bExecutionCount++;
  console.log(`\n=== E2B EXECUTION #${e2bExecutionCount} ===`);
  console.log('Code snippet:', code.substring(0, 200) + '...');
  const result = await originalExecuteWithFallback(code, options);
  console.log('Result:', result.results?.[0]?.text || 'No output');
  return result;
};

async function testQuestion() {
  const startTime = Date.now();
  
  try {
    const result = await unifiedAgent.process(formattedQuestion, {}, {
      requestId: 'gpqa-q3-' + Date.now()
    });
    
    const elapsedTime = ((Date.now() - startTime) / 1000).toFixed(1);
    
    // Save full response
    const filename = `gpqa-q3-response-${Date.now()}.txt`;
    fs.writeFileSync(filename, result.content);
    
    console.log('\n=== RESULTS ===');
    console.log('Response time:', elapsedTime, 'seconds');
    console.log('Response saved to:', filename);
    console.log('Response length:', result.content.length, 'characters');
    console.log('E2B executions:', e2bExecutionCount);
    console.log('Tools used:', result.metadata?.toolsUsed?.join(', ') || 'None');
    
    // Show the complete response for manual checking
    console.log('\n=== LAST 1000 CHARACTERS ===');
    const last1000 = result.content.substring(Math.max(0, result.content.length - 1000));
    console.log(last1000);
    
    // Extract answer manually
    console.log('\n=== MANUAL ANSWER CHECK ===');
    console.log('The correct answer is:', q3.correct_answer, `(${q3.options[q3.correct_answer]})`);
    console.log('Look for the model\'s final answer in the response above.');
    
    // Look for answer patterns
    const answerPatterns = [
      /final answer is[:\s]*\*?\*?\(?([A-D])\)?/i,
      /answer is[:\s]*\*?\*?\(?([A-D])\)?/i,
      /correct answer[:\s]*\*?\*?\(?([A-D])\)?/i,
      /Therefore[,\s]*\(?([A-D])\)?/i,
      /\*\*\(?([A-D])\)?\*\*/,
      /\boxed\{([A-D])\}/
    ];
    
    let extractedAnswer = null;
    for (const pattern of answerPatterns) {
      const match = result.content.match(pattern);
      if (match) {
        extractedAnswer = match[1].toUpperCase();
        break;
      }
    }
    
    if (extractedAnswer) {
      console.log('Extracted answer:', extractedAnswer);
      console.log('Status:', extractedAnswer === q3.correct_answer ? 'CORRECT ✓' : 'INCORRECT ✗');
    } else {
      console.log('Could not extract answer automatically - check manually');
    }
    
    return {
      question_id: q3.id,
      category: q3.category,
      correct_answer: q3.correct_answer,
      extracted_answer: extractedAnswer,
      response_time: elapsedTime,
      e2b_executions: e2bExecutionCount,
      response_file: filename,
      tools_used: result.metadata?.toolsUsed || [],
      match: extractedAnswer === q3.correct_answer
    };
    
  } catch (error) {
    console.error('\n=== ERROR ===');
    console.error('Error:', error.message);
    return null;
  } finally {
    await e2bManager.shutdown();
  }
}

testQuestion().then(result => {
  if (result) {
    console.log('\n=== SUMMARY ===');
    console.log(JSON.stringify(result, null, 2));
    
    // Show explanation for reference
    console.log('\n=== EXPLANATION (for reference) ===');
    console.log(q3.metadata.explanation);
  }
}).catch(console.error);