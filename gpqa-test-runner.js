const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const unifiedAgent = require('./src/unified-agent');
const e2bManager = require('./src/e2b-manager');
const fs = require('fs');

// Load GPQA questions
const gpqaData = JSON.parse(fs.readFileSync('gpqa-10-questions-clean.json', 'utf-8'));

async function testQuestion(question, index) {
  console.log(`\n${'='.repeat(70)}`);
  console.log(`TESTING GPQA Q${index + 1}: ${question.category.toUpperCase()}`);
  console.log(`${'='.repeat(70)}\n`);
  
  console.log('Model:', process.env.OPENROUTER_MODEL);
  console.log('Difficulty:', question.metadata.difficulty);
  console.log('Correct Answer:', question.correct_answer, `(${question.options[question.correct_answer]})`);
  
  // Format question with options
  const formattedQuestion = `${question.question}

Options:
A) ${question.options.A}
B) ${question.options.B}
C) ${question.options.C}
D) ${question.options.D}`;
  
  console.log('\nProcessing...\n');
  
  // Monitor E2B execution
  let e2bExecutionCount = 0;
  const originalExecute = e2bManager.executeWithFallback.bind(e2bManager);
  e2bManager.executeWithFallback = async function(code, options) {
    e2bExecutionCount++;
    console.log(`E2B Execution #${e2bExecutionCount}`);
    return await originalExecute(code, options);
  };
  
  const startTime = Date.now();
  
  try {
    const result = await unifiedAgent.process(formattedQuestion, {}, {
      requestId: `gpqa-q${index + 1}-${Date.now()}`
    });
    
    const elapsedTime = ((Date.now() - startTime) / 1000).toFixed(1);
    
    // Save response
    const filename = `gpqa-q${index + 1}-response-${Date.now()}.txt`;
    fs.writeFileSync(filename, result.content);
    
    console.log('\nResponse time:', elapsedTime, 's');
    console.log('E2B executions:', e2bExecutionCount);
    console.log('Tools used:', result.metadata?.toolsUsed?.join(', ') || 'None');
    console.log('Response saved to:', filename);
    
    // Show last part of response
    console.log('\n--- LAST 800 CHARACTERS ---');
    const last800 = result.content.substring(Math.max(0, result.content.length - 800));
    console.log(last800);
    
    // Extract answer
    const answerPatterns = [
      /final answer is[:\s]*\*?\*?\(?([A-D])\)?/i,
      /answer is[:\s]*\*?\*?\(?([A-D])\)?/i,
      /correct answer[:\s]*\*?\*?\(?([A-D])\)?/i,
      /\\boxed\{([A-D])\}/,
      /\*\*\(?([A-D])\)?\*\*/,
      /Therefore[,\s]*\(?([A-D])\)?/i
    ];
    
    let extractedAnswer = null;
    for (const pattern of answerPatterns) {
      const match = result.content.match(pattern);
      if (match) {
        extractedAnswer = match[1].toUpperCase();
        break;
      }
    }
    
    console.log('\\n--- ANSWER CHECK ---');
    console.log('Correct answer:', question.correct_answer, `(${question.options[question.correct_answer]})`);
    console.log('Extracted answer:', extractedAnswer || 'MANUAL CHECK NEEDED');
    console.log('Status:', extractedAnswer === question.correct_answer ? 'CORRECT ✓' : 'INCORRECT ✗');
    
    // Restore original execute function
    e2bManager.executeWithFallback = originalExecute;
    
    return {
      question_id: question.id,
      category: question.category,
      correct_answer: question.correct_answer,
      extracted_answer: extractedAnswer,
      response_time: elapsedTime,
      e2b_executions: e2bExecutionCount,
      filename: filename,
      match: extractedAnswer === question.correct_answer
    };
    
  } catch (error) {
    console.error('ERROR:', error.message);
    e2bManager.executeWithFallback = originalExecute;
    return {
      question_id: question.id,
      category: question.category,
      error: error.message
    };
  }
}

async function runTests(startIndex = 4, endIndex = 10) {
  console.log(`\\n=== TESTING GPQA QUESTIONS ${startIndex + 1}-${endIndex} ===\\n`);
  
  const results = [];
  
  for (let i = startIndex; i < endIndex && i < gpqaData.questions.length; i++) {
    const result = await testQuestion(gpqaData.questions[i], i);
    results.push(result);
    
    // Small delay between questions
    await new Promise(resolve => setTimeout(resolve, 2000));
  }
  
  // Summary
  console.log('\\n\\n' + '='.repeat(70));
  console.log('SUMMARY');
  console.log('='.repeat(70) + '\\n');
  
  const correct = results.filter(r => r.match).length;
  const errors = results.filter(r => r.error).length;
  
  console.log(`Total questions: ${results.length}`);
  console.log(`Correct: ${correct}`);
  console.log(`Incorrect: ${results.length - correct - errors}`);
  console.log(`Errors: ${errors}`);
  console.log(`\\nDetailed results:`);
  
  results.forEach(r => {
    if (r.error) {
      console.log(`${r.question_id}: ERROR - ${r.error}`);
    } else {
      const status = r.match ? '✓' : '✗';
      console.log(`${r.question_id}: ${status} Model: ${r.extracted_answer || '?'}, Correct: ${r.correct_answer} (${r.response_time}s, E2B: ${r.e2b_executions})`);
    }
  });
  
  // Save results
  const summaryFile = `gpqa-q${startIndex + 1}-q${endIndex}-summary.json`;
  fs.writeFileSync(summaryFile, JSON.stringify(results, null, 2));
  console.log(`\\nResults saved to: ${summaryFile}`);
  
  await e2bManager.shutdown();
}

// Run tests for Q5-Q10
runTests(4, 10).catch(console.error);