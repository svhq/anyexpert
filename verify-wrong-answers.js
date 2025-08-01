// Verify all wrong answers from both models
const fs = require('fs');
const workflowEngine = require('./src/workflow-engine');

// Load test results and questions
const questions = JSON.parse(fs.readFileSync('./mmlu-10-new-questions.json', 'utf8')).questions;
const o4Results = JSON.parse(fs.readFileSync('./o4-mini-high-10-results-incremental.json', 'utf8'));
const glmResults = JSON.parse(fs.readFileSync('./glm-4.5-air-10-results-incremental.json', 'utf8'));

// Find all wrong answers
const wrongAnswers = [];

// Check o4-mini-high wrong answers
o4Results.results.forEach((result, idx) => {
  if (!result.is_correct) {
    wrongAnswers.push({
      model: 'o4-mini-high',
      question_num: result.question_number,
      question: questions[idx],
      result: result
    });
  }
});

// Check GLM wrong answers
glmResults.results.forEach((result, idx) => {
  if (!result.is_correct) {
    wrongAnswers.push({
      model: 'GLM-4.5-air',
      question_num: result.question_number,
      question: questions[idx],
      result: result
    });
  }
});

console.log(`\n🔍 Found ${wrongAnswers.length} wrong answers to verify\n`);

async function verifyWrongAnswer(item) {
  console.log('='.repeat(80));
  console.log(`\n🤖 Model: ${item.model}`);
  console.log(`📚 Question ${item.question_num}: ${item.question.category.toUpperCase()}`);
  console.log(`\n❓ Question: ${item.question.question}`);
  console.log(`\n✅ Correct Answer: ${item.question.correct_answer}`);
  console.log(`❌ Model Selected: ${item.result.selected_answer || 'NULL'}`);
  
  // Temporarily switch to the model that gave the wrong answer
  const originalModel = process.env.OPENROUTER_MODEL;
  process.env.OPENROUTER_MODEL = item.model === 'o4-mini-high' ? 'openai/o4-mini-high' : 'z-ai/glm-4.5-air:free';
  
  console.log(`\n🔄 Re-testing with ${item.model}...`);
  
  const query = `${item.question.question}\n\nOptions:\n${item.question.options.join('\n')}\n\nPlease select the correct answer and explain your reasoning.`;
  
  try {
    const response = await workflowEngine.answer(query, {});
    
    console.log('\n📤 FULL RESPONSE:');
    console.log(response.content);
    console.log('\n' + '-'.repeat(40));
    
    // Comprehensive answer search
    const correctLetter = item.question.correct_answer;
    const correctOption = item.question.options.find(opt => opt.startsWith(correctLetter + ')'));
    
    // Check various ways the answer might appear
    const checks = {
      'Exact match': response.content.includes(`${correctLetter})`),
      'Answer statement': response.content.match(new RegExp(`answer is.*${correctLetter}`, 'i')),
      'Correct option text': correctOption ? response.content.includes(correctOption.substring(3)) : false,
      'Bold format': response.content.includes(`**${correctLetter})**`),
      'Final answer': response.content.match(new RegExp(`final answer.*${correctLetter}`, 'i'))
    };
    
    console.log('\n🔎 VERIFICATION CHECKS:');
    Object.entries(checks).forEach(([check, found]) => {
      console.log(`- ${check}: ${found ? '✅ YES' : '❌ NO'}`);
    });
    
    // Try to extract what the model actually selected
    const patterns = [
      /(?:answer is|Answer:|correct answer is|select|choose)[\s:]*([A-J])\)/i,
      /\*\*([A-J])\)\*\*/,
      /^([A-J])\)/m,
      /Option ([A-J])/i
    ];
    
    let extractedAnswer = null;
    for (const pattern of patterns) {
      const match = response.content.match(pattern);
      if (match) {
        extractedAnswer = match[1];
        break;
      }
    }
    
    console.log(`\n📊 ANALYSIS:`);
    console.log(`- Extracted Answer: ${extractedAnswer || 'Could not extract'}`);
    console.log(`- Contains Correct Answer: ${Object.values(checks).some(v => v) ? '✅ YES' : '❌ NO'}`);
    console.log(`- Original Extraction: ${item.result.selected_answer || 'NULL'}`);
    console.log(`- Verdict: ${extractedAnswer === item.question.correct_answer ? '🟢 Actually CORRECT!' : '🔴 Truly WRONG'}`);
    
    return {
      model: item.model,
      question_num: item.question_num,
      correct_answer: item.question.correct_answer,
      original_extraction: item.result.selected_answer,
      new_extraction: extractedAnswer,
      actually_correct: extractedAnswer === item.question.correct_answer,
      response_length: response.content.length
    };
    
  } catch (error) {
    console.error('\n❌ Error:', error.message);
    return {
      model: item.model,
      question_num: item.question_num,
      error: error.message
    };
  } finally {
    // Restore original model
    process.env.OPENROUTER_MODEL = originalModel;
  }
}

async function runVerification() {
  const results = [];
  
  for (const wrongAnswer of wrongAnswers) {
    const result = await verifyWrongAnswer(wrongAnswer);
    results.push(result);
    
    // Wait between tests
    await new Promise(resolve => setTimeout(resolve, 3000));
  }
  
  // Summary
  console.log('\n' + '='.repeat(80));
  console.log('\n📊 VERIFICATION SUMMARY\n');
  
  const actuallyCorrect = results.filter(r => r.actually_correct);
  console.log(`🟢 Actually Correct (extraction error): ${actuallyCorrect.length}`);
  console.log(`🔴 Truly Wrong: ${results.length - actuallyCorrect.length}`);
  
  console.log('\n📋 Detailed Results:');
  results.forEach(r => {
    if (r.error) {
      console.log(`\n${r.model} Q${r.question_num}: ❌ ERROR - ${r.error}`);
    } else {
      console.log(`\n${r.model} Q${r.question_num}:`);
      console.log(`  Correct: ${r.correct_answer}, Original: ${r.original_extraction}, New: ${r.new_extraction}`);
      console.log(`  Status: ${r.actually_correct ? '🟢 EXTRACTION ERROR' : '🔴 TRULY WRONG'}`);
    }
  });
  
  // Save verification results
  fs.writeFileSync('wrong-answers-verification.json', JSON.stringify({
    timestamp: new Date().toISOString(),
    total_wrong_answers: wrongAnswers.length,
    extraction_errors: actuallyCorrect.length,
    truly_wrong: results.length - actuallyCorrect.length,
    details: results
  }, null, 2));
  
  console.log('\n💾 Saved verification to wrong-answers-verification.json');
}

runVerification().catch(console.error);