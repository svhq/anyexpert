// Retest the 2 error questions from GLM-4.5-air test
const workflowEngine = require('./src/workflow-engine');
const fs = require('fs');

// Load the questions
const questionsData = JSON.parse(fs.readFileSync('mmlu-30-questions-clean.json', 'utf8'));
const allQuestions = questionsData.questions;

// Find the error questions
const q105 = allQuestions.find(q => q.question_id === 105);
const q1986 = allQuestions.find(q => q.question_id === 1986);

const errorQuestions = [q105, q1986];

async function testQuestion(questionData) {
  console.log('\n' + '='.repeat(100));
  console.log(`📝 Testing Question ID ${questionData.question_id} (${questionData.category})`);
  console.log('='.repeat(100));
  
  console.log(`\n📋 QUESTION:`);
  console.log(questionData.question);
  
  console.log(`\n📋 OPTIONS:`);
  questionData.options.forEach(opt => console.log(`  ${opt}`));
  
  console.log(`\n📌 Expected Answer: ${questionData.correct_answer}`);
  
  const query = `${questionData.question}\n\nOptions:\n${questionData.options.join('\n')}\n\nPlease select the correct answer and explain your reasoning. Make sure to clearly state which option (A, B, C, etc.) you are selecting as your final answer.`;
  
  console.log(`\n⏳ Sending to GLM-4.5-air...`);
  const startTime = Date.now();
  
  try {
    const response = await workflowEngine.answer(query, { timeout: 600000 }); // 10 minutes timeout
    const duration = Date.now() - startTime;
    
    console.log(`✅ Response received in ${Math.round(duration/1000)}s (${response.content.length} chars)`);
    
    // Show full response
    console.log(`\n📄 FULL RESPONSE:`);
    console.log('─'.repeat(100));
    console.log(response.content);
    console.log('─'.repeat(100));
    
    return {
      question_id: questionData.question_id,
      category: questionData.category,
      question: questionData.question,
      options: questionData.options,
      expected_answer: questionData.correct_answer,
      model_response: response.content,
      response_length: response.content.length,
      response_time_ms: duration,
      timestamp: new Date().toISOString(),
      success: true
    };
    
  } catch (error) {
    console.log(`❌ ERROR: ${error.message}`);
    
    return {
      question_id: questionData.question_id,
      category: questionData.category,
      expected_answer: questionData.correct_answer,
      error: error.message,
      timestamp: new Date().toISOString(),
      success: false
    };
  }
}

async function runRetest() {
  console.log('\n🔄 RETESTING ERROR QUESTIONS WITH GLM-4.5-AIR');
  console.log('Model:', process.env.OPENROUTER_MODEL || 'z-ai/glm-4.5-air:free');
  console.log('Date:', new Date().toISOString());
  console.log('Questions to retest:', errorQuestions.map(q => `ID ${q.question_id}`).join(', '));
  
  const results = [];
  
  for (const question of errorQuestions) {
    const result = await testQuestion(question);
    results.push(result);
    
    if (errorQuestions.indexOf(question) < errorQuestions.length - 1) {
      console.log(`\n⏳ Waiting 10 seconds before next question...`);
      await new Promise(resolve => setTimeout(resolve, 10000));
    }
  }
  
  // Save results
  const retestResults = {
    test_name: 'GLM-4.5-air Error Questions Retest',
    model: process.env.OPENROUTER_MODEL || 'z-ai/glm-4.5-air:free',
    timestamp: new Date().toISOString(),
    total_questions: results.length,
    results: results
  };
  
  fs.writeFileSync('glm45air-error-retest.json', JSON.stringify(retestResults, null, 2));
  console.log('\n💾 Retest results saved to: glm45air-error-retest.json');
  
  // Show summary
  console.log('\n' + '📊'.repeat(50));
  console.log('\n📈 RETEST SUMMARY');
  console.log('─'.repeat(100));
  
  for (const result of results) {
    if (result.success) {
      console.log(`Question ${result.question_id}: Response length ${result.response_length} chars`);
    } else {
      console.log(`Question ${result.question_id}: FAILED - ${result.error}`);
    }
  }
}

runRetest().catch(console.error);