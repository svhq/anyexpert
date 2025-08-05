const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const unifiedAgent = require('./src/unified-agent');
const fs = require('fs');

console.log('=== RETESTING GPQA Q11 & Q13 WITH MANUAL VERIFICATION ===\n');
console.log('Date:', new Date().toISOString());
console.log('Model:', process.env.OPENROUTER_MODEL);

// Load the questions
const data = JSON.parse(fs.readFileSync(path.join(__dirname, 'gpqa-q11-q20-randomized.json'), 'utf8'));

// Get Q11 and Q13
const q11 = data.questions[0];  // Q11 is at index 0
const q13 = data.questions[2];  // Q13 is at index 2

async function testQuestion(question, questionNum) {
  console.log(`\n${'='.repeat(80)}`);
  console.log(`TESTING Q${questionNum}: ${question.category}`);
  console.log('='.repeat(80));
  console.log('Correct Answer:', question.correct_answer);
  
  // Format with very clear structure
  const formattedQuery = `QUESTION:
${question.question}

OPTIONS:
A) ${question.options.A}
B) ${question.options.B}
C) ${question.options.C}
D) ${question.options.D}

Please analyze this question step by step, show your reasoning, and then clearly state your final answer as one of the letters: A, B, C, or D.`;

  const startTime = Date.now();
  
  try {
    console.log('\nProcessing with improved agent...\n');
    console.log('Waiting for FULL response to complete...\n');
    
    const result = await unifiedAgent.process(formattedQuery, {}, {
      requestId: `gpqa-q${questionNum}-retest-${Date.now()}`,
      timeout: 300000  // 5 minute timeout
    });
    
    const duration = Date.now() - startTime;
    
    console.log('\n--- FULL RESPONSE ---');
    console.log(result.content);
    console.log('\n--- END OF RESPONSE ---\n');
    
    // Manual extraction guidance
    console.log('=== MANUAL ANSWER EXTRACTION ===');
    console.log('Please manually read the response above and identify the answer.');
    console.log('Look for phrases like:');
    console.log('- "The answer is..."');
    console.log('- "Therefore,..."');
    console.log('- "The correct answer is..."');
    console.log('- "Answer: ..."');
    console.log('- Or the final letter mentioned\n');
    
    // Save full response for manual review
    const fullResult = {
      question_id: question.id,
      question_number: questionNum,
      category: question.category,
      test_date: new Date().toISOString(),
      model: process.env.OPENROUTER_MODEL,
      duration_ms: duration,
      duration_seconds: (duration / 1000).toFixed(1),
      steps: result.metadata.steps,
      tools_used: result.metadata.toolsUsed,
      correct_answer: question.correct_answer,
      full_response: result.content,
      response_length: result.content.length,
      manual_extraction_needed: true
    };
    
    const filename = `retest-gpqa-q${questionNum}-manual-${Date.now()}.json`;
    fs.writeFileSync(filename, JSON.stringify(fullResult, null, 2));
    
    console.log(`\n✓ Full response saved to ${filename}`);
    
    // Analysis
    console.log('\n--- ANALYSIS ---');
    console.log(`Duration: ${(duration / 1000).toFixed(1)}s`);
    console.log(`Response length: ${result.content.length} characters`);
    console.log(`Steps taken: ${result.metadata.steps.length}`);
    console.log(`Tools used: ${result.metadata.toolsUsed.join(' → ')}`);
    console.log(`Correct answer should be: ${question.correct_answer}`);
    
    return {
      questionNum,
      duration,
      tools: result.metadata.toolsUsed,
      responseLength: result.content.length
    };
    
  } catch (error) {
    console.error(`\n❌ ERROR testing Q${questionNum}:`, error.message);
    return {
      questionNum,
      error: error.message
    };
  }
}

async function main() {
  // Test Q11
  console.log('\n' + '*'.repeat(80));
  console.log('STARTING Q11 RETEST');
  console.log('*'.repeat(80));
  
  const q11Result = await testQuestion(q11, 11);
  
  console.log('\n\nWaiting 15 seconds before next question...\n');
  await new Promise(resolve => setTimeout(resolve, 15000));
  
  // Test Q13
  console.log('\n' + '*'.repeat(80));
  console.log('STARTING Q13 RETEST');
  console.log('*'.repeat(80));
  
  const q13Result = await testQuestion(q13, 13);
  
  // Summary
  console.log('\n' + '='.repeat(80));
  console.log('RETEST COMPLETE');
  console.log('='.repeat(80));
  
  console.log('\nPlease manually review the full responses above to extract the answers.');
  console.log('The correct answers are:');
  console.log('- Q11: C');
  console.log('- Q13: A');
}

main().catch(console.error);