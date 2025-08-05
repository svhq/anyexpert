const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const unifiedAgent = require('./src/unified-agent');
const e2bManager = require('./src/e2b-manager');
const fs = require('fs');

// Q36 - Polyomaviruses question
const question36 = {
  id: "mmlu-36",
  category: "health",
  question: `Which disease do polyomaviruses predominantly cause?

Options:
(A) Tumours
(B) Brain pathology
(C) No disease at all
(D) Kidney infections`,
  correct_answer: "C"
};

console.log('=== RETESTING MMLU QUESTION 36: HEALTH (POLYOMAVIRUSES) ===\n');
console.log('Correct Answer:', question36.correct_answer, '(No disease at all)');
console.log('\nQuestion:', question36.question);
console.log('\nNote: Most polyomaviruses are harmless in immunocompetent individuals.');
console.log('\nProcessing...\n');

async function testQuestion36() {
  const startTime = Date.now();
  
  try {
    const result = await unifiedAgent.process(question36.question, {}, {
      requestId: 'mmlu-q36-retest-' + Date.now()
    });
    
    const elapsedTime = ((Date.now() - startTime) / 1000).toFixed(1);
    
    // Save full response
    const filename = `mmlu-q36-retest-${Date.now()}.txt`;
    fs.writeFileSync(filename, result.content);
    
    console.log('\n=== RESULTS ===');
    console.log('Response time:', elapsedTime, 'seconds');
    console.log('Response saved to:', filename);
    
    // Show the ending
    console.log('\n=== LAST 800 CHARACTERS ===');
    console.log(result.content.substring(Math.max(0, result.content.length - 800)));
    
    console.log('\n=== VERIFICATION ===');
    console.log('The correct answer is C (No disease at all)');
    console.log('This is because most polyomaviruses are asymptomatic in healthy individuals.');
    
  } catch (error) {
    console.error('Error:', error.message);
  } finally {
    await e2bManager.shutdown();
  }
}

testQuestion36().catch(console.error);