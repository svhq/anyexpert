const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const unifiedAgent = require('./src/unified-agent');
const e2bManager = require('./src/e2b-manager');
const fs = require('fs');

console.log('=== TESTING MMLU Q36 (HEALTH - POLYOMAVIRUSES) WITH GROK-3-MINI ===\n');
console.log('Model:', process.env.OPENROUTER_MODEL);
console.log('Date:', new Date().toISOString());

const q36 = {
  id: "mmlu-q36",
  category: "Health",
  question: "Which disease do polyomaviruses predominantly cause?",
  options: {
    A: "Tumours",
    B: "Brain pathology",
    C: "No disease at all",
    D: "Kidney infections"
  },
  correct_answer: "C"
};

console.log('Correct Answer:', q36.correct_answer, `(${q36.options[q36.correct_answer]})`);

const formattedQuestion = `${q36.question}

Options:
A) ${q36.options.A}
B) ${q36.options.B}
C) ${q36.options.C}
D) ${q36.options.D}`;

console.log('\nProcessing...\n');

// Monitor E2B
let e2bExecutionCount = 0;
const originalExecute = e2bManager.executeWithFallback.bind(e2bManager);
e2bManager.executeWithFallback = async function(code, options) {
  e2bExecutionCount++;
  console.log(`\n=== E2B EXECUTION #${e2bExecutionCount} ===`);
  const result = await originalExecute(code, options);
  return result;
};

async function testQuestion() {
  const startTime = Date.now();
  
  try {
    const result = await unifiedAgent.process(formattedQuestion, {}, {
      requestId: 'mmlu-q36-grok3mini-' + Date.now()
    });
    
    const elapsedTime = ((Date.now() - startTime) / 1000).toFixed(1);
    
    // Save response
    const filename = `mmlu-q36-grok3mini-response-${Date.now()}.txt`;
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
      /final answer is[:?\s]*\*?\*?\(?([A-D])\)?/i,
      /answer is[:?\s]*\*?\*?\(?([A-D])\)?/i,
      /correct answer[:?\s]*\*?\*?\(?([A-D])\)?/i,
      /\\boxed\{([A-D])\}/,
      /\*\*\(?([A-D])\)?\*\*/,
      /Therefore[,\s]*\(?([A-D])\)?/i,
      /predominantly.*?\(?([A-D])\)?/i
    ];
    
    let extractedAnswer = null;
    for (const pattern of answerPatterns) {
      const match = result.content.match(pattern);
      if (match) {
        extractedAnswer = match[1].toUpperCase();
        break;
      }
    }
    
    console.log('\n--- ANSWER CHECK ---');
    console.log('Correct answer:', q36.correct_answer, `(${q36.options[q36.correct_answer]})`);
    console.log('Extracted answer:', extractedAnswer || 'MANUAL CHECK NEEDED');
    console.log('Status:', extractedAnswer === q36.correct_answer ? 'CORRECT ✓' : 'INCORRECT ✗');
    
    // Reference previous results
    console.log('\n--- PREVIOUS RESULTS ---');
    console.log('Gemini-2.5-flash-lite: B (Brain pathology)');
    console.log('O4-mini-high: B (Brain pathology)');
    console.log('Both were wrong (B=Brain pathology, correct is C=No disease at all)');
    
    return {
      question_id: q36.id,
      category: q36.category,
      correct_answer: q36.correct_answer,
      extracted_answer: extractedAnswer,
      response_time: elapsedTime,
      e2b_executions: e2bExecutionCount,
      filename: filename,
      match: extractedAnswer === q36.correct_answer
    };
    
  } catch (error) {
    console.error('ERROR:', error.message);
    return { error: error.message };
  } finally {
    await e2bManager.shutdown();
  }
}

testQuestion().then(result => {
  console.log('\n=== SUMMARY ===');
  console.log(JSON.stringify(result, null, 2));
}).catch(console.error);