const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const unifiedAgent = require('./src/unified-agent');
const e2bManager = require('./src/e2b-manager');
const fs = require('fs');

console.log('=== TESTING GPQA Q6 WITH GROK-3-MINI ===\n');
console.log('Model:', process.env.OPENROUTER_MODEL);
console.log('Date:', new Date().toISOString());

// Q6 data
const q6 = {
  id: "gpqa-6",
  category: "Chemistry (general)",
  question: "A coating is applied to a substrate resulting in a perfectly smooth surface. The measured contact angles of this smooth coating are 132° and 102° for water and hexadecane respectively. The coating formulation is then modified and when now applied to the same type of substrate, a rough surface is produced. When a droplet of water or oil sits on the rough surface, the wettability of the surface can now be described by the Cassie-Baxter state. The water contact angle on the rough surface is now 148°. What would be the best estimate of the contact angle of a droplet of octane on the rough surface?",
  options: {
    A: "139°",
    B: "124°",
    C: "134°",
    D: "129°"
  },
  correct_answer: "B"
};

console.log('Correct Answer:', q6.correct_answer, `(${q6.options[q6.correct_answer]})`);

const formattedQuestion = `${q6.question}

Options:
A) ${q6.options.A}
B) ${q6.options.B}
C) ${q6.options.C}
D) ${q6.options.D}`;

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
      requestId: 'gpqa-q6-grok3mini-' + Date.now()
    });
    
    const elapsedTime = ((Date.now() - startTime) / 1000).toFixed(1);
    
    // Save response
    const filename = `gpqa-q6-grok3mini-response-${Date.now()}.txt`;
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
      /best estimate.*?([A-D])\)/i  // For this specific question
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
    console.log('Correct answer:', q6.correct_answer, `(${q6.options[q6.correct_answer]})`);
    console.log('Extracted answer:', extractedAnswer || 'MANUAL CHECK NEEDED');
    console.log('Status:', extractedAnswer === q6.correct_answer ? 'CORRECT ✓' : 'INCORRECT ✗');
    
    return {
      question_id: q6.id,
      category: q6.category,
      correct_answer: q6.correct_answer,
      extracted_answer: extractedAnswer,
      response_time: elapsedTime,
      e2b_executions: e2bExecutionCount,
      filename: filename,
      match: extractedAnswer === q6.correct_answer
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