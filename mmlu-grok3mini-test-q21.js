const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const unifiedAgent = require('./src/unified-agent');
const e2bManager = require('./src/e2b-manager');
const fs = require('fs');

console.log('=== TESTING MMLU Q21 (MATH - COMPOSITE FUNCTIONS) WITH GROK-3-MINI ===\n');
console.log('Model:', process.env.OPENROUTER_MODEL);
console.log('Date:', new Date().toISOString());

const q21 = {
  id: "mmlu-q21",
  category: "Mathematics",
  question: "Let f(x, y) = x² + y² and g(x, y) = x² - y². What is the value of f(g(x, y), g(f(x, y)))?",
  options: {
    A: "2x⁴ - 2y⁴",
    B: "x⁴ + y⁴",
    C: "x⁴ - y⁴",
    D: "2x⁴ + 2y⁴",
    E: "x² + y²",
    F: "2x⁴ + 2y⁴",
    G: "x⁴ + 2x²y² + y⁴",
    H: "x⁴ - 2x²y² + y⁴",
    I: "4x²y²",
    J: "0"
  },
  correct_answer: "B"
};

console.log('Correct Answer:', q21.correct_answer, `(${q21.options[q21.correct_answer]})`);

const formattedQuestion = `${q21.question}

Options:
A) ${q21.options.A}
B) ${q21.options.B}
C) ${q21.options.C}
D) ${q21.options.D}
E) ${q21.options.E}
F) ${q21.options.F}
G) ${q21.options.G}
H) ${q21.options.H}
I) ${q21.options.I}
J) ${q21.options.J}`;

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
      requestId: 'mmlu-q21-grok3mini-' + Date.now()
    });
    
    const elapsedTime = ((Date.now() - startTime) / 1000).toFixed(1);
    
    // Save response
    const filename = `mmlu-q21-grok3mini-response-${Date.now()}.txt`;
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
      /final answer is[:?\s]*\*?\*?\(?([A-J])\)?/i,
      /answer is[:?\s]*\*?\*?\(?([A-J])\)?/i,
      /correct answer[:?\s]*\*?\*?\(?([A-J])\)?/i,
      /\\boxed\{([A-J])\}/,
      /\*\*\(?([A-J])\)?\*\*/,
      /Therefore[,\s]*\(?([A-J])\)?/i
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
    console.log('Correct answer:', q21.correct_answer, `(${q21.options[q21.correct_answer]})`);
    console.log('Extracted answer:', extractedAnswer || 'MANUAL CHECK NEEDED');
    console.log('Status:', extractedAnswer === q21.correct_answer ? 'CORRECT ✓' : 'INCORRECT ✗');
    
    // Reference previous results
    console.log('\n--- PREVIOUS RESULTS ---');
    console.log('Gemini-2.5-flash-lite: D');
    console.log('O4-mini-high: D');
    console.log('Both were wrong (D=2x⁴ + 2y⁴, correct is B=x⁴ + y⁴)');
    
    return {
      question_id: q21.id,
      category: q21.category,
      correct_answer: q21.correct_answer,
      extracted_answer: extractedAnswer,
      response_time: elapsedTime,
      e2b_executions: e2bExecutionCount,
      filename: filename,
      match: extractedAnswer === q21.correct_answer
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