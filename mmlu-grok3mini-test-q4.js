const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const unifiedAgent = require('./src/unified-agent');
const e2bManager = require('./src/e2b-manager');
const fs = require('fs');

console.log('=== TESTING MMLU Q4 (HISTORY - HANSEATIC LEAGUE) WITH GROK-3-MINI ===\n');
console.log('Model:', process.env.OPENROUTER_MODEL);
console.log('Date:', new Date().toISOString());

const q4 = {
  id: "mmlu-q4",
  category: "History",
  question: "The Hanseatic League was a commercial and defensive confederation of merchant guilds and their towns that dominated the Baltic maritime trade from the 13th to 17th centuries. Which of the following cities was NOT a major member of the Hanseatic League?",
  options: {
    A: "Lübeck",
    B: "Hamburg",
    C: "Bergen",
    D: "Novgorod",
    E: "Amsterdam",
    F: "Riga",
    G: "Danzig (Gdansk)",
    H: "Visby",
    I: "Cologne",
    J: "Stockholm"
  },
  correct_answer: "J"
};

console.log('Correct Answer:', q4.correct_answer, `(${q4.options[q4.correct_answer]})`);

const formattedQuestion = `${q4.question}

Options:
A) ${q4.options.A}
B) ${q4.options.B}
C) ${q4.options.C}
D) ${q4.options.D}
E) ${q4.options.E}
F) ${q4.options.F}
G) ${q4.options.G}
H) ${q4.options.H}
I) ${q4.options.I}
J) ${q4.options.J}`;

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
      requestId: 'mmlu-q4-grok3mini-' + Date.now()
    });
    
    const elapsedTime = ((Date.now() - startTime) / 1000).toFixed(1);
    
    // Save response
    const filename = `mmlu-q4-grok3mini-response-${Date.now()}.txt`;
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
      /Therefore[,\s]*\(?([A-J])\)?/i,
      /NOT a.*?member.*?\(?([A-J])\)?/i
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
    console.log('Correct answer:', q4.correct_answer, `(${q4.options[q4.correct_answer]})`);
    console.log('Extracted answer:', extractedAnswer || 'MANUAL CHECK NEEDED');
    console.log('Status:', extractedAnswer === q4.correct_answer ? 'CORRECT ✓' : 'INCORRECT ✗');
    
    // Reference previous results
    console.log('\n--- PREVIOUS RESULTS ---');
    console.log('Gemini-2.5-flash-lite: E');
    console.log('O4-mini-high: E');
    console.log('Both were wrong (E=Amsterdam, correct is J=Stockholm)');
    
    return {
      question_id: q4.id,
      category: q4.category,
      correct_answer: q4.correct_answer,
      extracted_answer: extractedAnswer,
      response_time: elapsedTime,
      e2b_executions: e2bExecutionCount,
      filename: filename,
      match: extractedAnswer === q4.correct_answer
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