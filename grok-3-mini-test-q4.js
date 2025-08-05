const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const unifiedAgent = require('./src/unified-agent');
const e2bManager = require('./src/e2b-manager');
const fs = require('fs');

console.log('=== TESTING GPQA Q4 WITH GROK-3-MINI ===\n');
console.log('Model:', process.env.OPENROUTER_MODEL);
console.log('Date:', new Date().toISOString());

// Q4 data
const q4 = {
  id: "gpqa-4",
  category: "Molecular Biology",
  question: "To investigate the causes of a complex genetic disease, you culture patient cells and carry out DNA sequencing to detect mutations in candidate genes. This revealed a mutation in the gene HOXB2 that is only present in the patient cells and not the healthy controls. To learn more about the role of this mutation in the disease, you want to explore the relationship between chromatin structure and gene expression in patient cells and compare your results to healthy cells. Which of the following combinations of methods would provide you with results that would help your investigations?",
  options: {
    A: "CHIP-seq, RNA-seq, and qRT PCR",
    B: "CHIP-seq, chromosome conformation capture, and qRT-PCR",
    C: "ChIP-seq and RNA-seq",
    D: "Chromosome conformation capture and RNA-seq"
  },
  correct_answer: "B"
};

console.log('Correct Answer:', q4.correct_answer, `(${q4.options[q4.correct_answer]})`);

const formattedQuestion = `${q4.question}

Options:
A) ${q4.options.A}
B) ${q4.options.B}
C) ${q4.options.C}
D) ${q4.options.D}`;

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
      requestId: 'gpqa-q4-grok3mini-' + Date.now()
    });
    
    const elapsedTime = ((Date.now() - startTime) / 1000).toFixed(1);
    
    // Save response
    const filename = `gpqa-q4-grok3mini-response-${Date.now()}.txt`;
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
    
    console.log('\n--- ANSWER CHECK ---');
    console.log('Correct answer:', q4.correct_answer, `(${q4.options[q4.correct_answer]})`);
    console.log('Extracted answer:', extractedAnswer || 'MANUAL CHECK NEEDED');
    console.log('Status:', extractedAnswer === q4.correct_answer ? 'CORRECT ✓' : 'INCORRECT ✗');
    
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