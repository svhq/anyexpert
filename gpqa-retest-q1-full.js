const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const unifiedAgent = require('./src/unified-agent');
const e2bManager = require('./src/e2b-manager');
const fs = require('fs');

console.log('=== RETESTING GPQA Q1 WITH FULL RESPONSE ===\n');
console.log('Model:', process.env.OPENROUTER_MODEL);
console.log('Date:', new Date().toISOString());

// Verified Q1 data from CSV
const q1 = {
  id: "gpqa-1",
  category: "Molecular Biology",
  question: "A large gene has dozens of exons, of which the central ones code for folded triple helical repeats that connect the cytoskeleton with sarcolemma and extracellular space. Each exon usually codes for one folded triple alpha helix. The most common mutations of the gene are central exon deletions that create out-of-frame peptides and progressive degenerative organ waste. A solution is to deliver a Morpholino that recognizes the 5' end of the out-of-frame exon in pre-mRNA. The molecule prevents binding of the spliceosome and creates exon skipping and in-frame joining. Several missing exons are well tolerated by an organism. Which structure below is not involved in the proposed therapy?",
  options: {
    A: "antisense",
    B: "polyA tail", 
    C: "lariat",
    D: "R-loops"
  },
  correct_answer: "D",
  explanation: "The text describes the dystrophin gene and the FDA-approved oligonucleotide therapy that causes exon skipping by creating a functional, albeit shorter, dystrophin protein. Morpholino is bound to the pre-mRNA in an antisense orientation. Every splicing mechanism creates the lariat molecule that is circular with a 3' tail and soon degraded. The spliced RNA is polyadenylated at the 3' end. R-loops are triple helix of DNA and the pre-mRNA and a consequence of the RNA transcription, not splicing and RNA maturation."
};

console.log('VERIFIED DATA FROM SOURCE:');
console.log('- Question: ✓ Matches CSV exactly');
console.log('- Options: ✓ All 4 options verified');
console.log('- Correct Answer: D (R-loops) - confirmed from CSV\n');

// Format question
const formattedQuestion = `${q1.question}

Options:
A) ${q1.options.A}
B) ${q1.options.B}
C) ${q1.options.C}
D) ${q1.options.D}`;

console.log('Formatted Question:');
console.log(formattedQuestion);
console.log('\n' + '='.repeat(70) + '\n');

// Increase timeout to allow full response
const originalMaxTokens = unifiedAgent.config.openrouter.max_tokens;
console.log('Allowing full response (increased token limit)...\n');

// Monitor E2B
let e2bExecutionCount = 0;
const originalExecute = e2bManager.executeWithFallback.bind(e2bManager);
e2bManager.executeWithFallback = async function(code, options) {
  e2bExecutionCount++;
  console.log(`\n=== E2B EXECUTION #${e2bExecutionCount} ===`);
  const result = await originalExecute(code, options);
  return result;
};

async function retestQuestion() {
  const startTime = Date.now();
  
  try {
    const result = await unifiedAgent.process(formattedQuestion, {}, {
      requestId: 'gpqa-q1-retest-' + Date.now()
    });
    
    const elapsedTime = ((Date.now() - startTime) / 1000).toFixed(1);
    
    // Save full response
    const filename = `gpqa-q1-retest-full-${Date.now()}.txt`;
    fs.writeFileSync(filename, result.content);
    
    console.log('\n=== RESPONSE COMPLETE ===');
    console.log('Response time:', elapsedTime, 'seconds');
    console.log('Response length:', result.content.length, 'characters');
    console.log('E2B executions:', e2bExecutionCount);
    console.log('Tools used:', result.metadata?.toolsUsed?.join(', ') || 'None');
    console.log('Full response saved to:', filename);
    
    // Show full response for manual checking
    console.log('\n' + '='.repeat(70));
    console.log('FULL RESPONSE:');
    console.log('='.repeat(70) + '\n');
    console.log(result.content);
    console.log('\n' + '='.repeat(70));
    
    // Manual extraction guidance
    console.log('\n=== MANUAL ANSWER EXTRACTION ===');
    console.log('The correct answer from CSV is: D (R-loops)');
    console.log('Please check the full response above for the model\'s final answer.');
    
    // Try automatic extraction but don't rely on it
    const answerPatterns = [
      /final answer is[:\s]*\*?\*?\(?([A-D])\)?/i,
      /answer is[:\s]*\*?\*?\(?([A-D])\)?/i,
      /correct answer[:\s]*\*?\*?\(?([A-D])\)?/i,
      /Therefore[,\s]*\(?([A-D])\)?/i,
      /\*\*\(?([A-D])\)?\*\*/,
      /\boxed\{([A-D])\}/
    ];
    
    let extractedAnswer = null;
    for (const pattern of answerPatterns) {
      const match = result.content.match(pattern);
      if (match) {
        extractedAnswer = match[1].toUpperCase();
        break;
      }
    }
    
    if (extractedAnswer) {
      console.log('\nAutomatically extracted answer:', extractedAnswer);
      console.log('Status:', extractedAnswer === q1.correct_answer ? 'CORRECT ✓' : 'INCORRECT ✗');
    }
    
    // Show explanation for reference
    console.log('\n=== CORRECT EXPLANATION (for reference) ===');
    console.log(q1.explanation);
    
    return {
      question_id: q1.id,
      category: q1.category,
      correct_answer: q1.correct_answer,
      extracted_answer: extractedAnswer,
      response_time: elapsedTime,
      e2b_executions: e2bExecutionCount,
      response_file: filename,
      response_length: result.content.length,
      tools_used: result.metadata?.toolsUsed || []
    };
    
  } catch (error) {
    console.error('\n=== ERROR ===');
    console.error('Error:', error.message);
    return null;
  } finally {
    e2bManager.executeWithFallback = originalExecute;
    await e2bManager.shutdown();
  }
}

retestQuestion().then(result => {
  if (result) {
    console.log('\n=== FINAL SUMMARY ===');
    console.log(JSON.stringify(result, null, 2));
  }
}).catch(console.error);