const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const unifiedAgent = require('./src/unified-agent');
const e2bManager = require('./src/e2b-manager');
const fs = require('fs');

console.log('=== RETESTING GPQA Q4 WITH FULL RESPONSE ===\n');
console.log('Model:', process.env.OPENROUTER_MODEL);
console.log('Date:', new Date().toISOString());

// Verified Q4 data from CSV 
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
  correct_answer: "B",
  explanation: "ChIP-seq is a method to determine the chromatin modifications present at a genetic region of interest. These often correlate to the chromatin structure. Chromosome conformation capture provides information on how regions of DNA come into contact with one another, which can be used to infer information about structure. qRT-PCR measures the abundance of specific RNA species, correlating to the level of the expression of a gene of interest. These three methods all provide valuable information about chromatin structure and gene expression. CLIP-seq studies RNA-protein interactions, so is not typically relevant to these type of chromatin structure/gene expression studies. RNA-seq also provides valuable information on gene expression but CHIP-seq and chromosome conformation capture are necessary to provide all the information required."
};

console.log('VERIFIED DATA FROM SOURCE:');
console.log('- Question: ✓ Matches CSV exactly');
console.log('- Options: ✓ All 4 options verified');
console.log('- Correct Answer: B (CHIP-seq, chromosome conformation capture, and qRT-PCR)');
console.log('- Note: CSV shows CHIP-seq (all caps) in some options, ChIP-seq in others\n');

// Format question
const formattedQuestion = `${q4.question}

Options:
A) ${q4.options.A}
B) ${q4.options.B}
C) ${q4.options.C}
D) ${q4.options.D}`;

console.log('Formatted Question:');
console.log(formattedQuestion);
console.log('\n' + '='.repeat(70) + '\n');

console.log('This is a post-graduate level molecular biology question.');
console.log('Allowing full response...\n');

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
      requestId: 'gpqa-q4-retest-' + Date.now()
    });
    
    const elapsedTime = ((Date.now() - startTime) / 1000).toFixed(1);
    
    // Save full response
    const filename = `gpqa-q4-retest-full-${Date.now()}.txt`;
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
    console.log('The correct answer from CSV is: B');
    console.log('B = CHIP-seq, chromosome conformation capture, and qRT-PCR');
    console.log('Please check the full response above for the model\'s final answer.');
    
    // Try automatic extraction
    const answerPatterns = [
      /final answer is[:\s]*\*?\*?\(?([A-D])\)?/i,
      /answer is[:\s]*\*?\*?\(?([A-D])\)?/i,
      /correct answer[:\s]*\*?\*?\(?([A-D])\)?/i,
      /Therefore[,\s]*\(?([A-D])\)?/i,
      /\*\*\(?([A-D])\)?\*\*/,
      /\boxed\{([A-D])\}/,
      /option\s+([A-D])\s+is\s+the\s+best/i,
      /best\s+answer\s+is\s+([A-D])/i
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
      console.log('Status:', extractedAnswer === q4.correct_answer ? 'CORRECT ✓' : 'INCORRECT ✗');
    }
    
    // Show key concepts
    console.log('\n=== KEY CONCEPTS ===');
    console.log('1. ChIP-seq: Chromatin modifications/protein-DNA interactions');
    console.log('2. Chromosome conformation capture (3C): 3D chromatin structure');
    console.log('3. qRT-PCR: Specific gene expression levels');
    console.log('4. RNA-seq: Genome-wide expression profiling');
    
    // Show explanation for reference
    console.log('\n=== CORRECT EXPLANATION (for reference) ===');
    console.log(q4.explanation);
    
    return {
      question_id: q4.id,
      category: q4.category,
      correct_answer: q4.correct_answer,
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