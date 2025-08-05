const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const unifiedAgent = require('./src/unified-agent');
const e2bManager = require('./src/e2b-manager');
const fs = require('fs');

console.log('=== RETESTING WRONG GPQA QUESTIONS WITH GROK-3-MINI ===\n');
console.log('Model:', process.env.OPENROUTER_MODEL);
console.log('Date:', new Date().toISOString());

// The 4 wrong questions
const questions = [
  {
    id: "gpqa-1",
    category: "Molecular Biology",
    question: "A large gene has dozens of exons, of which the central ones code for folded triple helical repeats that connect the cytoskeleton with sarcolemma and extracellular space. Each exon usually codes for one folded triple alpha helix. The most common mutations of the gene are central exon deletions that create out-of-frame peptides and progressive degenerative organ waste. A solution is to deliver a Morpholino that recognizes the 5' end of the out-of-frame exon in pre-mRNA. The molecule prevents binding of the spliceosome and creates exon skipping and in-frame joining. Several missing exons are well tolerated by an organism. Which structure below is not involved in the proposed therapy?",
    options: {
      A: "antisense",
      B: "polyA tail", 
      C: "lariat",
      D: "R-loops"
    },
    correct_answer: "D"
  },
  {
    id: "gpqa-3",
    category: "Organic Chemistry",
    question: "trans-cinnamaldehyde was treated with methylmagnesium bromide, forming product 1. 1 was treated with pyridinium chlorochromate, forming product 2. 3 was treated with (dimethyl(oxo)-l6-sulfaneylidene)methane in DMSO at elevated temperature, forming product 3. how many carbon atoms are there in product 3?",
    options: {
      A: "11",
      B: "10",
      C: "12",
      D: "14"
    },
    correct_answer: "A"
  },
  {
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
  },
  {
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
  }
];

// Monitor E2B globally
let globalE2bCount = 0;
const originalExecute = e2bManager.executeWithFallback.bind(e2bManager);
e2bManager.executeWithFallback = async function(code, options) {
  globalE2bCount++;
  console.log(`\\n=== E2B EXECUTION #${globalE2bCount} (Global) ===`);
  const result = await originalExecute(code, options);
  return result;
};

async function testQuestion(q) {
  console.log(`\\n${'='.repeat(70)}`);
  console.log(`TESTING ${q.id.toUpperCase()}: ${q.category}`);
  console.log(`${'='.repeat(70)}\\n`);
  
  console.log('Correct Answer:', q.correct_answer, `(${q.options[q.correct_answer]})`);
  
  const formattedQuestion = `${q.question}

Options:
A) ${q.options.A}
B) ${q.options.B}
C) ${q.options.C}
D) ${q.options.D}`;
  
  console.log('\\nProcessing...\\n');
  
  const startTime = Date.now();
  let questionE2bCount = 0;
  const startE2bCount = globalE2bCount;
  
  try {
    const result = await unifiedAgent.process(formattedQuestion, {}, {
      requestId: `${q.id}-grok3mini-${Date.now()}`
    });
    
    questionE2bCount = globalE2bCount - startE2bCount;
    const elapsedTime = ((Date.now() - startTime) / 1000).toFixed(1);
    
    // Save response
    const filename = `${q.id}-grok3mini-response-${Date.now()}.txt`;
    fs.writeFileSync(filename, result.content);
    
    console.log('\\nResponse time:', elapsedTime, 's');
    console.log('E2B executions for this question:', questionE2bCount);
    console.log('Tools used:', result.metadata?.toolsUsed?.join(', ') || 'None');
    console.log('Response saved to:', filename);
    
    // Show last part of response
    console.log('\\n--- LAST 800 CHARACTERS ---');
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
    
    // For numeric answers, also check
    if (q.id === 'gpqa-3') {
      answerPatterns.push(
        /(\\d+)\\s*carbon atoms.*final/i,
        /final.*?(\\d+)\\s*carbon/i,
        /answer.*?(\\d+)/i
      );
    }
    
    let extractedAnswer = null;
    for (const pattern of answerPatterns) {
      const match = result.content.match(pattern);
      if (match) {
        extractedAnswer = match[1];
        // Map numbers to letters for Q3
        if (q.id === 'gpqa-3') {
          if (extractedAnswer === "11") extractedAnswer = "A";
          else if (extractedAnswer === "10") extractedAnswer = "B";
          else if (extractedAnswer === "12") extractedAnswer = "C";
          else if (extractedAnswer === "14") extractedAnswer = "D";
        }
        extractedAnswer = extractedAnswer.toUpperCase();
        break;
      }
    }
    
    console.log('\\n--- ANSWER CHECK ---');
    console.log('Correct answer:', q.correct_answer, `(${q.options[q.correct_answer]})`);
    console.log('Extracted answer:', extractedAnswer || 'MANUAL CHECK NEEDED');
    console.log('Status:', extractedAnswer === q.correct_answer ? 'CORRECT ✓' : 'INCORRECT ✗');
    
    return {
      question_id: q.id,
      category: q.category,
      correct_answer: q.correct_answer,
      extracted_answer: extractedAnswer,
      response_time: elapsedTime,
      e2b_executions: questionE2bCount,
      filename: filename,
      match: extractedAnswer === q.correct_answer,
      response_length: result.content.length,
      tools_used: result.metadata?.toolsUsed || []
    };
    
  } catch (error) {
    console.error('ERROR:', error.message);
    return {
      question_id: q.id,
      category: q.category,
      error: error.message,
      e2b_executions: questionE2bCount
    };
  }
}

async function runTests() {
  const results = [];
  
  for (const q of questions) {
    const result = await testQuestion(q);
    results.push(result);
    
    // Small delay between questions
    await new Promise(resolve => setTimeout(resolve, 2000));
  }
  
  // Summary
  console.log('\\n\\n' + '='.repeat(70));
  console.log('GROK-3-MINI RETEST SUMMARY');
  console.log('='.repeat(70) + '\\n');
  
  const correct = results.filter(r => r.match).length;
  const errors = results.filter(r => r.error).length;
  
  console.log(`Total questions: ${results.length}`);
  console.log(`Correct: ${correct}/${results.length} (${(correct/results.length*100).toFixed(0)}%)`);
  console.log(`Errors: ${errors}`);
  console.log(`Total E2B executions: ${globalE2bCount}`);
  console.log(`\\nDetailed results:`);
  
  results.forEach(r => {
    if (r.error) {
      console.log(`${r.question_id}: ERROR - ${r.error}`);
    } else {
      const status = r.match ? '✓' : '✗';
      console.log(`${r.question_id}: ${status} Model: ${r.extracted_answer || '?'}, Correct: ${r.correct_answer} (${r.response_time}s, E2B: ${r.e2b_executions})`);
    }
  });
  
  // Compare with Gemini results
  console.log('\\n--- COMPARISON WITH GEMINI-2.5-FLASH-LITE ---');
  console.log('Q1 (Molecular Bio): Gemini ✗ (B), Correct: D');
  console.log('Q3 (Organic Chem): Gemini ✓ (A on retest), Correct: A'); 
  console.log('Q4 (Molecular Bio): Gemini ✗ (A), Correct: B');
  console.log('Q6 (Chemistry): Gemini ✗ (C/134°), Correct: B/124°');
  
  // Save results
  const summaryFile = 'grok-3-mini-wrong-questions-summary.json';
  fs.writeFileSync(summaryFile, JSON.stringify(results, null, 2));
  console.log(`\\nResults saved to: ${summaryFile}`);
  
  await e2bManager.shutdown();
}

runTests().catch(console.error);