const workflowEngine = require('./src/workflow-engine');

// Expected answers
const expectedAnswers = ['D', 'D', 'A', 'B', 'A'];
const questions = [
  "A large gene has dozens of exons, of which the central ones code for folded triple helical repeats that connect the cytoskeleton with sarcolemma and extracellular space. Each exon usually codes for one folded triple alpha helix. The most common mutations of the gene are central exon deletions that create out-of-frame peptides and progressive degenerative organ waste. A solution is to deliver a Morpholino that recognizes the 5' end of the out-of-frame exon in pre-mRNA. The molecule prevents binding of the spliceosome and creates exon skipping and in-frame joining. Several missing exons are well tolerated by an organism. Which structure below is not involved in the proposed therapy?",
  
  "Two quantum states with energies E1 and E2 have a lifetime of 10^-9 sec and 10^-8 sec, respectively. We want to clearly distinguish these two energy levels. Which one of the following options could be their energy difference so that they can be clearly resolved?",
  
  "trans-cinnamaldehyde was treated with methylmagnesium bromide, forming product 1. 1 was treated with pyridinium chlorochromate, forming product 2. 3 was treated with (dimethyl(oxo)-l6-sulfaneylidene)methane in DMSO at elevated temperature, forming product 3. how many carbon atoms are there in product 3?",
  
  "To investigate the causes of a complex genetic disease, you culture patient cells and carry out DNA sequencing to detect mutations in candidate genes. This revealed a mutation in the gene HOXB2 that is only present in the patient cells and not the healthy controls. To learn more about the role of this mutation in the disease, you want to explore the relationship between chromatin structure and gene expression in patient cells and compare your results to healthy cells. Which of the following combinations of methods would provide you with results that would help your investigations?",
  
  "A spin-half particle is in a linear superposition 0.5|up⟩+√3/2|down⟩ of its spin-up and spin-down states. If |up⟩ and |down⟩ are the eigenstates of σ_z , then what is the expectation value up to one decimal place, of the operator 10σ_z+5σ_x ? Here, symbols have their usual meanings"
];

async function analyzeAnswer(questionNum, question, expectedAnswer) {
  console.log(`\n${'='.repeat(50)}`);
  console.log(`Q${questionNum}: Expected Answer = ${expectedAnswer}`);
  console.log(`${'='.repeat(50)}`);
  
  try {
    const response = await workflowEngine.answer(question);
    const content = response.content || '';
    
    // Manual analysis based on question type
    let actualAnswer = 'UNCLEAR';
    let reasoning = '';
    
    if (questionNum === 1) {
      // Q1: Looking for "R-loops" not involved
      if (content.toLowerCase().includes('r-loop') && 
          (content.toLowerCase().includes('not involved') || 
           content.toLowerCase().includes('not required') ||
           content.toLowerCase().includes('not directly'))) {
        actualAnswer = 'D';
        reasoning = 'Mentions R-loops as not involved in therapy';
      } else if (content.toLowerCase().includes('cytoskeleton') && 
                 content.toLowerCase().includes('not involved')) {
        actualAnswer = 'WRONG - mentioned cytoskeleton instead of R-loops';
        reasoning = 'Incorrectly focused on cytoskeleton rather than R-loops';
      }
    }
    
    if (questionNum === 2) {
      // Q2: Physics - should find 10^-4 eV
      if (content.includes('10^-4') || content.includes('1e-04') || content.includes('10⁻⁴')) {
        if (content.toLowerCase().includes('option d') || 
            content.toLowerCase().includes('answer d') ||
            content.includes('D (1e-04') ||
            content.includes('D)')) {
          actualAnswer = 'D';
          reasoning = 'Correctly calculated and identified 10^-4 eV as answer D';
        } else {
          actualAnswer = 'PARTIAL';
          reasoning = 'Calculated 10^-4 eV but unclear if identified as option D';
        }
      }
    }
    
    if (questionNum === 3) {
      // Q3: Chemistry - should count 11 carbons
      if (content.includes('11') && content.toLowerCase().includes('carbon')) {
        actualAnswer = 'A';
        reasoning = 'Correctly counted 11 carbons in final product';
      } else if (content.includes('10') && content.toLowerCase().includes('carbon')) {
        actualAnswer = 'B';
        reasoning = 'Counted 10 carbons (incorrect)';
      }
    }
    
    if (questionNum === 4) {
      // Q4: Should mention ChIP-seq, chromosome conformation, qRT-PCR
      const hasChipSeq = content.toLowerCase().includes('chip-seq');
      const hasChromConform = content.toLowerCase().includes('chromosome conformation');
      const hasQrtPcr = content.toLowerCase().includes('qrt-pcr') || content.toLowerCase().includes('qrt pcr');
      
      if (hasChipSeq && hasChromConform && hasQrtPcr) {
        actualAnswer = 'B';
        reasoning = 'Mentions all three: ChIP-seq, chromosome conformation, qRT-PCR';
      } else {
        reasoning = `Mentions: ChIP-seq=${hasChipSeq}, ChromConform=${hasChromConform}, qRT-PCR=${hasQrtPcr}`;
      }
    }
    
    if (questionNum === 5) {
      // Q5: Should calculate expectation value = -0.7
      if (content.includes('-0.7')) {
        actualAnswer = 'A';
        reasoning = 'Calculated expectation value -0.7';
      } else {
        reasoning = 'Need to check for mathematical calculation';
      }
    }
    
    console.log(`Extracted Answer: ${actualAnswer}`);
    console.log(`Reasoning: ${reasoning}`);
    console.log(`Match: ${actualAnswer === expectedAnswer ? '✅ CORRECT' : '❌ INCORRECT'}`);
    
    // Show key parts of response
    const lines = content.split('\n');
    console.log('\nKey response lines:');
    lines.forEach((line, i) => {
      const lowerLine = line.toLowerCase();
      if (lowerLine.includes('answer') || lowerLine.includes('option') || 
          lowerLine.includes('conclusion') || lowerLine.includes('therefore') ||
          (questionNum === 1 && lowerLine.includes('r-loop')) ||
          (questionNum === 2 && (line.includes('10^-4') || line.includes('1e-04'))) ||
          (questionNum === 3 && line.includes('11') && lowerLine.includes('carbon')) ||
          (questionNum === 5 && line.includes('-0.7'))) {
        console.log(`  Line ${i}: ${line.trim()}`);
      }
    });
    
    return actualAnswer === expectedAnswer;
    
  } catch (error) {
    console.error(`Error: ${error.message}`);
    return false;
  }
}

async function runAnalysis() {
  console.log('🔍 MANUAL ANSWER EXTRACTION AND VERIFICATION\n');
  
  let correctCount = 0;
  
  for (let i = 0; i < questions.length; i++) {
    const isCorrect = await analyzeAnswer(i + 1, questions[i], expectedAnswers[i]);
    if (isCorrect) correctCount++;
    
    // Brief pause
    await new Promise(resolve => setTimeout(resolve, 1000));
  }
  
  console.log(`\n${'='.repeat(50)}`);
  console.log(`FINAL SCORE: ${correctCount}/5 (${Math.round(correctCount/5*100)}%)`);
  console.log(`${'='.repeat(50)}`);
}

runAnalysis();