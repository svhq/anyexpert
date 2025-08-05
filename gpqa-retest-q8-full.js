const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const unifiedAgent = require('./src/unified-agent');
const e2bManager = require('./src/e2b-manager');
const fs = require('fs');

console.log('=== RETESTING GPQA Q8 WITH FULL RESPONSE ===\n');
console.log('Model:', process.env.OPENROUTER_MODEL);
console.log('Date:', new Date().toISOString());

// Verified Q8 data from CSV 
const q8 = {
  id: "gpqa-8",
  category: "Electromagnetism and Photonics",
  question: "In a parallel universe where a magnet can have an isolated North or South pole, Maxwell's equations look different. But, specifically, which of those equations are different?",
  options: {
    A: "The one related to the divergence of the magnetic field.",
    B: "The ones related to the divergence and the curl of the magnetic field.",
    C: "The ones related to the circulation of the electric field and the divergence of the magnetic field.",
    D: "The one related to the circulation of the magnetic field and the flux of the electric field."
  },
  correct_answer: "C",
  explanation: "Let's call E and B the electric and magnetic fields, respectively: The ones related to the circulation of the electric field and the divergence of the magnetic field is correct, since knowing that magnets can have an isolated pole means that magnetic monopoles exist and, thus, the contributions of magnetic charges and magnetic currents must be included in the equations. The way to include them is to \"symmetry-copy\" the other equations, with the following dictionary: E <-> B; electric charge <-> magnetic charge; electric current <-> magnetic current. In this way, the equations that become modified, with added terms, are the ones related to the circulation (or curl, in differential form) of E, and to the divergence (or flux in integral form) of B. The ones related to the divergence and the curl of the magnetic field is incorrect, because the one with the curl does not change, since it already includes all symmetric contributions appearing in its symmetric equation (curl of electric field). The one related to the divergence of the magnetic field is incorrect because that equation does get changed, but it's not the only one; the equation for the curl (or circulation) of E also changes. The one related to the circulation of the magnetic field and the flux of the electric field is incorrect because none of those equations are changed, since they already include the symmetric terms appearing in their symmetric equations (circulation of E and flux of B)."
};

console.log('VERIFIED DATA FROM SOURCE:');
console.log('- Question: ✓ Matches CSV exactly');
console.log('- Options: ✓ All 4 options verified');
console.log('- Correct Answer: C (circulation of E field and divergence of B field)');
console.log('- Key concept: Magnetic monopoles require symmetry in Maxwell equations\n');

// Format question
const formattedQuestion = `${q8.question}

Options:
A) ${q8.options.A}
B) ${q8.options.B}
C) ${q8.options.C}
D) ${q8.options.D}`;

console.log('Formatted Question:');
console.log(formattedQuestion);
console.log('\n' + '='.repeat(70) + '\n');

console.log('This is a conceptual physics question about Maxwell\'s equations.');
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
      requestId: 'gpqa-q8-retest-' + Date.now()
    });
    
    const elapsedTime = ((Date.now() - startTime) / 1000).toFixed(1);
    
    // Save full response
    const filename = `gpqa-q8-retest-full-${Date.now()}.txt`;
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
    console.log('The correct answer from CSV is: C');
    console.log('C = The ones related to the circulation of the electric field');
    console.log('    and the divergence of the magnetic field');
    console.log('Please check the full response above for the model\'s final answer.');
    
    // Try automatic extraction
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
      console.log('Status:', extractedAnswer === q8.correct_answer ? 'CORRECT ✓' : 'INCORRECT ✗');
    }
    
    // Show Maxwell's equations reference
    console.log('\n=== MAXWELL\'S EQUATIONS REFERENCE ===');
    console.log('Standard (no monopoles):');
    console.log('1. ∇·E = ρ/ε₀           (Gauss\'s law)');
    console.log('2. ∇·B = 0              (No magnetic monopoles)');
    console.log('3. ∇×E = -∂B/∂t         (Faraday\'s law)');
    console.log('4. ∇×B = μ₀(J + ε₀∂E/∂t) (Ampère-Maxwell law)');
    console.log('\nWith monopoles:');
    console.log('1. ∇·E = ρₑ/ε₀          (unchanged)');
    console.log('2. ∇·B = μ₀ρₘ           (CHANGED - magnetic charge density)');
    console.log('3. ∇×E = -∂B/∂t - μ₀Jₘ  (CHANGED - magnetic current)');
    console.log('4. ∇×B = μ₀(Jₑ + ε₀∂E/∂t) (unchanged)');
    
    // Show key insight
    console.log('\n=== KEY INSIGHT ===');
    console.log('Magnetic monopoles introduce SYMMETRY between E and B fields.');
    console.log('The two equations that change are:');
    console.log('- Divergence of B (was 0, now includes magnetic charge)');
    console.log('- Curl of E (now includes magnetic current term)');
    console.log('This is option C.');
    
    // Show explanation for reference
    console.log('\n=== CORRECT EXPLANATION (for reference) ===');
    console.log(q8.explanation);
    
    return {
      question_id: q8.id,
      category: q8.category,
      correct_answer: q8.correct_answer,
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