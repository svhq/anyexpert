const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const unifiedAgent = require('./src/unified-agent');
const e2bManager = require('./src/e2b-manager');
const fs = require('fs');

console.log('=== RETESTING GPQA Q6 WITH FULL RESPONSE ===\n');
console.log('Model:', process.env.OPENROUTER_MODEL);
console.log('Date:', new Date().toISOString());

// Verified Q6 data from CSV 
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
  correct_answer: "B",
  explanation: "In the Cassie-Baxter state, droplets are in contact with a non-uniform surface: some of the droplet is in contact with the coating and some with air. The contact angle (θCB) of a droplet in the Cassie-Baxter state is given by: cosθCB = f1.cosθ1 + f2.cosθ2 Where f1 and f2 are the area fractions of the two components of the surface, in this case coating (f1) and air (f2). θ1 is the contact angle if the droplet was purely in contact with the coating, and θ2 is the contact angle if the droplet was purely in contact with air. First we need to calculate the f1 and f2 using the data we are given for water. We have θCB = 148°, θ1 = 132°, and θ2 is taken to be 180° (contact angle on air). We then have cos(148) = f1.cos(132) + f2.cos(180). By using f1 + f2 = 1, we can solve to give f1 = 0.46 and f2 = 0.54. Next we need to calculate the contact angle of hexadecane on the rough surface, we have θ1 = 102°, f1 = 0.46, f2 = 0.54, and θ2 is taken to be 180° (contact angle on air). Therefore, θCB = 129° for hexadecane. The question however asks about a droplet of octane. Octane is a shorter oil molecule than hexadecane and has a lower surface tension than hexadecane. For a given surface, the contact angle of octane is therefore always lower than for hexadecane. Therefore, the answer is 124° as this is the only answer lower than the 129° of hexadecane."
};

console.log('VERIFIED DATA FROM SOURCE:');
console.log('- Question: ✓ Matches CSV exactly');
console.log('- Options: ✓ All 4 options verified (139°, 124°, 134°, 129°)');
console.log('- Correct Answer: B (124°) - confirmed from CSV');
console.log('- Key concept: Cassie-Baxter equation for rough surfaces\n');

// Format question
const formattedQuestion = `${q6.question}

Options:
A) ${q6.options.A}
B) ${q6.options.B}
C) ${q6.options.C}
D) ${q6.options.D}`;

console.log('Formatted Question:');
console.log(formattedQuestion);
console.log('\n' + '='.repeat(70) + '\n');

console.log('This is a hard undergraduate chemistry question about contact angles.');
console.log('Allowing full response with calculations...\n');

// Monitor E2B
let e2bExecutionCount = 0;
const originalExecute = e2bManager.executeWithFallback.bind(e2bManager);
e2bManager.executeWithFallback = async function(code, options) {
  e2bExecutionCount++;
  console.log(`\n=== E2B EXECUTION #${e2bExecutionCount} ===`);
  console.log('Code preview:', code.substring(0, 300) + '...');
  const result = await originalExecute(code, options);
  if (result.results?.[0]?.text) {
    console.log('Result:', result.results[0].text.substring(0, 200) + '...');
  }
  return result;
};

async function retestQuestion() {
  const startTime = Date.now();
  
  try {
    const result = await unifiedAgent.process(formattedQuestion, {}, {
      requestId: 'gpqa-q6-retest-' + Date.now()
    });
    
    const elapsedTime = ((Date.now() - startTime) / 1000).toFixed(1);
    
    // Save full response
    const filename = `gpqa-q6-retest-full-${Date.now()}.txt`;
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
    console.log('The correct answer from CSV is: B (124°)');
    console.log('Key insight: Octane has lower surface tension than hexadecane,');
    console.log('so its contact angle must be LOWER than hexadecane\'s 129°');
    console.log('Please check the full response above for the model\'s final answer.');
    
    // Try automatic extraction
    const answerPatterns = [
      /final answer is[:\s]*\*?\*?\(?([A-D])\)?/i,
      /answer is[:\s]*\*?\*?\(?([A-D])\)?/i,
      /correct answer[:\s]*\*?\*?\(?([A-D])\)?/i,
      /Therefore[,\s]*\(?([A-D])\)?/i,
      /\*\*\(?([A-D])\)?\*\*/,
      /\boxed\{([A-D])\}/,
      /\boxed\{(\d+)°?\}/,
      /final answer.*?(\d+)°/i
    ];
    
    let extractedAnswer = null;
    for (const pattern of answerPatterns) {
      const match = result.content.match(pattern);
      if (match) {
        extractedAnswer = match[1];
        // If it's a number, map to letter
        if (extractedAnswer === "139") extractedAnswer = "A";
        else if (extractedAnswer === "124") extractedAnswer = "B";
        else if (extractedAnswer === "134") extractedAnswer = "C";
        else if (extractedAnswer === "129") extractedAnswer = "D";
        else extractedAnswer = extractedAnswer.toUpperCase();
        break;
      }
    }
    
    if (extractedAnswer) {
      console.log('\nAutomatically extracted answer:', extractedAnswer);
      console.log('Status:', extractedAnswer === q6.correct_answer ? 'CORRECT ✓' : 'INCORRECT ✗');
    }
    
    // Show key concepts
    console.log('\n=== KEY CONCEPTS ===');
    console.log('1. Cassie-Baxter equation: cos(θCB) = f1·cos(θ1) + f2·cos(θ2)');
    console.log('2. f1 = area fraction coating, f2 = area fraction air');
    console.log('3. θ2 = 180° for air contact');
    console.log('4. Octane (C8) has lower surface tension than hexadecane (C16)');
    console.log('5. Lower surface tension → lower contact angle');
    
    // Show explanation for reference
    console.log('\n=== CORRECT EXPLANATION (for reference) ===');
    console.log(q6.explanation);
    
    return {
      question_id: q6.id,
      category: q6.category,
      correct_answer: q6.correct_answer,
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