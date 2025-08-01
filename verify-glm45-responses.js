// Manual verification of all GLM-4.5 responses to ensure accurate extraction
const fs = require('fs');

// Load both original results and retest results
const originalResults = JSON.parse(fs.readFileSync('./glm-4.5-fixed-final.json', 'utf8'));
const retestResults = JSON.parse(fs.readFileSync('./glm-4.5-complete-final.json', 'utf8'));

console.log('🔍 MANUAL VERIFICATION OF ALL GLM-4.5 RESPONSES\n');
console.log('Checking each response to verify extraction accuracy...\n');

function manuallyExtractAnswer(responseText, questionId, correctAnswer) {
  console.log(`${'='.repeat(100)}`);
  console.log(`📋 VERIFYING: ${questionId}`);
  console.log(`Correct answer should be: ${correctAnswer}`);
  console.log(`${'='.repeat(100)}`);
  
  // Show key parts of response
  console.log('\n📄 RESPONSE ANALYSIS:');
  
  // Look for explicit answer statements
  const answerPatterns = [
    /(?:correct answer is|answer is|the answer is)[:\s]*\*?\*?([A-J])\)/gi,
    /therefore[,\s]*(?:the )?(?:correct )?answer is[:\s]*\*?\*?([A-J])\)/gi,
    /based on this[,\s]*(?:the )?(?:correct )?answer is[:\s]*\*?\*?([A-J])\)/gi,
    /option ([A-J])[)\s]*(?:is )?(?:the )?correct/gi,
    /\*\*([A-J])\)\*\*/g,
    /\*\*.*?([A-J])\).*?\*\*/g
  ];
  
  let foundAnswers = [];
  
  answerPatterns.forEach((pattern, i) => {
    const matches = [...responseText.matchAll(pattern)];
    if (matches.length > 0) {
      matches.forEach(match => {
        foundAnswers.push({
          pattern: i,
          letter: match[1],
          context: responseText.substring(Math.max(0, match.index - 50), match.index + 100)
        });
      });
    }
  });
  
  // Show what we found
  if (foundAnswers.length > 0) {
    console.log('✅ FOUND EXPLICIT ANSWERS:');
    foundAnswers.forEach((found, i) => {
      console.log(`${i + 1}. Letter "${found.letter}" - Context: "${found.context.trim()}"`);
    });
  } else {
    console.log('❌ NO EXPLICIT ANSWERS FOUND');
  }
  
  // Check for option analysis
  console.log('\n🔎 OPTION ANALYSIS IN RESPONSE:');
  const letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'];
  letters.forEach(letter => {
    const regex = new RegExp(`option ${letter}[^a-zA-Z]`, 'gi');
    const matches = responseText.match(regex);
    if (matches) {
      console.log(`- Option ${letter}: mentioned ${matches.length} times`);
    }
  });
  
  // Final determination
  let manualAnswer = null;
  if (foundAnswers.length > 0) {
    // Use the most frequent explicit answer
    const freq = {};
    foundAnswers.forEach(f => {
      freq[f.letter] = (freq[f.letter] || 0) + 1;
    });
    const mostFreq = Object.entries(freq).sort((a, b) => b[1] - a[1])[0];
    manualAnswer = mostFreq[0];
    console.log(`\n🎯 MANUAL DETERMINATION: ${manualAnswer} (appears ${mostFreq[1]} times)`);
  } else {
    console.log('\n⚠️ MANUAL REVIEW NEEDED - No clear answer found');
    // Show response snippet for manual review
    console.log('\n📝 RESPONSE SNIPPET (first 500 chars):');
    console.log(responseText.substring(0, 500) + '...');
    console.log('\n📝 RESPONSE SNIPPET (last 500 chars):');
    console.log('...' + responseText.substring(responseText.length - 500));
  }
  
  const isCorrect = manualAnswer === correctAnswer;
  console.log(`✅ CORRECT?: ${isCorrect ? 'YES' : 'NO'} (Manual: ${manualAnswer}, Correct: ${correctAnswer})`);
  
  return {
    manual_answer: manualAnswer,
    is_correct: isCorrect,
    found_answers: foundAnswers.length,
    needs_review: foundAnswers.length === 0
  };
}

// Combine all results
const allResults = [...originalResults.results];

// Add retest results, replacing original failed ones
retestResults.retest_results.forEach(retest => {
  const originalIndex = allResults.findIndex(r => r.id === retest.id);
  if (originalIndex >= 0) {
    allResults[originalIndex] = {
      ...retest,
      original_selected: allResults[originalIndex].selected_answer,
      retest: true
    };
  }
});

console.log(`Total results to verify: ${allResults.length}\n`);

let verificationResults = [];
let needsReview = [];
let extractionErrors = 0;
let correctedAnswers = 0;

async function runVerification() {
// Verify each result
for (let i = 0; i < allResults.length; i++) {
  const result = allResults[i];
  
  if (!result.full_response && !result.error) {
    console.log(`⚠️ Skipping ${result.id} - no full response available`);
    continue;
  }
  
  if (result.error) {
    console.log(`❌ Skipping ${result.id} - had error: ${result.error}`);
    continue;
  }
  
  const verification = manuallyExtractAnswer(
    result.full_response, 
    result.id, 
    result.correct_answer
  );
  
  verification.id = result.id;
  verification.category = result.category;
  verification.original_selected = result.selected_answer;
  verification.correct_answer = result.correct_answer;
  verification.original_correct = result.is_correct;
  
  verificationResults.push(verification);
  
  if (verification.needs_review) {
    needsReview.push(result.id);
  }
  
  // Check if extraction was wrong
  if (result.selected_answer !== verification.manual_answer) {
    extractionErrors++;
    console.log(`🔧 EXTRACTION ERROR: ${result.id} - Extracted: ${result.selected_answer}, Manual: ${verification.manual_answer}`);
  }
  
  // Check if we found a correct answer that was missed
  if (!result.is_correct && verification.is_correct) {
    correctedAnswers++;
    console.log(`✨ CORRECTION FOUND: ${result.id} - Was marked wrong, actually correct!`);
  }
  
  console.log('\n' + '-'.repeat(50) + '\n');
  
  // Small delay to make output readable
  await new Promise(resolve => setTimeout(resolve, 100));
}

// Final summary
console.log('\n' + '█'.repeat(100));
console.log('📊 VERIFICATION SUMMARY');
console.log('█'.repeat(100));

console.log(`\n🔍 Total responses verified: ${verificationResults.length}`);
console.log(`🔧 Extraction errors found: ${extractionErrors}`);
console.log(`✨ Corrections discovered: ${correctedAnswers}`);
console.log(`⚠️ Responses needing manual review: ${needsReview.length}`);

if (needsReview.length > 0) {
  console.log(`\nQuestions needing review: ${needsReview.join(', ')}`);
}

// Calculate corrected accuracy
const manuallyCorrect = verificationResults.filter(v => v.is_correct).length;
const totalVerified = verificationResults.length;
const correctedAccuracy = Math.round(manuallyCorrect / totalVerified * 100);

console.log(`\n📈 ACCURACY COMPARISON:`);
console.log(`Original extraction: ${originalResults.summary?.correct || 'unknown'}/29`);
console.log(`Manual verification: ${manuallyCorrect}/${totalVerified} (${correctedAccuracy}%)`);

if (correctedAnswers > 0) {
  console.log(`\n🎯 CORRECTED FINAL SCORE: ${manuallyCorrect}/30 (${Math.round(manuallyCorrect/30 * 100)}%)`);
}

// Save verification results
const verificationSummary = {
  model: 'z-ai/glm-4.5',
  verification_date: new Date().toISOString(),
  total_verified: totalVerified,
  extraction_errors: extractionErrors,
  corrections_found: correctedAnswers,
  manual_correct: manuallyCorrect,
  manual_accuracy: correctedAccuracy,
  needs_review: needsReview,
  verification_results: verificationResults,
  final_comparison: {
    'o4-mini-high': '29/30 (96.7%)',
    'GLM-4.5-air': '24/30 (80.0%)',
    'z-ai/glm-4.5-corrected': `${manuallyCorrect}/30 (${Math.round(manuallyCorrect/30 * 100)}%)`
  }
};

fs.writeFileSync('glm-4.5-verification-final.json', JSON.stringify(verificationSummary, null, 2));
console.log(`\n💾 Verification results saved to: glm-4.5-verification-final.json`);

return verificationSummary;
}

runVerification().catch(console.error);