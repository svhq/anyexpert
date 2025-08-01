// Check specific GLM-4.5 responses to confirm answers
const fs = require('fs');

console.log('🔍 CHECKING SPECIFIC GLM-4.5 RESPONSES\n');

// Load retest results which have full_response
const retestResults = JSON.parse(fs.readFileSync('./glm-4.5-complete-final.json', 'utf8'));

console.log('Found retest results with full responses:\n');

retestResults.retest_results.forEach((result, i) => {
  console.log(`${'='.repeat(80)}`);
  console.log(`${i + 1}. ${result.id} - ${result.category.toUpperCase()}`);
  console.log(`Correct Answer: ${result.correct_answer}`);
  console.log(`Extracted Answer: ${result.selected_answer}`);
  console.log(`Marked as: ${result.is_correct ? 'CORRECT' : 'WRONG'}`);
  console.log(`${'='.repeat(80)}`);
  
  // Look for explicit answer in response
  const response = result.full_response;
  
  // Search for answer patterns
  const patterns = [
    /(?:correct answer is|answer is|the answer is)[:\s]*\*?\*?([A-J])\)/gi,
    /therefore[,\s]*(?:the )?(?:correct )?answer is[:\s]*\*?\*?([A-J])\)/gi,
    /based on this[,\s]*(?:the )?(?:correct )?answer is[:\s]*\*?\*?([A-J])\)/gi,
    /\*\*([A-J])\)\*\*/g,
    /option ([A-J])[)\s]*(?:is )?(?:the )?correct/gi
  ];
  
  let foundAnswers = [];
  patterns.forEach(pattern => {
    const matches = [...response.matchAll(pattern)];
    matches.forEach(match => {
      foundAnswers.push(match[1]);
    });
  });
  
  if (foundAnswers.length > 0) {
    console.log(`✅ FOUND IN RESPONSE: ${[...new Set(foundAnswers)].join(', ')}`);
    
    // Check if the most common found answer is correct
    const freq = {};
    foundAnswers.forEach(ans => freq[ans] = (freq[ans] || 0) + 1);
    const mostCommon = Object.entries(freq).sort((a, b) => b[1] - a[1])[0][0];
    
    console.log(`🎯 MOST COMMON: ${mostCommon}`);
    console.log(`✅ ACTUALLY CORRECT?: ${mostCommon === result.correct_answer ? 'YES' : 'NO'}`);
    
    if (mostCommon !== result.selected_answer) {
      console.log(`🔧 EXTRACTION ERROR: Extracted "${result.selected_answer}" but response says "${mostCommon}"`);
    }
    
    if (mostCommon === result.correct_answer && !result.is_correct) {
      console.log(`✨ MISSED CORRECT ANSWER: Should be marked correct!`);
    }
  } else {
    // Show key phrases for manual inspection
    console.log(`❌ NO CLEAR ANSWER FOUND`);
    
    // Look for conclusion sections
    const conclusionMatch = response.match(/conclusion[:\s]*([\s\S]{0,200})/gi);
    if (conclusionMatch) {
      console.log(`📝 CONCLUSION: ${conclusionMatch[0].substring(0, 150)}...`);
    }
    
    // Look for "correct answer" mentions
    const answerMentions = response.match(/correct answer[\s\S]{0,100}/gi);
    if (answerMentions) {
      console.log(`📝 ANSWER MENTIONS:`);
      answerMentions.forEach(mention => {
        console.log(`   - ${mention.substring(0, 80)}...`);
      });
    }
  }
  
  console.log('\n');
});

// Summary of findings
console.log('\n' + '█'.repeat(80));
console.log('SUMMARY OF FINDINGS');
console.log('█'.repeat(80));

let extractionErrors = 0;
let missedCorrectAnswers = 0;
let actuallyCorrect = 0;

retestResults.retest_results.forEach(result => {
  const response = result.full_response;
  const patterns = [
    /(?:correct answer is|answer is|the answer is)[:\s]*\*?\*?([A-J])\)/gi,
    /therefore[,\s]*(?:the )?(?:correct )?answer is[:\s]*\*?\*?([A-J])\)/gi,
    /based on this[,\s]*(?:the )?(?:correct )?answer is[:\s]*\*?\*?([A-J])\)/gi,
    /\*\*([A-J])\)\*\*/g,
    /option ([A-J])[)\s]*(?:is )?(?:the )?correct/gi
  ];
  
  let foundAnswers = [];
  patterns.forEach(pattern => {
    const matches = [...response.matchAll(pattern)];
    matches.forEach(match => {
      foundAnswers.push(match[1]);
    });
  });
  
  if (foundAnswers.length > 0) {
    const freq = {};
    foundAnswers.forEach(ans => freq[ans] = (freq[ans] || 0) + 1);
    const mostCommon = Object.entries(freq).sort((a, b) => b[1] - a[1])[0][0];
    
    if (mostCommon !== result.selected_answer) {
      extractionErrors++;
    }
    
    if (mostCommon === result.correct_answer) {
      actuallyCorrect++;
      if (!result.is_correct) {
        missedCorrectAnswers++;
      }
    }
  }
});

console.log(`\n📊 ANALYSIS RESULTS:`);
console.log(`- Questions analyzed: ${retestResults.retest_results.length}`);
console.log(`- Extraction errors: ${extractionErrors}`);
console.log(`- Missed correct answers: ${missedCorrectAnswers}`);
console.log(`- Actually correct answers: ${actuallyCorrect}`);

const correctedScore = actuallyCorrect;
console.log(`\n🎯 CORRECTED SCORE FOR RETEST QUESTIONS: ${correctedScore}/${retestResults.retest_results.length}`);

// Now need to add the original questions that were marked correct
const originalCorrect = 16; // From the summary
const totalCorrected = originalCorrect + correctedScore - 3; // Remove the 3 that were double-counted

console.log(`\n📈 ESTIMATED TOTAL CORRECTED SCORE:`);
console.log(`- Original correct: ${originalCorrect}`);
console.log(`- Additional from retest corrections: ${missedCorrectAnswers}`);
console.log(`- TOTAL: ${originalCorrect + missedCorrectAnswers}/30 (${Math.round((originalCorrect + missedCorrectAnswers)/30 * 100)}%)`);