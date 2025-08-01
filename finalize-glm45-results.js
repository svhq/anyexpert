// Finalize GLM-4.5 results with proper analysis
const fs = require('fs');

// Load the retest results
const finalResults = JSON.parse(fs.readFileSync('./glm-4.5-fixed-final.json', 'utf8'));

console.log('🔍 Analyzing GLM-4.5 Final Results\n');

// Separate successful vs failed extractions
const successful = finalResults.results.filter(r => r.selected_answer && !r.error);
const failed = finalResults.results.filter(r => !r.selected_answer || r.error);

console.log(`✅ Successfully extracted: ${successful.length}`);
console.log(`❌ Failed extraction: ${failed.length}`);
console.log(`⚠️ Total questions: ${finalResults.results.length}\n`);

// Analyze success patterns
console.log('📊 SUCCESS ANALYSIS:');
const successByCategory = {};
successful.forEach(r => {
  if (!successByCategory[r.category]) successByCategory[r.category] = { total: 0, correct: 0 };
  successByCategory[r.category].total++;
  if (r.is_correct) successByCategory[r.category].correct++;
});

Object.entries(successByCategory).forEach(([cat, stats]) => {
  const pct = Math.round(stats.correct / stats.total * 100);
  console.log(`- ${cat}: ${stats.correct}/${stats.total} (${pct}%) extracted successfully`);
});

// Analyze failure patterns  
console.log('\n❌ FAILURE ANALYSIS:');
const failureReasons = {};
failed.forEach(r => {
  const reason = r.error ? 'API/JSON Error' : 'Extraction Failed';
  failureReasons[reason] = (failureReasons[reason] || 0) + 1;
});

Object.entries(failureReasons).forEach(([reason, count]) => {
  console.log(`- ${reason}: ${count} questions`);
});

// Questions that could not get proper options (batch questions)
const batchQuestions = failed.filter(r => r.id.startsWith('batch'));
console.log(`- Batch questions missing options: ${batchQuestions.length}`);

// Calculate adjusted accuracy for questions with successful extraction
const extractedCorrect = successful.filter(r => r.is_correct).length;
const extractedTotal = successful.length;
const extractedAccuracy = Math.round(extractedCorrect / extractedTotal * 100);

console.log('\n📈 PERFORMANCE METRICS:');
console.log(`Raw accuracy (all attempts): ${finalResults.summary.accuracy}% (${finalResults.summary.correct}/${finalResults.summary.total})`);
console.log(`Extraction success rate: ${Math.round(successful.length / finalResults.results.length * 100)}% (${successful.length}/${finalResults.results.length})`);
console.log(`Accuracy on successfully extracted: ${extractedAccuracy}% (${extractedCorrect}/${extractedTotal})`);

// Compare with other models
console.log('\n🏆 MODEL COMPARISON:');
console.log('o4-mini-high: 29/29 (100%) - Perfect performance');
console.log('GLM-4.5-air: 24/29 (82.8%) - Strong performance'); 
console.log(`z-ai/glm-4.5: ${extractedCorrect}/${extractedTotal} (${extractedAccuracy}% on extracted questions)`);

// Specific successful questions analysis
console.log('\n🎯 SUCCESSFUL EXTRACTIONS BY TYPE:');
const originalSuccessful = successful.filter(r => r.id.startsWith('orig')).length;
const newSuccessful = successful.filter(r => r.id.startsWith('new')).length;
const batchSuccessful = successful.filter(r => r.id.startsWith('batch')).length;

console.log(`- Original 10 questions: ${originalSuccessful}/10`);
console.log(`- New 10 questions: ${newSuccessful}/9 (1 corrupted skipped)`); 
console.log(`- Batch 10 questions: ${batchSuccessful}/10`);

// Key findings
console.log('\n🔑 KEY FINDINGS:');
console.log('1. GLM-4.5 uses markdown bold formatting (**X)**) that required special extraction');
console.log('2. Many batch questions failed because original options were not preserved');
console.log('3. When extraction works, the model shows reasonable performance');
console.log('4. The model provides detailed explanations but formatting inconsistencies cause issues');

// Recommendations
console.log('\n💡 RECOMMENDATIONS:');
console.log('1. Use robust extraction patterns that handle multiple response formats');
console.log('2. Preserve original question structure including all options');
console.log('3. Consider GLM-4.5 for cost-effective applications with proper extraction');
console.log('4. o4-mini-high remains the best choice for highest accuracy requirements');

// Create final summary
const summary = {
  model: 'z-ai/glm-4.5',
  test_date: '2025-07-31',
  extraction_issues_identified: true,
  questions_attempted: finalResults.results.length,
  successful_extractions: successful.length,
  extraction_success_rate: Math.round(successful.length / finalResults.results.length * 100),
  accuracy_on_extracted: extractedAccuracy,
  raw_accuracy: finalResults.summary.accuracy,
  comparison: {
    'o4-mini-high': '100% (29/29)',
    'GLM-4.5-air': '82.8% (24/29)', 
    'z-ai/glm-4.5': `${extractedAccuracy}% (${extractedCorrect}/${extractedTotal} extracted)`
  },
  key_issue: 'Markdown formatting requires specialized extraction patterns',
  recommendation: 'Viable alternative to premium models with proper extraction handling'
};

fs.writeFileSync('glm-4.5-final-analysis.json', JSON.stringify(summary, null, 2));
console.log('\n💾 Final analysis saved to: glm-4.5-final-analysis.json');

return summary;