// Complete retest of all failed GLM-4.5 questions with manual extraction
const fs = require('fs');
const workflowEngine = require('./src/workflow-engine');

// Get ALL 30 questions including the corrupted one
function getAllQuestions() {
  const allQuestions = [];
  
  // Original 10  
  const original = JSON.parse(fs.readFileSync('./mmlu_test_results_2025-07-30T16-45-40-316Z.json', 'utf8'));
  original.results.forEach((r, idx) => {
    allQuestions.push({
      id: `orig-${idx + 1}`,
      category: r.category,
      question: r.question,
      options: r.options,
      correct_answer: r.correct_answer,
      query_sent: r.query_sent
    });
  });
  
  // Get batch questions from the working GLM-4.5-air test that succeeded
  const airResults = JSON.parse(fs.readFileSync('./mmlu_batch_test_2025-07-31T02-07-13-672Z.json', 'utf8'));
  
  // Batch 1 - reconstruct from the 10-new-questions file which has the batch questions
  const newQuestions = JSON.parse(fs.readFileSync('./mmlu-10-new-questions.json', 'utf8'));
  
  // First 5 from new questions are batch1
  for (let i = 0; i < 5; i++) {
    const q = newQuestions.questions[i];
    allQuestions.push({
      id: `batch1-${i + 1}`,
      category: q.category,
      question: q.question,
      options: q.options,
      correct_answer: q.correct_answer,
      query_sent: `${q.question}\n\nOptions:\n${q.options.join('\n')}\n\nPlease select the correct answer and explain your reasoning.`
    });
  }
  
  // Next 5 from new questions are batch2  
  for (let i = 5; i < 10; i++) {
    const q = newQuestions.questions[i];
    allQuestions.push({
      id: `batch2-${i - 4}`,
      category: q.category,
      question: q.question,
      options: q.options,
      correct_answer: q.correct_answer,
      query_sent: `${q.question}\n\nOptions:\n${q.options.join('\n')}\n\nPlease select the correct answer and explain your reasoning.`
    });
  }
  
  // New 10 questions (including the corrupted one to test)
  newQuestions.questions.forEach((q, idx) => {
    allQuestions.push({
      id: `new-${idx + 1}`,
      category: q.category,
      question: q.question,
      options: q.options,
      correct_answer: q.correct_answer,
      query_sent: `${q.question}\n\nOptions:\n${q.options.join('\n')}\n\nPlease select the correct answer and explain your reasoning.`,
      is_corrupted: q.correct_answer === 'I' && !q.options.find(o => o.startsWith('I'))
    });
  });
  
  return allQuestions;
}

// Manual extraction with extensive logging
function manualExtractAnswer(responseContent, questionId, correctAnswer) {
  console.log(`\n${'='.repeat(80)}`);
  console.log(`🔍 MANUAL EXTRACTION for ${questionId}`);
  console.log(`Correct answer should be: ${correctAnswer}`);
  console.log(`Response length: ${responseContent.length} characters`);
  console.log(`${'='.repeat(80)}`);
  
  // Show full response for manual verification
  console.log('\n📄 FULL RESPONSE:');
  console.log(responseContent);
  console.log(`\n${'-'.repeat(80)}`);
  
  // Try multiple extraction patterns with logging
  const patterns = [
    { name: 'Bold Option', regex: /\*\*([A-J])\)/g },
    { name: 'Bold Answer', regex: /\*\*.*?([A-J])\).*?\*\*/g },
    { name: 'Correct Answer Is', regex: /correct answer is[:\s]*([A-J])\)/gi },
    { name: 'Answer:', regex: /Answer[:\s]+([A-J])\)/gi },
    { name: 'The Answer', regex: /The answer[:\s]+([A-J])\)/gi },
    { name: 'Select', regex: /select[:\s]+([A-J])\)/gi },
    { name: 'Choose', regex: /choose[:\s]+([A-J])\)/gi },
    { name: 'Option', regex: /option[:\s]+([A-J])\)/gi },
    { name: 'Letter Only', regex: /\b([A-J])\)/g }
  ];
  
  let extractedAnswer = null;
  let matchedPattern = 'none';
  let allMatches = [];
  
  console.log('\n🔎 PATTERN MATCHING RESULTS:');
  
  patterns.forEach(pattern => {
    const matches = [...responseContent.matchAll(pattern.regex)];
    if (matches.length > 0) {
      const letters = matches.map(m => m[1]);
      console.log(`- ${pattern.name}: Found ${letters.join(', ')}`);
      allMatches.push(...letters);
      
      if (!extractedAnswer) {
        extractedAnswer = letters[0];
        matchedPattern = pattern.name;
      }
    } else {
      console.log(`- ${pattern.name}: No matches`);
    }
  });
  
  // Frequency analysis
  if (allMatches.length > 0) {
    const freq = {};
    allMatches.forEach(letter => {
      freq[letter] = (freq[letter] || 0) + 1;
    });
    
    console.log('\n📊 FREQUENCY ANALYSIS:');
    Object.entries(freq).sort((a, b) => b[1] - a[1]).forEach(([letter, count]) => {
      console.log(`- ${letter}: ${count} occurrences`);
    });
    
    // Use most frequent if no clear pattern match
    if (!extractedAnswer) {
      const mostFrequent = Object.entries(freq).sort((a, b) => b[1] - a[1])[0];
      extractedAnswer = mostFrequent[0];
      matchedPattern = 'most frequent';
    }
  }
  
  console.log(`\n🎯 EXTRACTED ANSWER: ${extractedAnswer || 'NULL'}`);
  console.log(`🔧 MATCHED PATTERN: ${matchedPattern}`);
  console.log(`✅ CORRECT?: ${extractedAnswer === correctAnswer ? 'YES' : 'NO'}`);
  
  // Manual verification prompt
  if (!extractedAnswer) {
    console.log('\n⚠️ NO ANSWER EXTRACTED - MANUAL REVIEW NEEDED');
    console.log('Please review the full response above and determine the answer');
  }
  
  return {
    extracted: extractedAnswer,
    pattern: matchedPattern,
    all_matches: allMatches,
    is_correct: extractedAnswer === correctAnswer
  };
}

async function retestFailedQuestion(q) {
  console.log(`\n${'█'.repeat(100)}`);
  console.log(`🔄 RETESTING: ${q.id} - ${q.category.toUpperCase()}`);
  console.log(`Question: ${q.question}`);
  console.log(`Options: ${q.options.join(', ')}`);
  console.log(`Correct Answer: ${q.correct_answer}`);
  console.log(`Is Corrupted: ${q.is_corrupted || false}`);
  console.log(`${'█'.repeat(100)}`);
  
  // Log the exact query being sent
  console.log('\n📤 QUERY BEING SENT TO MODEL:');
  console.log(q.query_sent);
  console.log(`\n${'-'.repeat(80)}`);
  
  const startTime = Date.now();
  
  try {
    console.log('🚀 Sending request to GLM-4.5...');
    const response = await workflowEngine.answer(q.query_sent, {});
    const duration = Date.now() - startTime;
    
    console.log(`✅ Response received in ${duration}ms`);
    console.log(`📏 Response length: ${response.content.length} characters`);
    
    // Manual extraction with full logging
    const extraction = manualExtractAnswer(response.content, q.id, q.correct_answer);
    
    const result = {
      id: q.id,
      category: q.category,
      question: q.question.substring(0, 100) + '...',
      options: q.options,
      correct_answer: q.correct_answer,
      selected_answer: extraction.extracted,
      is_correct: extraction.is_correct,
      duration_ms: duration,
      extraction_pattern: extraction.pattern,
      all_matches: extraction.all_matches,
      is_corrupted: q.is_corrupted || false,
      full_response: response.content,
      retest: true,
      timestamp: new Date().toISOString()
    };
    
    console.log(`\n🏁 FINAL RESULT: ${extraction.extracted} (${extraction.is_correct ? 'CORRECT' : 'WRONG'})`);
    
    return result;
    
  } catch (error) {
    console.error(`\n❌ ERROR testing ${q.id}:`, error.message);
    console.error('Full error:', error);
    
    return {
      id: q.id,
      category: q.category,
      correct_answer: q.correct_answer,
      error: error.message,
      error_stack: error.stack,
      duration_ms: Date.now() - startTime,
      retest: true,
      timestamp: new Date().toISOString()
    };
  }
}

async function runCompleteRetest() {
  console.log('🧪 COMPLETE RETEST OF ALL FAILED GLM-4.5 QUESTIONS');
  console.log('Model:', process.env.OPENROUTER_MODEL || 'z-ai/glm-4.5');
  console.log('Date:', new Date().toISOString());
  
  const allQuestions = getAllQuestions();
  console.log(`\nTotal questions available: ${allQuestions.length}`);
  
  // Load previous results to identify failures
  const previousResults = JSON.parse(fs.readFileSync('./glm-4.5-fixed-final.json', 'utf8'));
  
  // Find all failed questions (extraction failed, errors, corrupted)
  const failedIds = previousResults.results
    .filter(r => !r.selected_answer || r.error || r.selected_answer === null)
    .map(r => r.id);
  
  // Add the corrupted question that was skipped
  const corruptedQuestion = allQuestions.find(q => q.is_corrupted);
  if (corruptedQuestion && !failedIds.includes(corruptedQuestion.id)) {
    failedIds.push(corruptedQuestion.id);
  }
  
  console.log(`\nFailed questions to retest: ${failedIds.join(', ')}`);
  console.log(`Total to retest: ${failedIds.length}`);
  
  const questionsToRetest = allQuestions.filter(q => failedIds.includes(q.id));
  
  console.log('\n📋 QUESTIONS TO RETEST:');
  questionsToRetest.forEach((q, i) => {
    console.log(`${i + 1}. ${q.id}: ${q.category} - ${q.question.substring(0, 60)}...`);
  });
  
  const retestResults = [];
  
  for (let i = 0; i < questionsToRetest.length; i++) {
    const q = questionsToRetest[i];
    console.log(`\n\n>>> TESTING ${i + 1}/${questionsToRetest.length} <<<`);
    
    const result = await retestFailedQuestion(q);
    retestResults.push(result);
    
    // Save incremental results
    const incrementalSave = {
      model: 'z-ai/glm-4.5',
      retest_timestamp: new Date().toISOString(),
      progress: `${i + 1}/${questionsToRetest.length}`,
      retest_results: retestResults,
      current_stats: {
        tested: retestResults.length,
        successful_extractions: retestResults.filter(r => r.selected_answer && !r.error).length,
        correct: retestResults.filter(r => r.is_correct).length,
        errors: retestResults.filter(r => r.error).length
      }
    };
    
    fs.writeFileSync('glm-4.5-complete-retest.json', JSON.stringify(incrementalSave, null, 2));
    
    // Delay between questions
    if (i < questionsToRetest.length - 1) {
      console.log('\n⏳ Waiting 3 seconds before next question...');
      await new Promise(resolve => setTimeout(resolve, 3000));
    }
  }
  
  // Final analysis
  const successful = retestResults.filter(r => r.selected_answer && !r.error);
  const correct = retestResults.filter(r => r.is_correct);
  const errors = retestResults.filter(r => r.error);
  
  console.log('\n' + '█'.repeat(100));
  console.log('🏁 COMPLETE RETEST RESULTS');
  console.log('█'.repeat(100));
  console.log(`\n✅ Successful extractions: ${successful.length}/${retestResults.length}`);
  console.log(`🎯 Correct answers: ${correct.length}/${successful.length} (${Math.round(correct.length/successful.length * 100)}%)`);
  console.log(`❌ Errors: ${errors.length}`);
  
  // Now calculate the true 30/30 score
  const originalSuccessful = previousResults.results.filter(r => r.selected_answer && !r.error && r.is_correct);
  const newCorrect = correct.length;
  const totalCorrect = originalSuccessful.length + newCorrect;
  
  console.log('\n📊 FINAL 30/30 SCORE:');
  console.log(`Previously correct: ${originalSuccessful.length}`);
  console.log(`Newly correct: ${newCorrect}`);
  console.log(`TOTAL CORRECT: ${totalCorrect}/30 (${Math.round(totalCorrect/30 * 100)}%)`);
  
  // Save final results
  const finalResults = {
    model: 'z-ai/glm-4.5',
    test_date: '2025-07-31',
    retest_completed: true,
    total_questions: 30,
    retest_results: retestResults,
    final_score: {
      correct: totalCorrect,
      total: 30,
      percentage: Math.round(totalCorrect/30 * 100)
    },
    comparison: {
      'o4-mini-high': '29/30 (96.7%)',
      'GLM-4.5-air': '24/30 (80.0%)',
      'z-ai/glm-4.5': `${totalCorrect}/30 (${Math.round(totalCorrect/30 * 100)}%)`
    }
  };
  
  fs.writeFileSync('glm-4.5-complete-final.json', JSON.stringify(finalResults, null, 2));
  console.log(`\n💾 Complete results saved to: glm-4.5-complete-final.json`);
  
  return finalResults;
}

runCompleteRetest().catch(console.error);