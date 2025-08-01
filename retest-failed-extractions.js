// Retest questions with failed extractions
const fs = require('fs');
const workflowEngine = require('./src/workflow-engine');

// Load failed results
const results = JSON.parse(fs.readFileSync('./glm-4.5-fixed-results.json', 'utf8'));

// Get all original questions for reference
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
  
  // Batch questions
  const batch = JSON.parse(fs.readFileSync('./mmlu_batch_test_2025-07-31T02-07-13-672Z.json', 'utf8'));
  
  // Batch 1
  batch.batch1.results.forEach((r, idx) => {
    const optionsMatch = r.response.content_preview.match(/Options:[\s\S]*?(?=\n\nPlease)/);
    const options = [];
    if (optionsMatch) {
      const optText = optionsMatch[0].replace('Options:', '').trim();
      const lines = optText.split('\n').filter(l => l.match(/^[A-J]\)/));
      options.push(...lines);
    }
    
    allQuestions.push({
      id: `batch1-${idx + 1}`,
      category: r.category,
      question: r.question,
      options: options,
      correct_answer: r.correct_answer,
      query_sent: `${r.question}\n\nOptions:\n${options.join('\n')}\n\nPlease select the correct answer and explain your reasoning.`
    });
  });
  
  // Batch 2
  batch.batch2.results.forEach((r, idx) => {
    const optionsMatch = r.response.content_preview.match(/Options:[\s\S]*?(?=\n\nPlease)/);
    const options = [];
    if (optionsMatch) {
      const optText = optionsMatch[0].replace('Options:', '').trim();
      const lines = optText.split('\n').filter(l => l.match(/^[A-J]\)/));
      options.push(...lines);
    }
    
    allQuestions.push({
      id: `batch2-${idx + 1}`,
      category: r.category,
      question: r.question,
      options: options,
      correct_answer: r.correct_answer,
      query_sent: `${r.question}\n\nOptions:\n${options.join('\n')}\n\nPlease select the correct answer and explain your reasoning.`
    });
  });
  
  // New 10 questions
  const newQuestions = JSON.parse(fs.readFileSync('./mmlu-10-new-questions.json', 'utf8'));
  newQuestions.questions.forEach((q, idx) => {
    if (q.correct_answer === 'I' && !q.options.find(o => o.startsWith('I'))) {
      return; // Skip corrupted
    }
    allQuestions.push({
      id: `new-${idx + 1}`,
      category: q.category,
      question: q.question,
      options: q.options,
      correct_answer: q.correct_answer,
      query_sent: `${q.question}\n\nOptions:\n${q.options.join('\n')}\n\nPlease select the correct answer and explain your reasoning.`
    });
  });
  
  return allQuestions;
}

// Enhanced extraction with debugging
async function extractAnswerWithDebug(responseContent, questionId) {
  console.log(`\n🔍 Debugging extraction for ${questionId}:`);
  console.log(`Response length: ${responseContent.length} chars`);
  
  let selectedAnswer = null;
  let matchedPattern = 'none';
  
  // Pattern 1: **X)** format (GLM-4.5 common format)
  let match = responseContent.match(/\*\*([A-J])\).*?\*\*/i);
  if (match) {
    selectedAnswer = match[1];
    matchedPattern = '**X)**';
  }
  
  // Pattern 2: **X)** without content after
  if (!match) {
    match = responseContent.match(/\*\*([A-J])\)\*\*/i);
    if (match) {
      selectedAnswer = match[1];
      matchedPattern = '**X)**';
    }
  }
  
  // Pattern 3: "correct answer is X)" with or without bold
  if (!match) {
    match = responseContent.match(/correct answer is[\s:]*\*?\*?([A-J])\)/i);
    if (match) {
      selectedAnswer = match[1];
      matchedPattern = 'correct answer is X)';
    }
  }
  
  // Pattern 4: "Answer: X)" or "Answer is X)"
  if (!match) {
    match = responseContent.match(/Answer[\s:]+(?:is[\s:]+)?\*?\*?([A-J])\)/i);
    if (match) {
      selectedAnswer = match[1];
      matchedPattern = 'Answer: X)';
    }
  }
  
  // Pattern 5: "The answer is X"
  if (!match) {
    match = responseContent.match(/The answer is[:\s]*\*?\*?([A-J])\)/i);
    if (match) {
      selectedAnswer = match[1];
      matchedPattern = 'The answer is X';
    }
  }
  
  // Pattern 6: Look for bold option letters
  if (!match) {
    match = responseContent.match(/\*\*.*?([A-J])\).*?\*\*/i);
    if (match) {
      selectedAnswer = match[1];
      matchedPattern = '**...X)...**';
    }
  }
  
  // Pattern 7: "select X)" or "choose X)"
  if (!match) {
    match = responseContent.match(/(?:select|choose)[\s:]*([A-J])\)/i);
    if (match) {
      selectedAnswer = match[1];
      matchedPattern = 'select/choose X)';
    }
  }
  
  // Pattern 8: Last resort - any option letter followed by )
  if (!match) {
    const allMatches = responseContent.match(/\b([A-J])\)/g);
    if (allMatches && allMatches.length > 0) {
      // Take the most frequent option, or last one
      const freq = {};
      allMatches.forEach(m => {
        const letter = m.charAt(0);
        freq[letter] = (freq[letter] || 0) + 1;
      });
      
      // Get most frequent
      let maxCount = 0;
      let mostFrequent = null;
      Object.entries(freq).forEach(([letter, count]) => {
        if (count > maxCount) {
          maxCount = count;
          mostFrequent = letter;
        }
      });
      
      selectedAnswer = mostFrequent;
      matchedPattern = 'most frequent option';
    }
  }
  
  console.log(`Extracted: ${selectedAnswer || 'NULL'} (pattern: ${matchedPattern})`);
  
  // Show response preview for debugging
  if (!selectedAnswer) {
    console.log(`Response preview: ${responseContent.substring(0, 500)}...`);
  }
  
  return selectedAnswer;
}

async function retestQuestion(q) {
  console.log(`\n🔄 Retesting: ${q.id} - ${q.category}`);
  console.log(`Question: ${q.question.substring(0, 80)}...`);
  
  const startTime = Date.now();
  
  try {
    const response = await workflowEngine.answer(q.query_sent, {});
    const duration = Date.now() - startTime;
    
    const selectedAnswer = await extractAnswerWithDebug(response.content, q.id);
    
    const result = {
      id: q.id,
      category: q.category,
      correct_answer: q.correct_answer,
      selected_answer: selectedAnswer,
      is_correct: selectedAnswer === q.correct_answer,
      duration_ms: duration,
      retest: true
    };
    
    console.log(`Result: ${selectedAnswer} (correct: ${q.correct_answer}) - ${result.is_correct ? '✅' : '❌'}`);
    
    return result;
    
  } catch (error) {
    console.error(`❌ Error retesting ${q.id}:`, error.message);
    return {
      id: q.id,
      category: q.category,
      correct_answer: q.correct_answer,
      error: error.message,
      retest: true
    };
  }
}

async function runRetests() {
  console.log('🧪 Retesting Failed Extractions for GLM-4.5');
  console.log('Model:', process.env.OPENROUTER_MODEL);
  
  const allQuestions = getAllQuestions();
  const questionMap = {};
  allQuestions.forEach(q => questionMap[q.id] = q);
  
  // Find failed extractions
  const failed = results.results.filter(r => 
    r.selected_answer === null || r.error || r.selected_answer === undefined
  );
  
  console.log(`\nFound ${failed.length} failed extractions to retest:`);
  failed.forEach(f => {
    console.log(`- ${f.id}: ${f.category} (${f.selected_answer || 'error'})`);
  });
  
  const retestResults = [];
  
  for (const failedResult of failed) {
    const question = questionMap[failedResult.id];
    if (!question) {
      console.log(`⚠️ Question ${failedResult.id} not found`);
      continue;
    }
    
    const result = await retestQuestion(question);
    retestResults.push(result);
    
    // Small delay
    await new Promise(resolve => setTimeout(resolve, 2000));
  }
  
  // Update original results with successful retests
  const updatedResults = [...results.results];
  retestResults.forEach(retest => {
    if (retest.selected_answer && !retest.error) {
      const originalIndex = updatedResults.findIndex(r => r.id === retest.id);
      if (originalIndex >= 0) {
        updatedResults[originalIndex] = {
          ...updatedResults[originalIndex],
          selected_answer: retest.selected_answer,
          is_correct: retest.is_correct,
          retest: true
        };
      }
    }
  });
  
  // Calculate new stats
  const correct = updatedResults.filter(r => r.is_correct).length;
  const errors = updatedResults.filter(r => r.error).length;
  const total = updatedResults.length - errors;
  
  console.log('\n' + '='.repeat(60));
  console.log('\n🏁 RETEST RESULTS');
  console.log(`✅ Successful retests: ${retestResults.filter(r => r.selected_answer && !r.error).length}`);
  console.log(`❌ Still failed: ${retestResults.filter(r => !r.selected_answer || r.error).length}`);
  
  console.log('\n📊 UPDATED STATS:');
  console.log(`✅ Total Correct: ${correct}/${total} (${Math.round(correct/total * 100)}%)`);
  console.log(`❌ Incorrect: ${total - correct}`);
  console.log(`⚠️ Errors: ${errors}`);
  
  // Save updated results
  const finalResults = {
    model: 'z-ai/glm-4.5',
    timestamp: new Date().toISOString(),
    results: updatedResults,
    summary: {
      total: total,
      correct: correct,
      incorrect: total - correct,
      errors: errors,
      accuracy: Math.round(correct / total * 100)
    },
    retest_summary: {
      failed_original: failed.length,
      retested: retestResults.length,
      successful_retests: retestResults.filter(r => r.selected_answer && !r.error).length
    }
  };
  
  fs.writeFileSync('glm-4.5-fixed-final.json', JSON.stringify(finalResults, null, 2));
  console.log(`\n💾 Updated results saved to: glm-4.5-fixed-final.json`);
  
  return finalResults;
}

runRetests().catch(console.error);