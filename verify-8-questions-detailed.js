// Detailed verification of 8 questionable GLM-4.5 responses
const fs = require('fs');
const workflowEngine = require('./src/workflow-engine');

// The 8 questions that need detailed verification
const questionableResults = [
  'batch1-2', // Business - said J but extracted null
  'batch1-3', // Law - said A but correct is D  
  'batch1-5', // Psychology - said E but extracted A
  'batch2-1', // Psychology - said H but correct is I (corrupted?)
  'batch2-3', // Biology - said D but extracted A
  'new-3',    // Law - said C but correct is D
  'new-6',    // Psychology - said H but correct is I (corrupted?)
  'batch2-4'  // Missing from retest - need to test
];

// Get all questions with full details
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
  
  // New 10 questions
  const newQuestions = JSON.parse(fs.readFileSync('./mmlu-10-new-questions.json', 'utf8'));
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
  
  // Batch questions - reconstruct from new questions
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
  
  return allQuestions;
}

async function verifyQuestionDetailed(questionId, questionData) {
  console.log('\n' + '█'.repeat(120));
  console.log(`🔍 DETAILED VERIFICATION: ${questionId}`);
  console.log('█'.repeat(120));
  
  console.log(`📋 QUESTION DETAILS:`);
  console.log(`- ID: ${questionData.id}`);
  console.log(`- Category: ${questionData.category}`);
  console.log(`- Correct Answer: ${questionData.correct_answer}`);
  console.log(`- Is Corrupted: ${questionData.is_corrupted || false}`);
  
  console.log(`\n📝 QUESTION TEXT:`);
  console.log(questionData.question);
  
  console.log(`\n📋 OPTIONS PROVIDED:`);
  questionData.options.forEach((option, i) => {
    console.log(`${String.fromCharCode(65 + i)}) ${option.replace(/^[A-J]\)\s*/, '')}`);
  });
  
  console.log(`\n📤 EXACT QUERY SENT TO MODEL:`);
  console.log(`"${questionData.query_sent}"`);
  
  console.log(`\n🚀 SENDING TO MODEL...`);
  const startTime = Date.now();
  
  try {
    const response = await workflowEngine.answer(questionData.query_sent, {});
    const duration = Date.now() - startTime;
    
    console.log(`\n✅ RESPONSE RECEIVED (${duration}ms)`);
    console.log(`📏 Response Length: ${response.content.length} characters`);
    
    console.log(`\n📄 FULL MODEL RESPONSE:`);
    console.log('─'.repeat(80));
    console.log(response.content);
    console.log('─'.repeat(80));
    
    // Multiple extraction attempts
    console.log(`\n🔍 EXTRACTION ANALYSIS:`);
    
    // Pattern 1: Direct answer statements
    console.log(`\n1️⃣ LOOKING FOR DIRECT ANSWER STATEMENTS:`);
    const directPatterns = [
      /(?:correct answer is|answer is|the answer is)[:\s]*\*?\*?([A-J])\)/gi,
      /therefore[,\s]*(?:the )?(?:correct )?answer is[:\s]*\*?\*?([A-J])\)/gi,
      /based on (?:this )?analysis[,\s]*(?:the )?(?:correct )?answer is[:\s]*\*?\*?([A-J])\)/gi
    ];
    
    let directAnswers = [];
    directPatterns.forEach((pattern, i) => {
      const matches = [...response.content.matchAll(pattern)];
      matches.forEach(match => {
        directAnswers.push({
          pattern: `Direct Pattern ${i + 1}`,
          letter: match[1],
          context: response.content.substring(Math.max(0, match.index - 30), match.index + 70)
        });
      });
    });
    
    if (directAnswers.length > 0) {
      console.log(`✅ Found ${directAnswers.length} direct answer(s):`);
      directAnswers.forEach(ans => {
        console.log(`   - ${ans.letter} via ${ans.pattern}: "${ans.context.trim()}"`);
      });
    } else {
      console.log(`❌ No direct answer statements found`);
    }
    
    // Pattern 2: Bold/emphasized answers
    console.log(`\n2️⃣ LOOKING FOR BOLD/EMPHASIZED ANSWERS:`);
    const boldPatterns = [
      /\*\*([A-J])\)\*\*/g,
      /\*\*.*?([A-J])\).*?\*\*/g,
      /\*([A-J])\)\*/g
    ];
    
    let boldAnswers = [];
    boldPatterns.forEach((pattern, i) => {
      const matches = [...response.content.matchAll(pattern)];
      matches.forEach(match => {
        boldAnswers.push({
          pattern: `Bold Pattern ${i + 1}`,
          letter: match[1],
          context: response.content.substring(Math.max(0, match.index - 30), match.index + 70)
        });
      });
    });
    
    if (boldAnswers.length > 0) {
      console.log(`✅ Found ${boldAnswers.length} bold answer(s):`);
      boldAnswers.forEach(ans => {
        console.log(`   - ${ans.letter} via ${ans.pattern}: "${ans.context.trim()}"`);
      });
    } else {
      console.log(`❌ No bold answers found`);
    }
    
    // Pattern 3: Option analysis
    console.log(`\n3️⃣ LOOKING FOR OPTION SELECTION:`);
    const optionPatterns = [
      /option ([A-J])[^a-zA-Z]*(?:is )?(?:the )?correct/gi,
      /select[:\s]*option[:\s]*([A-J])/gi,
      /choose[:\s]*option[:\s]*([A-J])/gi,
      /pick[:\s]*option[:\s]*([A-J])/gi
    ];
    
    let optionAnswers = [];
    optionPatterns.forEach((pattern, i) => {
      const matches = [...response.content.matchAll(pattern)];
      matches.forEach(match => {
        optionAnswers.push({
          pattern: `Option Pattern ${i + 1}`,
          letter: match[1],
          context: response.content.substring(Math.max(0, match.index - 30), match.index + 70)
        });
      });
    });
    
    if (optionAnswers.length > 0) {
      console.log(`✅ Found ${optionAnswers.length} option selection(s):`);
      optionAnswers.forEach(ans => {
        console.log(`   - ${ans.letter} via ${ans.pattern}: "${ans.context.trim()}"`);
      });
    } else {
      console.log(`❌ No option selections found`);
    }
    
    // Pattern 4: Final conclusion scanning
    console.log(`\n4️⃣ SCANNING FINAL CONCLUSIONS:`);
    const conclusionSections = response.content.split(/\n\n/).slice(-3); // Last 3 paragraphs
    let conclusionAnswers = [];
    
    conclusionSections.forEach((section, i) => {
      const letterMatches = section.match(/\b([A-J])\)/g);
      if (letterMatches) {
        letterMatches.forEach(match => {
          conclusionAnswers.push({
            section: `Final Section ${i + 1}`,
            letter: match[0],
            context: section.substring(0, 100) + '...'
          });
        });
      }
    });
    
    if (conclusionAnswers.length > 0) {
      console.log(`✅ Found letters in conclusions:`);
      conclusionAnswers.forEach(ans => {
        console.log(`   - ${ans.letter} in ${ans.section}: "${ans.context}"`);
      });
    } else {
      console.log(`❌ No letters found in final conclusions`);
    }
    
    // Determine final answer
    console.log(`\n🎯 FINAL DETERMINATION:`);
    
    let finalAnswer = null;
    let confidence = 'LOW';
    let method = 'NONE';
    
    // Prioritize direct answers
    if (directAnswers.length > 0) {
      const freq = {};
      directAnswers.forEach(ans => freq[ans.letter] = (freq[ans.letter] || 0) + 1);
      const mostFreq = Object.entries(freq).sort((a, b) => b[1] - a[1])[0];
      finalAnswer = mostFreq[0];
      confidence = directAnswers.length > 1 ? 'HIGH' : 'MEDIUM';
      method = 'DIRECT_STATEMENT';
    }
    // Then bold answers
    else if (boldAnswers.length > 0) {
      const freq = {};
      boldAnswers.forEach(ans => freq[ans.letter] = (freq[ans.letter] || 0) + 1);
      const mostFreq = Object.entries(freq).sort((a, b) => b[1] - a[1])[0];
      finalAnswer = mostFreq[0];
      confidence = 'MEDIUM';
      method = 'BOLD_FORMAT';
    }
    // Then option selections
    else if (optionAnswers.length > 0) {
      const freq = {};
      optionAnswers.forEach(ans => freq[ans.letter] = (freq[ans.letter] || 0) + 1);
      const mostFreq = Object.entries(freq).sort((a, b) => b[1] - a[1])[0];
      finalAnswer = mostFreq[0];
      confidence = 'LOW';
      method = 'OPTION_MENTION';
    }
    
    console.log(`📊 EXTRACTION RESULT:`);
    console.log(`   - Final Answer: ${finalAnswer || 'NULL'}`);
    console.log(`   - Confidence: ${confidence}`);
    console.log(`   - Method: ${method}`);
    console.log(`   - Correct Answer: ${questionData.correct_answer}`);
    console.log(`   - Is Correct: ${finalAnswer === questionData.correct_answer ? '✅ YES' : '❌ NO'}`);
    
    // Special check for corrupted questions
    if (questionData.is_corrupted) {
      console.log(`   - ⚠️ NOTE: This question is corrupted (option ${questionData.correct_answer} doesn't exist)`);
    }
    
    return {
      id: questionData.id,
      category: questionData.category,
      question_received: questionData.question,
      options_received: questionData.options, 
      correct_answer: questionData.correct_answer,
      model_response: response.content,
      extracted_answer: finalAnswer,
      extraction_method: method,
      extraction_confidence: confidence,
      is_correct: finalAnswer === questionData.correct_answer,
      is_corrupted: questionData.is_corrupted || false,
      response_length: response.content.length,
      duration_ms: duration,
      timestamp: new Date().toISOString()
    };
    
  } catch (error) {
    console.log(`\n❌ ERROR: ${error.message}`);
    return {
      id: questionData.id,
      error: error.message,
      timestamp: new Date().toISOString()
    };
  }
}

async function runDetailedVerification() {
  console.log('🔍 DETAILED VERIFICATION OF 8 QUESTIONABLE GLM-4.5 RESULTS');
  console.log('Model:', process.env.OPENROUTER_MODEL || 'z-ai/glm-4.5');
  console.log('Date:', new Date().toISOString());
  console.log('\n' + '═'.repeat(120));
  
  const allQuestions = getAllQuestions();
  const questionsToTest = [];
  
  // Find the 8 questions
  questionableResults.forEach(qId => {
    const question = allQuestions.find(q => q.id === qId);
    if (question) {
      questionsToTest.push(question);
    } else {
      console.log(`⚠️ Question ${qId} not found in dataset`);
    }
  });
  
  console.log(`\nFound ${questionsToTest.length} questions to verify:\n`);
  questionsToTest.forEach((q, i) => {
    console.log(`${i + 1}. ${q.id}: ${q.category} - Correct: ${q.correct_answer}`);
  });
  
  const verificationResults = [];
  
  // Test each question with detailed logging
  for (let i = 0; i < questionsToTest.length; i++) {
    const question = questionsToTest[i];
    
    console.log(`\n\n>>> TESTING ${i + 1}/${questionsToTest.length} <<<`);
    
    const result = await verifyQuestionDetailed(question.id, question);
    verificationResults.push(result);
    
    // Save incremental results
    const summary = {
      model: 'z-ai/glm-4.5',
      verification_type: 'detailed_8_questions',
      timestamp: new Date().toISOString(),
      progress: `${i + 1}/${questionsToTest.length}`,
      results: verificationResults
    };
    fs.writeFileSync('detailed-verification-8.json', JSON.stringify(summary, null, 2));
    
    // Wait between questions
    if (i < questionsToTest.length - 1) {
      console.log('\n⏳ Waiting 5 seconds before next question...');
      await new Promise(resolve => setTimeout(resolve, 5000));
    }
  }
  
  // Final summary
  console.log('\n\n' + '█'.repeat(120));
  console.log('📊 DETAILED VERIFICATION SUMMARY');
  console.log('█'.repeat(120));
  
  const successful = verificationResults.filter(r => !r.error);
  const correct = verificationResults.filter(r => r.is_correct);
  const errors = verificationResults.filter(r => r.error);
  const corrupted = verificationResults.filter(r => r.is_corrupted);
  
  console.log(`\n📈 RESULTS:`);
  console.log(`- Questions tested: ${verificationResults.length}`);
  console.log(`- Successful responses: ${successful.length}`);
  console.log(`- Correct answers: ${correct.length}`);
  console.log(`- Errors: ${errors.length}`);
  console.log(`- Corrupted questions: ${corrupted.length}`);
  
  console.log(`\n📋 DETAILED BREAKDOWN:`);
  verificationResults.forEach(result => {
    if (result.error) {
      console.log(`❌ ${result.id}: ERROR - ${result.error}`);
    } else {
      const status = result.is_correct ? '✅' : '❌';
      const corrupt = result.is_corrupted ? ' (CORRUPTED)' : '';
      console.log(`${status} ${result.id}: ${result.extracted_answer} → ${result.correct_answer}${corrupt}`);
    }
  });
  
  // Save final results
  const finalSummary = {
    model: 'z-ai/glm-4.5',
    verification_type: 'detailed_8_questions',
    completed_timestamp: new Date().toISOString(),
    summary: {
      total_tested: verificationResults.length,
      successful_responses: successful.length,
      correct_answers: correct.length,
      errors: errors.length,
      corrupted_questions: corrupted.length,
      accuracy_on_valid: Math.round(correct.filter(r => !r.is_corrupted).length / successful.filter(r => !r.is_corrupted).length * 100)
    },
    results: verificationResults
  };
  
  fs.writeFileSync('detailed-verification-8-final.json', JSON.stringify(finalSummary, null, 2));
  console.log(`\n💾 Detailed verification saved to: detailed-verification-8-final.json`);
  
  return finalSummary;
}

runDetailedVerification().catch(console.error);