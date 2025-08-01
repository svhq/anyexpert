const fs = require('fs');
const workflowEngine = require('./src/workflow-engine');

// Load the cleaned GPQA questions
const data = JSON.parse(fs.readFileSync('gpqa-10-questions-clean.json', 'utf8'));

// Test a single question
async function testQuestion(index) {
  const question = data.questions[index];
  
  console.log(`\n${'='.repeat(70)}`);
  console.log(`RETESTING QUESTION ${index + 1} WITH E2B SERVICE`);
  console.log('='.repeat(70));
  console.log('\nQuestion ID:', question.id);
  console.log('Category:', question.category);
  console.log('Domain:', question.metadata.domain);
  console.log('Correct Answer:', question.correct_answer);
  console.log('\n' + '-'.repeat(70) + '\n');
  
  // Format question
  const prompt = `Question: ${question.question}

Options:

A) ${question.options.A}
B) ${question.options.B}
C) ${question.options.C}
D) ${question.options.D}`;

  console.log('SENDING TO MODEL WITH E2B AVAILABLE...\n');
  
  try {
    const startTime = Date.now();
    
    // Use 5 minute timeout
    const response = await workflowEngine.answer(prompt, { timeout: 300000 });
    const duration = Date.now() - startTime;
    
    const modelResponse = response.content;
    
    console.log('MODEL RESPONSE:');
    console.log('=' + '='.repeat(70));
    console.log(modelResponse);
    console.log('=' + '='.repeat(70));
    
    // Extract answer - enhanced patterns
    let extractedAnswer = null;
    const patterns = [
      /answer\s*(?:is|:)\s*([A-D])/i,
      /correct\s*answer\s*(?:is|:)\s*([A-D])/i,
      /Therefore,?\s*([A-D])/i,
      /choose\s*([A-D])/i,
      /select\s*([A-D])/i,
      /\b([A-D])\)\s*is\s*(?:the\s*)?correct/i,
      /The answer is\s*([A-D])/i,
      /Option\s*([A-D])\s*is\s*correct/i,
      /correct\s*option\s*is\s*([A-D])/i
    ];
    
    for (const pattern of patterns) {
      const match = modelResponse.match(pattern);
      if (match) {
        extractedAnswer = match[1].toUpperCase();
        break;
      }
    }
    
    // Also check for bold/emphasized answers
    if (!extractedAnswer) {
      const boldPattern = /\*\*\s*([A-D])\s*\*\*|Answer:\s*\*\*\s*([A-D])\s*\*\*/i;
      const boldMatch = modelResponse.match(boldPattern);
      if (boldMatch) {
        extractedAnswer = (boldMatch[1] || boldMatch[2]).toUpperCase();
      }
    }
    
    const isCorrect = extractedAnswer === question.correct_answer;
    
    console.log('\nRESULTS:');
    console.log('-'.repeat(70));
    console.log('Extracted Answer:', extractedAnswer || 'Could not extract');
    console.log('Correct Answer:', question.correct_answer);
    console.log('Is Correct:', isCorrect ? '✅ YES' : '❌ NO');
    console.log('Response Time:', `${(duration / 1000).toFixed(2)} seconds`);
    console.log('Search Used:', response.searchPerformed ? 'Yes' : 'No');
    console.log('Code Executed:', response.codeExecuted ? 'Yes' : 'No');
    
    // Save result
    const result = {
      test_date: new Date().toISOString(),
      model: 'z-ai/glm-4.5-air:free',
      question_id: question.id,
      question_category: question.category,
      question_domain: question.metadata.domain,
      e2b_service: 'running',
      prompt_sent: prompt,
      model_response: modelResponse,
      extracted_answer: extractedAnswer,
      correct_answer: question.correct_answer,
      is_correct: isCorrect,
      response_time_ms: duration,
      search_used: response.searchPerformed,
      code_executed: response.codeExecuted
    };
    
    fs.writeFileSync(
      `gpqa-glm45air-q${index+1}-retest.json`, 
      JSON.stringify(result, null, 2)
    );
    
    console.log(`\n✓ Result saved to gpqa-glm45air-q${index+1}-retest.json`);
    
    return result;
    
  } catch (error) {
    console.error('Error testing question:', error.message);
    return {
      question_id: question.id,
      error: error.message,
      timestamp: new Date().toISOString()
    };
  }
}

// Main function to test questions 2 and 3
async function main() {
  console.log('RETESTING GPQA QUESTIONS 2 & 3 WITH E2B SERVICE\n');
  
  // Test question 2 (Physics)
  const result2 = await testQuestion(1);
  
  console.log('\nPausing 3 seconds...');
  await new Promise(resolve => setTimeout(resolve, 3000));
  
  // Test question 3 (Chemistry)
  const result3 = await testQuestion(2);
  
  console.log('\n' + '='.repeat(70));
  console.log('RETESTING COMPLETE');
  console.log('='.repeat(70));
}

main().catch(console.error);