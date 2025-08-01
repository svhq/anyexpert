// Verify the 3 questions with answer extraction issues
const fs = require('fs');
const workflowEngine = require('./src/workflow-engine');

// Load the batch test data to get the actual questions
const batchData = JSON.parse(fs.readFileSync('./mmlu_batch_test_2025-07-31T02-07-13-672Z.json', 'utf8'));

// Questions to verify
const questionsToVerify = [
  {
    name: 'Business Question (Original Test #3)',
    question: 'Typical advertising regulatory bodies suggest, for example that adverts must not: encourage _________, cause unnecessary ________ or _____, and must not cause _______ offence.',
    options: [
      'A) Safe practices, Fear, Jealousy, Trivial',
      'B) Unsafe practices, Distress, Joy, Trivial',
      'C) Safe practices, Wants, Jealousy, Trivial',
      'D) Safe practices, Distress, Fear, Trivial',
      'E) Unsafe practices, Wants, Jealousy, Serious',
      'F) Safe practices, Distress, Jealousy, Serious',
      'G) Safe practices, Wants, Fear, Serious',
      'H) Unsafe practices, Wants, Fear, Trivial',
      'I) Unsafe practices, Distress, Fear, Serious'
    ],
    correct_answer: 'I',
    query_sent: 'Typical advertising regulatory bodies suggest, for example that adverts must not: encourage _________, cause unnecessary ________ or _____, and must not cause _______ offence.\n\nOptions:\nA) Safe practices, Fear, Jealousy, Trivial\nB) Unsafe practices, Distress, Joy, Trivial\nC) Safe practices, Wants, Jealousy, Trivial\nD) Safe practices, Distress, Fear, Trivial\nE) Unsafe practices, Wants, Jealousy, Serious\nF) Safe practices, Distress, Jealousy, Serious\nG) Safe practices, Wants, Fear, Serious\nH) Unsafe practices, Wants, Fear, Trivial\nI) Unsafe practices, Distress, Fear, Serious\n\nPlease select the correct answer and explain your reasoning.'
  },
  {
    name: 'Chemistry Question (Batch 1 #1)',
    question: 'Which of the following compounds would be expected to have the highest boiling point?',
    correct_answer: 'B',
    query_sent: 'Which of the following compounds would be expected to have the highest boiling point?\n\nOptions:\nA) CH₃OCH₃\nB) CH₃CH₂CH₂OH\nC) CH₃CH₂CH₃\nD) CH₃CH₂F\nE) CH₃CHO\n\nPlease select the correct answer and explain your reasoning.'
  },
  {
    name: 'Physics Question (Batch 1 #5)', 
    question: 'A particle of mass m moves in a one-dimensional potential V(x) = ½kx². In quantum mechanics, the energy levels of this harmonic oscillator are:',
    correct_answer: 'B',
    query_sent: 'A particle of mass m moves in a one-dimensional potential V(x) = ½kx². In quantum mechanics, the energy levels of this harmonic oscillator are:\n\nOptions:\nA) En = nℏω\nB) En = (n + ½)ℏω\nC) En = n²ℏω\nD) En = (2n + 1)ℏω\nE) En = nℏω/2\n\nPlease select the correct answer and explain your reasoning.'
  }
];

async function verifyQuestion(questionData) {
  console.log('\n' + '='.repeat(80));
  console.log(`\n🔍 Verifying: ${questionData.name}`);
  console.log('Correct Answer:', questionData.correct_answer);
  console.log('\nSending query to o4-mini-high...');
  
  try {
    const response = await workflowEngine.answer(questionData.query_sent, {});
    
    console.log('\n📤 FULL RESPONSE:');
    console.log(response.content);
    console.log('\n' + '-'.repeat(40));
    
    // Look for the answer in various formats
    const answerPatterns = [
      new RegExp(`\\b${questionData.correct_answer}\\)`, 'i'),
      new RegExp(`answer is.*${questionData.correct_answer}`, 'i'),
      new RegExp(`Option ${questionData.correct_answer}`, 'i'),
      new RegExp(`\\(${questionData.correct_answer}\\)`, 'i'),
      new RegExp(`\\*\\*${questionData.correct_answer}\\)`, 'i')
    ];
    
    let foundCorrectAnswer = false;
    for (const pattern of answerPatterns) {
      if (pattern.test(response.content)) {
        foundCorrectAnswer = true;
        break;
      }
    }
    
    console.log('\n📊 ANALYSIS:');
    console.log('- Contains correct answer letter:', foundCorrectAnswer ? '✅ YES' : '❌ NO');
    
    // For fill-in-the-blank questions, check if correct words are in right order
    if (questionData.name.includes('Business')) {
      const correctWords = ['Unsafe practices', 'Distress', 'Fear', 'Serious'];
      const hasAllWords = correctWords.every(word => 
        response.content.toLowerCase().includes(word.toLowerCase())
      );
      console.log('- Contains all correct words:', hasAllWords ? '✅ YES' : '❌ NO');
    }
    
    // Extract what the model thinks is the answer
    const finalAnswerMatch = response.content.match(/(?:final answer|answer is|correct answer is|Answer:)[\s:]*([A-J])/i);
    if (finalAnswerMatch) {
      console.log('- Model\'s selected answer:', finalAnswerMatch[1]);
      console.log('- Is correct:', finalAnswerMatch[1] === questionData.correct_answer ? '✅ YES' : '❌ NO');
    }
    
    return {
      question: questionData.name,
      correct_answer: questionData.correct_answer,
      found_correct_answer: foundCorrectAnswer,
      response_length: response.content.length
    };
    
  } catch (error) {
    console.error('❌ Error:', error.message);
    return { question: questionData.name, error: error.message };
  }
}

async function runVerification() {
  console.log('🧪 Verifying Answer Extraction Issues');
  console.log('Model:', process.env.OPENROUTER_MODEL);
  console.log('\nChecking 3 questions that had "answer detection" issues...');
  
  const results = [];
  for (const question of questionsToVerify) {
    const result = await verifyQuestion(question);
    results.push(result);
    
    // Wait between queries
    await new Promise(resolve => setTimeout(resolve, 2000));
  }
  
  console.log('\n' + '='.repeat(80));
  console.log('\n📋 VERIFICATION SUMMARY\n');
  results.forEach(r => {
    console.log(`${r.question}:`);
    console.log(`  - Correct answer: ${r.correct_answer}`);
    console.log(`  - Found in response: ${r.found_correct_answer ? '✅ YES' : '❌ NO'}`);
  });
}

runVerification().catch(console.error);