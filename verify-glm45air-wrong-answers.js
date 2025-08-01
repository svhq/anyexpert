// Verify GLM-4.5-air wrong answers
const fs = require('fs');

// Load the results and questions
const results = JSON.parse(fs.readFileSync('./glm-4.5-air-10-results-incremental.json', 'utf8'));
const newQuestions = JSON.parse(fs.readFileSync('./mmlu-10-new-questions.json', 'utf8'));

// Map results to questions
const wrongAnswers = [];
results.results.forEach((result, idx) => {
  if (!result.is_correct) {
    const question = newQuestions.questions[idx];
    wrongAnswers.push({
      question_number: result.question_number,
      category: result.category,
      question: question.question,
      options: question.options,
      correct_answer: result.correct_answer,
      selected_answer: result.selected_answer,
      question_id: question.question_id
    });
  }
});

console.log('GLM-4.5-AIR WRONG ANSWERS ANALYSIS');
console.log('='.repeat(80));
console.log(`Total wrong: ${wrongAnswers.length}/10`);
console.log('='.repeat(80));

wrongAnswers.forEach((item, idx) => {
  console.log(`\n${idx + 1}. Question ${item.question_number} (${item.category}) - ID: ${item.question_id}`);
  console.log('─'.repeat(80));
  console.log('QUESTION:', item.question);
  console.log('\nOPTIONS:');
  item.options.forEach(opt => console.log(opt));
  console.log(`\nCORRECT: ${item.correct_answer}`);
  console.log(`GLM SELECTED: ${item.selected_answer}`);
  console.log('─'.repeat(80));
});

// Save for reference
fs.writeFileSync('glm45air-wrong-answers.json', JSON.stringify(wrongAnswers, null, 2));
console.log('\n\nSaved to glm45air-wrong-answers.json');