// Find the corrupted psychology question from Excel
const XLSX = require('xlsx');
const fs = require('fs');

// Load Excel file
const workbook = XLSX.readFile('./MMLU-PRO Benchmark Questions.xlsx');
const sheet = workbook.Sheets[workbook.SheetNames[0]];
const data = XLSX.utils.sheet_to_json(sheet);

// Find question ID 1987
const targetId = 1987;
const question = data.find(row => row.question_id === targetId);

if (question) {
  console.log('Found question ID 1987:');
  console.log('Category:', question.category);
  console.log('Question:', question.question);
  console.log('Options:', question.options);
  console.log('Correct Answer:', question.correct_answer);
  console.log('Answer Index:', question.answer_index);
  console.log('\nParsing options...');
  
  // Try to parse the options
  const optionsStr = question.options;
  const matches = optionsStr.match(/'([^']+)'/g);
  if (matches) {
    const options = matches.map((m, i) => {
      const letter = String.fromCharCode(65 + i); // A, B, C, etc.
      return `${letter}) ${m.replace(/'/g, '').trim()}`;
    });
    
    console.log('\nFormatted options:');
    options.forEach(opt => console.log(opt));
    
    // Save the corrected question
    const correctedQuestion = {
      question_id: question.question_id,
      category: question.category,
      question: question.question,
      options: options,
      correct_answer: String.fromCharCode(65 + question.answer_index), // Convert index to letter
      answer_index: question.answer_index,
      original_correct_answer: question.correct_answer
    };
    
    fs.writeFileSync('corrected-question-1987.json', JSON.stringify(correctedQuestion, null, 2));
    console.log('\nCorrected question saved to corrected-question-1987.json');
  }
} else {
  console.log('Question ID 1987 not found');
}