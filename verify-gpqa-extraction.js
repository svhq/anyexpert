const fs = require('fs');
const csv = require('csv-parse/sync');

// Load the extracted questions
const extractedQuestions = JSON.parse(fs.readFileSync('gpqa-10-questions.json', 'utf-8'));

// Load and parse the CSV
const csvContent = fs.readFileSync('GPQA dataset/gpqa_main.csv', 'utf-8');
const records = csv.parse(csvContent, {
  columns: true,
  skip_empty_lines: true
});

console.log('=== VERIFYING GPQA QUESTION EXTRACTION ===\n');
console.log(`Total records in CSV: ${records.length}`);
console.log(`Extracted questions: ${extractedQuestions.questions.length}\n`);

// Function to clean text for comparison
function cleanText(text) {
  if (!text) return '';
  return text
    .replace(/\s+/g, ' ')
    .replace(/[""]/g, '"')
    .replace(/'/g, "'")
    .trim();
}

// Check first 10 questions
for (let i = 0; i < Math.min(10, extractedQuestions.questions.length); i++) {
  const extracted = extractedQuestions.questions[i];
  const csvRow = records[i];
  
  console.log(`\n${'='.repeat(60)}`);
  console.log(`QUESTION ${i + 1} (${extracted.id})`);
  console.log(`${'='.repeat(60)}`);
  
  // Check question text
  const csvQuestion = cleanText(csvRow['Question']);
  const extractedQuestion = cleanText(extracted.question);
  const questionMatch = csvQuestion === extractedQuestion;
  
  console.log('\nQUESTION TEXT:');
  console.log(`Match: ${questionMatch ? '✓' : '✗'}`);
  if (!questionMatch) {
    console.log('CSV:', csvQuestion.substring(0, 100) + '...');
    console.log('Extracted:', extractedQuestion.substring(0, 100) + '...');
  }
  
  // Check correct answer
  const csvCorrect = cleanText(csvRow['Correct Answer']);
  const extractedCorrect = extracted.options[extracted.correct_answer];
  const correctMatch = csvCorrect === cleanText(extractedCorrect);
  
  console.log('\nCORRECT ANSWER:');
  console.log(`CSV: "${csvCorrect}"`);
  console.log(`Extracted: ${extracted.correct_answer} = "${extractedCorrect}"`);
  console.log(`Match: ${correctMatch ? '✓' : '✗'}`);
  
  // Check all options
  const csvOptions = [
    csvRow['Correct Answer'],
    csvRow['Incorrect Answer 1'],
    csvRow['Incorrect Answer 2'],
    csvRow['Incorrect Answer 3']
  ].map(cleanText).sort();
  
  const extractedOptions = Object.values(extracted.options).map(cleanText).sort();
  
  console.log('\nOPTIONS:');
  console.log('CSV Options:', csvOptions);
  console.log('Extracted Options:', extractedOptions);
  
  const allOptionsMatch = csvOptions.every((opt, idx) => opt === extractedOptions[idx]);
  console.log(`All options match: ${allOptionsMatch ? '✓' : '✗'}`);
  
  // Check explanation
  const csvExplanation = cleanText(csvRow['Explanation']);
  const extractedExplanation = cleanText(extracted.metadata.explanation);
  const explanationMatch = csvExplanation === extractedExplanation;
  
  console.log('\nEXPLANATION:');
  console.log(`Match: ${explanationMatch ? '✓' : '✗'}`);
  if (!explanationMatch) {
    console.log('Length difference:', Math.abs(csvExplanation.length - extractedExplanation.length));
  }
  
  // Overall status
  const allMatch = questionMatch && correctMatch && allOptionsMatch && explanationMatch;
  console.log(`\nOVERALL: ${allMatch ? '✓ ALL GOOD' : '✗ ISSUES FOUND'}`);
}