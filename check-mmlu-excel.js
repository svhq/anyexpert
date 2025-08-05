const XLSX = require('xlsx');
const fs = require('fs');

// Read the Excel file
const workbook = XLSX.readFile('MMLU-PRO Benchmark Questions.xlsx');

// Get the first sheet name
const sheetName = workbook.SheetNames[0];
const sheet = workbook.Sheets[sheetName];

// Convert to JSON
const data = XLSX.utils.sheet_to_json(sheet);

console.log('Total questions:', data.length);
console.log('First few columns:', Object.keys(data[0] || {}));

// Search for our math question about f(g(x,y), g(f(x,y)))
const searchTerms = ['f(g(x,y), g(f(x,y)))', 'Let f(x,y) = x² + y²', 'x² + y²', 'x² - y²'];

console.log('\nSearching for question 21 or similar math questions...\n');

let found = false;
data.forEach((row, index) => {
  const questionText = JSON.stringify(row).toLowerCase();
  
  // Check if this might be our question
  if (questionText.includes('f(g(x,y)') || 
      questionText.includes('x² + y²') || 
      questionText.includes('x^2 + y^2') ||
      (questionText.includes('composite') && questionText.includes('function'))) {
    
    console.log(`\nRow ${index + 1}:`);
    console.log('Question:', row.question || row.Question || row.text || JSON.stringify(row).substring(0, 200));
    console.log('Answer:', row.answer || row.Answer || row.correct_answer || row.correct);
    console.log('Full row keys:', Object.keys(row));
    found = true;
  }
});

if (!found) {
  console.log('\nQuestion not found with initial search. Let me check the structure...');
  
  // Show first 5 rows to understand structure
  console.log('\nFirst 5 rows of data:');
  data.slice(0, 5).forEach((row, i) => {
    console.log(`\nRow ${i + 1}:`);
    Object.entries(row).forEach(([key, value]) => {
      if (typeof value === 'string' && value.length > 100) {
        console.log(`${key}: ${value.substring(0, 100)}...`);
      } else {
        console.log(`${key}: ${value}`);
      }
    });
  });
}

// Also check our JSON file for comparison
console.log('\n\nChecking mmlu-30-questions-formatted.json for question 21...');
const jsonData = JSON.parse(fs.readFileSync('mmlu-30-questions-formatted.json', 'utf-8'));
const q21 = jsonData.questions.find(q => q.id === 'mmlu-21');
if (q21) {
  console.log('\nQuestion 21 from JSON:');
  console.log('Question:', q21.question.substring(0, 150) + '...');
  console.log('Correct Answer:', q21.correct_answer);
}