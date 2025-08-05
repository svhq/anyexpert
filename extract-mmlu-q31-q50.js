const XLSX = require('xlsx');
const fs = require('fs');

// Read the Excel file
const workbook = XLSX.readFile('MMLU-PRO Benchmark Questions.xlsx');
const sheetName = workbook.SheetNames[0];
const sheet = workbook.Sheets[sheetName];
const data = XLSX.utils.sheet_to_json(sheet);

console.log('=== EXTRACTING MMLU QUESTIONS 31-50 ===\n');
console.log(`Total questions available: ${data.length}`);

// Extract questions 31-50 (indices 30-49)
const questions31to50 = data.slice(30, 50).map((row, index) => {
  const questionId = `mmlu-${31 + index}`;
  
  // Parse options - they might be in a single string or separate columns
  let options = row.options;
  if (typeof options === 'string') {
    // Options might be formatted as "A) ... B) ... C) ..." etc.
    options = options.trim();
  }
  
  return {
    id: questionId,
    category: row.category || 'Unknown',
    question: `${row.question}\n\nOptions:\n${options}`,
    correct_answer: row.answer || row.correct_answer || row.answer_index
  };
});

// Create the formatted JSON structure
const formattedData = {
  test_name: "MMLU Questions 31-50",
  created_date: new Date().toISOString(),
  total_questions: 20,
  questions: questions31to50
};

// Save to file
const outputFile = 'mmlu-q31-q50-extracted.json';
fs.writeFileSync(outputFile, JSON.stringify(formattedData, null, 2));

console.log(`\nExtracted ${questions31to50.length} questions`);
console.log(`Saved to: ${outputFile}\n`);

// Show first few questions as preview
console.log('Preview of extracted questions:\n');
questions31to50.slice(0, 3).forEach(q => {
  console.log(`${q.id} (${q.category}):`);
  console.log(`Question: ${q.question.substring(0, 150)}...`);
  console.log(`Answer: ${q.correct_answer}\n`);
});

// Category breakdown
const categories = {};
questions31to50.forEach(q => {
  categories[q.category] = (categories[q.category] || 0) + 1;
});

console.log('Category breakdown:');
Object.entries(categories).forEach(([cat, count]) => {
  console.log(`  ${cat}: ${count} questions`);
});