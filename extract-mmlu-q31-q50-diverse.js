const XLSX = require('xlsx');
const fs = require('fs');

// Read the Excel file
const workbook = XLSX.readFile('MMLU-PRO Benchmark Questions.xlsx');
const sheetName = workbook.SheetNames[0];
const sheet = workbook.Sheets[sheetName];
const data = XLSX.utils.sheet_to_json(sheet);

console.log('=== EXTRACTING DIVERSE MMLU QUESTIONS (Q31-Q50) ===\n');

// Get all unique categories
const allCategories = [...new Set(data.map(q => q.category))];
console.log(`Available categories: ${allCategories.join(', ')}\n`);

// Try to get a diverse set - different categories
const categoriesNeeded = [
  'chemistry', 'law', 'health', 'business', 'physics', 
  'math', 'biology', 'psychology', 'computer science', 'economics',
  'engineering', 'history', 'philosophy', 'other'
];

const diverseQuestions = [];
let questionNumber = 31;

// First, try to get 2 questions from each of our preferred categories
for (const category of categoriesNeeded) {
  const categoryQuestions = data.filter(q => 
    q.category && q.category.toLowerCase().includes(category.toLowerCase())
  );
  
  if (categoryQuestions.length > 0 && diverseQuestions.length < 20) {
    // Take up to 2 questions from this category
    const questionsToTake = Math.min(2, 20 - diverseQuestions.length);
    for (let i = 0; i < questionsToTake && i < categoryQuestions.length; i++) {
      const q = categoryQuestions[i];
      diverseQuestions.push({
        id: `mmlu-${questionNumber++}`,
        category: q.category,
        question: `${q.question}\n\nOptions:\n${q.options}`,
        correct_answer: q.answer || q.correct_answer || q.answer_index
      });
    }
  }
}

// If we still need more questions, take from any category
if (diverseQuestions.length < 20) {
  const startIndex = 100; // Start from a different section
  for (let i = startIndex; i < data.length && diverseQuestions.length < 20; i++) {
    const q = data[i];
    if (!diverseQuestions.find(dq => dq.question === q.question)) {
      diverseQuestions.push({
        id: `mmlu-${questionNumber++}`,
        category: q.category,
        question: `${q.question}\n\nOptions:\n${q.options}`,
        correct_answer: q.answer || q.correct_answer || q.answer_index
      });
    }
  }
}

// Create the formatted JSON structure
const formattedData = {
  test_name: "MMLU Questions 31-50 (Diverse Categories)",
  created_date: new Date().toISOString(),
  total_questions: diverseQuestions.length,
  questions: diverseQuestions.slice(0, 20) // Ensure exactly 20
};

// Save to file
const outputFile = 'mmlu-q31-q50-diverse.json';
fs.writeFileSync(outputFile, JSON.stringify(formattedData, null, 2));

console.log(`Extracted ${formattedData.questions.length} questions`);
console.log(`Saved to: ${outputFile}\n`);

// Show category breakdown
const categoryBreakdown = {};
formattedData.questions.forEach(q => {
  categoryBreakdown[q.category] = (categoryBreakdown[q.category] || 0) + 1;
});

console.log('Category breakdown:');
Object.entries(categoryBreakdown).sort().forEach(([cat, count]) => {
  console.log(`  ${cat}: ${count} questions`);
});

// Preview first few questions
console.log('\nPreview of questions:');
formattedData.questions.slice(0, 5).forEach(q => {
  console.log(`\n${q.id} (${q.category}):`);
  console.log(`Question: ${q.question.substring(0, 100)}...`);
  console.log(`Answer: ${q.correct_answer}`);
});