// Manual extraction of 20 diverse MMLU questions with careful formatting
const XLSX = require('xlsx');
const fs = require('fs');

// Load the Excel file
const workbook = XLSX.readFile('MMLU-PRO Benchmark Questions.xlsx');
const sheetName = workbook.SheetNames[0];
const worksheet = workbook.Sheets[sheetName];

// Convert to JSON
const data = XLSX.utils.sheet_to_json(worksheet);

// Previously used IDs (from your summary)
const usedIds = new Set([71, 72, 866, 867, 1986, 1987, 2804, 2805, 3526, 3527, 
                         70, 73, 868, 869, 1988, 1991, 2807, 2808, 3528, 3530]);

// Target diverse categories
const categories = [
  'physics', 'chemistry', 'biology', 'math',
  'computer science', 'engineering', 'psychology',
  'economics', 'business', 'philosophy', 'history',
  'law', 'sociology', 'medicine', 'geography'
];

const selectedQuestions = [];
const categoryCount = {};

// Helper to parse options more carefully
function parseOptionsCarefully(optionsStr, questionId) {
  console.log(`\nParsing options for question ${questionId}:`);
  console.log(`Raw: ${optionsStr.substring(0, 200)}...`);
  
  // Clean up the string
  let clean = optionsStr.trim();
  
  // Remove outer quotes if present
  if (clean.startsWith('"') && clean.endsWith('"')) {
    clean = clean.slice(1, -1);
  }
  if (clean.startsWith("'") && clean.endsWith("'")) {
    clean = clean.slice(1, -1);
  }
  
  // Replace escaped quotes
  clean = clean.replace(/\\"/g, '"').replace(/\\'/g, "'");
  
  // Try different parsing strategies
  let options = [];
  
  // Strategy 1: Look for A) B) C) pattern
  const letterPattern = /([A-J])\)\s*([^A-J]+?)(?=(?:[A-J]\)|$))/g;
  const matches = [...clean.matchAll(letterPattern)];
  
  if (matches.length >= 2) {
    options = matches.map(m => `${m[1]}) ${m[2].trim()}`);
    console.log(`Strategy 1 worked: ${options.length} options`);
  } else {
    // Strategy 2: Split by newlines
    const lines = clean.split(/\\n|\n/).filter(line => line.trim());
    if (lines.length >= 2) {
      options = lines.map((line, idx) => {
        const letter = String.fromCharCode(65 + idx);
        if (line.match(/^[A-J]\)/)) {
          return line.trim();
        } else {
          return `${letter}) ${line.trim()}`;
        }
      });
      console.log(`Strategy 2 worked: ${options.length} options`);
    }
  }
  
  return options;
}

// Manual selection of good questions
const manualSelection = [
  // Math
  { id: 7695, cat: 'math' },
  { id: 7701, cat: 'math' },
  // Physics  
  { id: 8209, cat: 'physics' },
  { id: 8215, cat: 'physics' },
  // Chemistry
  { id: 3531, cat: 'chemistry' },
  { id: 3535, cat: 'chemistry' },
  // Biology
  { id: 2811, cat: 'biology' },
  { id: 2815, cat: 'biology' },
  // Computer Science
  { id: 2989, cat: 'computer science' },
  { id: 2995, cat: 'computer science' },
  // Psychology
  { id: 1995, cat: 'psychology' },
  { id: 2001, cat: 'psychology' },
  // Economics
  { id: 6831, cat: 'economics' },
  { id: 6835, cat: 'economics' },
  // History
  { id: 4675, cat: 'history' },
  { id: 4681, cat: 'history' },
  // Philosophy
  { id: 4350, cat: 'philosophy' },
  { id: 4355, cat: 'philosophy' },
  // Law
  { id: 875, cat: 'law' },
  { id: 881, cat: 'law' }
];

// Extract the manually selected questions
for (const selection of manualSelection) {
  const row = data.find(r => r.question_id === selection.id);
  if (!row) {
    console.log(`⚠️  Could not find question ${selection.id}`);
    continue;
  }
  
  const question = row.question || row.Question;
  const optionsStr = row.options || row.Options;
  const answer = row.answer || row.Answer;
  
  if (!question || !optionsStr) {
    console.log(`⚠️  Missing data for question ${selection.id}`);
    continue;
  }
  
  const options = parseOptionsCarefully(optionsStr, selection.id);
  
  if (options.length < 2) {
    console.log(`⚠️  Too few options for question ${selection.id}`);
    continue;
  }
  
  // Determine answer letter
  let answerLetter;
  if (typeof answer === 'string' && answer.match(/^[A-J]$/)) {
    answerLetter = answer;
  } else if (typeof answer === 'number') {
    answerLetter = String.fromCharCode(65 + answer);
  } else {
    console.log(`⚠️  Invalid answer format for question ${selection.id}: ${answer}`);
    continue;
  }
  
  selectedQuestions.push({
    question_id: row.question_id,
    category: selection.cat,
    question: question.trim(),
    options: options,
    correct_answer: answerLetter,
    answer_index: answerLetter.charCodeAt(0) - 65
  });
  
  console.log(`✅ Added ${selection.cat} question ${selection.id}`);
}

// Save the cleaned questions
const output = {
  timestamp: new Date().toISOString(),
  total_questions: selectedQuestions.length,
  questions: selectedQuestions
};

fs.writeFileSync('mmlu-20-manual-clean.json', JSON.stringify(output, null, 2));

console.log(`\n✅ Successfully extracted ${selectedQuestions.length} questions`);
console.log('💾 Saved to: mmlu-20-manual-clean.json');

// Display distribution
const dist = {};
selectedQuestions.forEach(q => {
  dist[q.category] = (dist[q.category] || 0) + 1;
});

console.log('\n📊 Category Distribution:');
Object.entries(dist).forEach(([cat, count]) => {
  console.log(`  ${cat}: ${count}`);
});