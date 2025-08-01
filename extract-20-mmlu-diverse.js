// Extract 20 diverse MMLU questions from Excel
const XLSX = require('xlsx');
const fs = require('fs');

// Load the Excel file
const workbook = XLSX.readFile('MMLU-PRO Benchmark Questions.xlsx');
const sheetName = workbook.SheetNames[0];
const worksheet = workbook.Sheets[sheetName];

// Convert to JSON
const data = XLSX.utils.sheet_to_json(worksheet);

// Track already used question IDs
const usedQuestionIds = new Set();

// Load previously used questions to avoid duplicates
const previousFiles = [
  'mmlu-batch1-questions.json',
  'mmlu-10-new-questions.json',
  'mmlu_batch2_o4_mini_high_test_2025-07-31T04-04-14.json'
];

previousFiles.forEach(file => {
  if (fs.existsSync(file)) {
    try {
      const content = JSON.parse(fs.readFileSync(file, 'utf8'));
      if (content.questions) {
        content.questions.forEach(q => usedQuestionIds.add(q.question_id));
      }
    } catch (e) {
      console.log(`Could not load ${file}: ${e.message}`);
    }
  }
});

console.log(`Already used question IDs: ${usedQuestionIds.size}`);

// Target domains for diversity
const targetDomains = [
  'physics', 'chemistry', 'biology', 'math',
  'computer science', 'engineering', 'medicine', 'psychology',
  'economics', 'business', 'philosophy', 'history',
  'law', 'sociology', 'political science', 'geography',
  'literature', 'art', 'music', 'education'
];

// Function to clean and parse options
function parseOptions(optionsStr) {
  if (!optionsStr) return [];
  
  // Remove quotes and split by newlines or commas
  const cleanStr = optionsStr.replace(/['"]/g, '').trim();
  
  // Try to split by common patterns
  let options = [];
  
  // Pattern 1: Letter followed by parenthesis
  const letterPattern = /([A-J])\)\s*([^A-J]+?)(?=(?:[A-J]\)|$))/g;
  const matches = [...cleanStr.matchAll(letterPattern)];
  
  if (matches.length > 0) {
    options = matches.map(match => `${match[1]}) ${match[2].trim()}`);
  } else {
    // Pattern 2: Try splitting by newlines
    options = cleanStr.split(/\n/).filter(opt => opt.trim());
    
    // If no newlines, try other delimiters
    if (options.length < 2) {
      options = cleanStr.split(/(?=[A-J]\))/).filter(opt => opt.trim());
    }
  }
  
  // Ensure each option starts with a letter and parenthesis
  options = options.map((opt, idx) => {
    const letterMatch = opt.match(/^([A-J])\)/);
    if (letterMatch) {
      return opt.trim();
    } else {
      // Add letter if missing
      const letter = String.fromCharCode(65 + idx); // A, B, C...
      return `${letter}) ${opt.trim()}`;
    }
  });
  
  return options.filter(opt => opt && opt.length > 3); // Filter out empty options
}

// Extract questions
const extractedQuestions = [];
const domainCounts = {};

for (const row of data) {
  // Skip if already used
  if (usedQuestionIds.has(row.question_id)) continue;
  
  // Get domain/category
  const category = (row.category || '').toLowerCase().trim();
  if (!category) continue;
  
  // Check if we need more from this domain
  const currentCount = domainCounts[category] || 0;
  if (currentCount >= 2) continue; // Max 2 per domain for diversity
  
  // Extract question data
  const questionText = row.question || row.Question;
  const optionsStr = row.options || row.Options;
  const correctAnswer = row.answer || row.Answer || row.correct_answer;
  
  if (!questionText || !optionsStr || correctAnswer === undefined) continue;
  
  // Parse options
  const options = parseOptions(optionsStr);
  if (options.length < 2) continue; // Skip if too few options
  
  // Determine correct answer letter
  let answerLetter;
  if (typeof correctAnswer === 'string' && correctAnswer.match(/^[A-J]$/)) {
    answerLetter = correctAnswer;
  } else if (typeof correctAnswer === 'number') {
    answerLetter = String.fromCharCode(65 + correctAnswer); // Convert index to letter
  } else {
    continue; // Skip if we can't determine answer
  }
  
  // Create clean question object
  const questionObj = {
    question_id: row.question_id,
    category: category,
    question: questionText.trim(),
    options: options,
    correct_answer: answerLetter,
    answer_index: answerLetter.charCodeAt(0) - 65,
    cot_content: row.cot_content || ""
  };
  
  // Validate question
  if (questionObj.options.length >= 2 && 
      questionObj.answer_index < questionObj.options.length) {
    extractedQuestions.push(questionObj);
    domainCounts[category] = currentCount + 1;
    
    console.log(`✓ Added question ${questionObj.question_id} from ${category} (${extractedQuestions.length}/20)`);
    
    if (extractedQuestions.length >= 20) break;
  }
}

// If we don't have enough, add more from any domain
if (extractedQuestions.length < 20) {
  for (const row of data) {
    if (usedQuestionIds.has(row.question_id)) continue;
    if (extractedQuestions.some(q => q.question_id === row.question_id)) continue;
    
    const category = (row.category || '').toLowerCase().trim();
    if (!category) continue;
    
    const questionText = row.question || row.Question;
    const optionsStr = row.options || row.Options;
    const correctAnswer = row.answer || row.Answer || row.correct_answer;
    
    if (!questionText || !optionsStr || correctAnswer === undefined) continue;
    
    const options = parseOptions(optionsStr);
    if (options.length < 2) continue;
    
    let answerLetter;
    if (typeof correctAnswer === 'string' && correctAnswer.match(/^[A-J]$/)) {
      answerLetter = correctAnswer;
    } else if (typeof correctAnswer === 'number') {
      answerLetter = String.fromCharCode(65 + correctAnswer);
    } else {
      continue;
    }
    
    const questionObj = {
      question_id: row.question_id,
      category: category,
      question: questionText.trim(),
      options: options,
      correct_answer: answerLetter,
      answer_index: answerLetter.charCodeAt(0) - 65,
      cot_content: row.cot_content || ""
    };
    
    if (questionObj.options.length >= 2 && 
        questionObj.answer_index < questionObj.options.length) {
      extractedQuestions.push(questionObj);
      console.log(`✓ Added question ${questionObj.question_id} from ${category} (${extractedQuestions.length}/20)`);
      
      if (extractedQuestions.length >= 20) break;
    }
  }
}

// Display domain distribution
console.log('\n📊 Domain Distribution:');
const finalDomainCounts = {};
extractedQuestions.forEach(q => {
  finalDomainCounts[q.category] = (finalDomainCounts[q.category] || 0) + 1;
});

Object.entries(finalDomainCounts)
  .sort((a, b) => b[1] - a[1])
  .forEach(([domain, count]) => {
    console.log(`  ${domain}: ${count} questions`);
  });

// Save to file
const output = {
  timestamp: new Date().toISOString(),
  total_questions: extractedQuestions.length,
  domain_distribution: finalDomainCounts,
  questions: extractedQuestions
};

fs.writeFileSync('mmlu-20-diverse-questions.json', JSON.stringify(output, null, 2));

console.log(`\n✅ Extracted ${extractedQuestions.length} questions`);
console.log('💾 Saved to: mmlu-20-diverse-questions.json');

// Display sample question for verification
if (extractedQuestions.length > 0) {
  console.log('\n📋 Sample Question:');
  const sample = extractedQuestions[0];
  console.log(`ID: ${sample.question_id}`);
  console.log(`Category: ${sample.category}`);
  console.log(`Question: ${sample.question}`);
  console.log('Options:');
  sample.options.forEach(opt => console.log(`  ${opt}`));
  console.log(`Correct Answer: ${sample.correct_answer}`);
}