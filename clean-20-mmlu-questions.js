// Clean and format the 20 MMLU questions
const fs = require('fs');

// Load the raw extracted questions
const rawData = JSON.parse(fs.readFileSync('mmlu-20-diverse-questions.json', 'utf8'));

// Function to clean options
function cleanOptions(options) {
  // Join all options into one string first
  const fullText = options.join(' ');
  
  // Remove brackets
  const cleaned = fullText.replace(/[\[\]]/g, '');
  
  // Extract individual options using regex
  const optionPattern = /([A-J])\)\s*([^A-J]+?)(?=(?:[A-J]\)|$))/g;
  const matches = [...cleaned.matchAll(optionPattern)];
  
  if (matches.length > 0) {
    return matches.map(match => `${match[1]}) ${match[2].trim()}`);
  }
  
  // Fallback: try to split by letter patterns
  const parts = cleaned.split(/(?=[A-J]\))/);
  return parts
    .filter(part => part.trim())
    .map(part => part.trim());
}

// Clean all questions
const cleanedQuestions = rawData.questions.map((q, idx) => {
  console.log(`\nProcessing question ${idx + 1}: ${q.question_id}`);
  console.log(`Original options: ${q.options.length}`);
  
  // Clean options
  const cleanedOptions = cleanOptions(q.options);
  console.log(`Cleaned options: ${cleanedOptions.length}`);
  
  // Validate answer index
  const answerIndex = q.correct_answer.charCodeAt(0) - 65;
  if (answerIndex >= cleanedOptions.length) {
    console.warn(`⚠️  Answer index ${answerIndex} out of range for ${cleanedOptions.length} options`);
  }
  
  return {
    question_id: q.question_id,
    category: q.category,
    question: q.question,
    options: cleanedOptions,
    correct_answer: q.correct_answer,
    answer_index: answerIndex
  };
});

// Check for any issues
console.log('\n📊 Validation Summary:');
let issueCount = 0;
cleanedQuestions.forEach(q => {
  if (q.options.length < 2) {
    console.log(`❌ Question ${q.question_id}: Only ${q.options.length} options`);
    issueCount++;
  }
  if (q.answer_index >= q.options.length) {
    console.log(`❌ Question ${q.question_id}: Answer ${q.correct_answer} out of range`);
    issueCount++;
  }
  // Check if options are properly formatted
  q.options.forEach((opt, idx) => {
    if (!opt.match(/^[A-J]\)/)) {
      console.log(`❌ Question ${q.question_id}: Option ${idx} not properly formatted: "${opt}"`);
      issueCount++;
    }
  });
});

if (issueCount === 0) {
  console.log('✅ All questions validated successfully!');
} else {
  console.log(`⚠️  Found ${issueCount} issues that need manual review`);
}

// Save cleaned questions
const output = {
  timestamp: new Date().toISOString(),
  total_questions: cleanedQuestions.length,
  domain_distribution: rawData.domain_distribution,
  questions: cleanedQuestions
};

fs.writeFileSync('mmlu-20-clean-questions.json', JSON.stringify(output, null, 2));
console.log('\n💾 Saved cleaned questions to: mmlu-20-clean-questions.json');

// Display sample cleaned questions
console.log('\n📋 Sample Cleaned Questions:');
cleanedQuestions.slice(0, 3).forEach((q, idx) => {
  console.log(`\n--- Question ${idx + 1} ---`);
  console.log(`ID: ${q.question_id}`);
  console.log(`Category: ${q.category}`);
  console.log(`Question: ${q.question.substring(0, 100)}...`);
  console.log('Options:');
  q.options.forEach(opt => console.log(`  ${opt}`));
  console.log(`Correct Answer: ${q.correct_answer}`);
});