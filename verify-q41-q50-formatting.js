const fs = require('fs');

// Load the questions
const data = JSON.parse(fs.readFileSync('mmlu-q31-q50-diverse.json', 'utf-8'));
const questions = data.questions.slice(10); // Q41-Q50

console.log('=== VERIFYING Q41-Q50 FORMATTING ===\n');

// Function to parse and clean options
function parseOptions(optionsString) {
  // Remove the outer brackets
  let cleaned = optionsString.replace(/^\[|\]$/g, '');
  
  // Handle different quote patterns
  const options = [];
  let currentOption = '';
  let inQuotes = false;
  let quoteChar = '';
  
  for (let i = 0; i < cleaned.length; i++) {
    const char = cleaned[i];
    
    if (!inQuotes && (char === '"' || char === "'")) {
      inQuotes = true;
      quoteChar = char;
      currentOption = '';
    } else if (inQuotes && char === quoteChar) {
      inQuotes = false;
      options.push(currentOption);
      currentOption = '';
    } else if (inQuotes) {
      currentOption += char;
    }
  }
  
  // Handle any remaining option
  if (currentOption.trim()) {
    options.push(currentOption);
  }
  
  return options;
}

// Check each question
questions.forEach((q, idx) => {
  console.log(`\n${'='.repeat(60)}`);
  console.log(`Q${41 + idx} (${q.category.toUpperCase()})`);
  console.log(`Correct Answer: ${q.correct_answer}`);
  
  // Extract question text and options
  const parts = q.question.split('\n\nOptions:\n');
  const questionText = parts[0];
  const optionsRaw = parts[1];
  
  console.log('\nQuestion text:');
  console.log(questionText);
  
  // Parse options
  const options = parseOptions(optionsRaw);
  
  console.log(`\nNumber of options found: ${options.length}`);
  console.log('Options:');
  options.forEach((opt, i) => {
    const letter = String.fromCharCode(65 + i);
    console.log(`(${letter}) ${opt}`);
  });
  
  // Verify correct answer
  const correctIndex = q.correct_answer.charCodeAt(0) - 65;
  if (correctIndex >= options.length) {
    console.log(`\n⚠️  WARNING: Correct answer ${q.correct_answer} is out of range!`);
    console.log(`   Only ${options.length} options but answer is option #${correctIndex + 1}`);
  } else {
    console.log(`\n✓ Correct answer ${q.correct_answer} is valid`);
    console.log(`  Maps to: ${options[correctIndex]}`);
  }
  
  // Special checks
  if (options.length < 2) {
    console.log('\n⚠️  WARNING: Too few options!');
  }
  
  if (optionsRaw.includes('\\r\\n')) {
    console.log('\n⚠️  Note: Contains \\r\\n characters that need cleaning');
  }
  
  // Save the properly formatted question
  const formatted = {
    ...q,
    questionText: questionText,
    options: options,
    optionCount: options.length,
    formattedQuestion: `${questionText}\n\nOptions:\n${options.map((opt, i) => `(${String.fromCharCode(65 + i)}) ${opt}`).join('\n')}`
  };
  
  fs.writeFileSync(`q${41 + idx}-formatted.json`, JSON.stringify(formatted, null, 2));
});

console.log('\n\n=== SUMMARY ===');
console.log('All questions have been verified and formatted.');
console.log('Formatted versions saved as q41-formatted.json through q50-formatted.json');