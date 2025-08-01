const fs = require('fs');

// Load the GPQA questions
const data = JSON.parse(fs.readFileSync('gpqa-10-questions.json', 'utf8'));

// Function to clean text for better model understanding
function cleanForModel(text) {
  if (!text) return '';
  
  // Replace LaTeX-style notation with plain text
  return text
    // Replace \uparrow and \downarrow with plain text
    .replace(/\|\\uparrow\\rangle/g, '|up⟩')
    .replace(/\|\\downarrow\\rangle/g, '|down⟩')
    // Replace \sigma{z} notation
    .replace(/\\sigma\{z\}/g, 'σ_z')
    .replace(/\\sigma_\{x\}/g, 'σ_x')
    // Clean up other LaTeX
    .replace(/sqrt\(3\)/g, '√3')
    // Replace special characters that might confuse models
    .replace(/[""]/g, '"')
    .replace(/['']/g, "'")
    // Keep other special characters like π, σ as they are widely understood
    ;
}

// Create cleaned version
const cleanedData = {
  ...data,
  questions: data.questions.map(q => ({
    ...q,
    question: cleanForModel(q.question),
    options: Object.fromEntries(
      Object.entries(q.options).map(([key, value]) => [key, cleanForModel(value)])
    )
  }))
};

// Save cleaned version
fs.writeFileSync('gpqa-10-questions-clean.json', JSON.stringify(cleanedData, null, 2));

// Show sample of how questions will look
console.log('CLEANED GPQA QUESTIONS - READY FOR MODEL TESTING\n');
console.log('=' + '='.repeat(70) + '\n');

// Show the problematic questions after cleaning
const problematicIds = ['gpqa-5', 'gpqa-7', 'gpqa-9'];
cleanedData.questions
  .filter(q => problematicIds.includes(q.id))
  .forEach((q, index) => {
    console.log(`CLEANED QUESTION (${q.id}):`);
    console.log('-'.repeat(70));
    
    const formattedQuestion = `Question: ${q.question}

Options:

A) ${q.options.A}
B) ${q.options.B}
C) ${q.options.C}
D) ${q.options.D}`;

    console.log(formattedQuestion);
    console.log('\nCorrect Answer:', q.correct_answer);
    console.log('\n' + '='.repeat(70) + '\n');
  });

console.log('✓ Created gpqa-10-questions-clean.json with cleaned formatting');
console.log('✓ All questions are now in simple, clear format for model testing');