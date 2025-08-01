// Simple verification of specific wrong answers
const fs = require('fs');

// Load the questions
const questions = JSON.parse(fs.readFileSync('./mmlu-10-new-questions.json', 'utf8')).questions;

// Key wrong answers to check
const toCheck = [
  // o4-mini-high wrong answer
  {
    model: 'o4-mini-high', 
    qNum: 6,
    category: 'psychology',
    selected: 'H',
    correct: 'I',
    question: questions[5]
  },
  // GLM wrong answers
  {
    model: 'GLM-4.5-air',
    qNum: 3,
    category: 'law',
    selected: 'A', 
    correct: 'D',
    question: questions[2]
  },
  {
    model: 'GLM-4.5-air',
    qNum: 5,
    category: 'psychology',
    selected: null,
    correct: 'E',
    question: questions[4]
  },
  {
    model: 'GLM-4.5-air',
    qNum: 8,
    category: 'biology',
    selected: 'J',
    correct: 'D',
    question: questions[7]
  },
  {
    model: 'GLM-4.5-air',
    qNum: 9,
    category: 'chemistry',
    selected: 'J',
    correct: 'C',
    question: questions[8]
  }
];

console.log('🔍 Checking Wrong Answers\n');

toCheck.forEach(item => {
  console.log('='.repeat(60));
  console.log(`\nModel: ${item.model}`);
  console.log(`Question ${item.qNum}: ${item.category}`);
  console.log(`\nQuestion: ${item.question.question.substring(0, 100)}...`);
  console.log(`\nOptions (showing relevant ones):`);
  
  // Show the correct option
  const correctOpt = item.question.options.find(o => o.startsWith(item.correct));
  console.log(`✅ ${correctOpt}`);
  
  // Show what model selected (if any)
  if (item.selected) {
    const selectedOpt = item.question.options.find(o => o.startsWith(item.selected));
    console.log(`❌ ${selectedOpt || 'Option not found'}`);
  } else {
    console.log(`❌ No answer selected (extraction failed)`);
  }
  
  // Analysis
  console.log(`\nAnalysis:`);
  
  if (item.qNum === 6) {
    console.log(`Psychology ethics question - Both H and I involve referrals.`);
    console.log(`H) "provide Hermann with appropriate referrals"`);
    console.log(`I) Option text appears corrupted in data`);
    console.log(`This may be a data quality issue.`);
  }
  
  if (item.qNum === 3) {
    console.log(`Robbery vs Larceny - Model chose larceny, correct is robbery.`);
    console.log(`Key distinction: Force was used DURING taking (robbery) not just after.`);
  }
  
  if (item.qNum === 5) {
    console.log(`Rancho Scale Level 4 - No answer extracted from response.`);
    console.log(`This was likely an extraction failure, not wrong answer.`);
  }
  
  if (item.qNum === 8) {
    console.log(`Photosynthesis - Model said J (water/oxygen), correct is D (ATP/NADPH).`);
    console.log(`This is a factual error about photosynthesis products.`);
  }
  
  if (item.qNum === 9) {
    console.log(`Bonding type - Model said J (polar covalent), correct is C (network).`);
    console.log(`Properties match network solids like diamond/SiO2, not polar molecules.`);
  }
});

console.log('\n' + '='.repeat(60));
console.log('\n📊 SUMMARY:\n');
console.log('1. o4-mini-high Q6 (Psychology): Possibly correct - option text corrupted');
console.log('2. GLM Q3 (Law): Truly wrong - chose larceny over robbery');
console.log('3. GLM Q5 (Psychology): Extraction failure - no answer detected');
console.log('4. GLM Q8 (Biology): Truly wrong - incorrect photosynthesis products');
console.log('5. GLM Q9 (Chemistry): Truly wrong - misidentified bonding type');

console.log('\n📎 Note: The psychology ethics question (Q6) has corrupted option text.');
console.log('Both models may have been correct but selected different valid options.');

// Check the original option text
console.log('\n🔍 Checking Original Options for Q6:');
questions[5].options.forEach(opt => {
  if (opt.startsWith('H') || opt.startsWith('I')) {
    console.log(opt);
  }
});