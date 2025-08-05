const fs = require('fs');

// Let's manually compute what the answer should be
console.log('=== MANUAL VERIFICATION OF QUESTION 21 ===\n');

console.log('Given:');
console.log('f(x,y) = x² + y²');
console.log('g(x,y) = x² - y²');
console.log('Find: f(g(x,y), g(f(x,y)))\n');

console.log('Step 1: Calculate g(x,y)');
console.log('g(x,y) = x² - y²\n');

console.log('Step 2: Calculate f(x,y)');
console.log('f(x,y) = x² + y²\n');

console.log('Step 3: Calculate g(f(x,y))');
console.log('This is where the ambiguity lies.');
console.log('g is a function of two variables: g(a,b) = a² - b²');
console.log('But f(x,y) gives us a single value: x² + y²\n');

console.log('Common interpretations:');
console.log('1. g(f(x,y), f(x,y)) = (x²+y²)² - (x²+y²)² = 0');
console.log('2. g(f(x,y), 0) = (x²+y²)² - 0² = (x²+y²)²');
console.log('3. The problem has a typo and meant f(g(x,y), f(x,y))\n');

console.log('If interpretation 1 (most common):');
console.log('f(g(x,y), g(f(x,y))) = f(x²-y², 0) = (x²-y²)² + 0²');
console.log('= x⁴ - 2x²y² + y⁴');
console.log('This is Option H\n');

console.log('If the problem meant f(g(x,y), f(x,y)):');
console.log('f(x²-y², x²+y²) = (x²-y²)² + (x²+y²)²');
console.log('= (x⁴ - 2x²y² + y⁴) + (x⁴ + 2x²y² + y⁴)');
console.log('= 2x⁴ + 2y⁴');
console.log('This is Option D\n');

console.log('Option B (x⁴ + y⁴) would require:');
console.log('Some calculation that gives x⁴ + y⁴...\n');

// Check all test result files
console.log('\n=== CHECKING ALL TEST RESULTS ===\n');

const testFiles = fs.readdirSync('.').filter(f => 
  f.includes('mmlu') && f.includes('test') && f.endsWith('.json')
);

testFiles.forEach(file => {
  try {
    const data = JSON.parse(fs.readFileSync(file, 'utf-8'));
    
    // Search for question 21 results
    const searchInObject = (obj, path = '') => {
      if (typeof obj === 'object' && obj !== null) {
        for (const [key, value] of Object.entries(obj)) {
          if (key.includes('21') || (typeof value === 'string' && value.includes('f(g(x,y)'))) {
            console.log(`Found in ${file} at ${path}.${key}:`, 
              typeof value === 'string' ? value.substring(0, 100) : value);
          }
          if (typeof value === 'object') {
            searchInObject(value, `${path}.${key}`);
          }
        }
      }
    };
    
    searchInObject(data);
  } catch (e) {
    // Skip
  }
});

// Final check of our formatted questions
console.log('\n=== FINAL VERIFICATION ===\n');
const formatted = JSON.parse(fs.readFileSync('mmlu-30-questions-formatted.json', 'utf-8'));
const q21 = formatted.questions.find(q => q.id === 'mmlu-21');

console.log('Source file: mmlu-30-questions-formatted.json');
console.log('Question ID:', q21.id);
console.log('Category:', q21.category);
console.log('Stated correct answer:', q21.correct_answer);

// Extract just the options
const options = q21.question.split('\n').filter(line => line.match(/^\([A-J]\)/));
console.log('\nAll options:');
options.forEach(opt => console.log(opt));

console.log(`\nThe file states the correct answer is ${q21.correct_answer}, which corresponds to:`);
console.log(options.find(opt => opt.startsWith(`(${q21.correct_answer})`)));