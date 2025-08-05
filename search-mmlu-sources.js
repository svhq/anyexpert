const fs = require('fs');
const path = require('path');

// Check all our MMLU JSON files
const mmluFiles = [
  'mmlu-30-questions-formatted.json',
  'mmlu-30-questions-clean.json',
  'mmlu-20-final-clean.json',
  'mmlu-batch1-questions.json'
];

console.log('Searching for Question 21 across all MMLU files...\n');

mmluFiles.forEach(filename => {
  if (fs.existsSync(filename)) {
    console.log(`\nChecking ${filename}:`);
    try {
      const data = JSON.parse(fs.readFileSync(filename, 'utf-8'));
      
      // Look for question 21
      let q21 = null;
      if (data.questions) {
        q21 = data.questions.find(q => q.id === 'mmlu-21' || q.question_id === 'mmlu-21');
      } else if (Array.isArray(data)) {
        q21 = data.find(q => q.id === 'mmlu-21' || q.question_id === 'mmlu-21');
      }
      
      if (q21) {
        console.log('Found Question 21!');
        console.log('Category:', q21.category);
        console.log('Question:', q21.question.substring(0, 100) + '...');
        console.log('Correct Answer:', q21.correct_answer);
        
        // Check if it mentions the specific formulas
        if (q21.question.includes('f(g(x,y), g(f(x,y)))')) {
          console.log('✓ Contains f(g(x,y), g(f(x,y)))');
        }
        if (q21.question.includes('x² + y²') || q21.question.includes('x^2 + y^2')) {
          console.log('✓ Contains f(x,y) = x² + y²');
        }
        if (q21.question.includes('x² - y²') || q21.question.includes('x^2 - y^2')) {
          console.log('✓ Contains g(x,y) = x² - y²');
        }
        
        // Show options B and D specifically
        const lines = q21.question.split('\n');
        lines.forEach(line => {
          if (line.includes('(B)')) console.log('Option B:', line.trim());
          if (line.includes('(D)')) console.log('Option D:', line.trim());
        });
      } else {
        console.log('Question 21 not found in this file');
      }
    } catch (e) {
      console.log('Error reading file:', e.message);
    }
  }
});

// Also check if there's any documentation about where these questions came from
console.log('\n\nLooking for any README or documentation files...');
const readmeFiles = fs.readdirSync('.').filter(f => 
  f.toLowerCase().includes('readme') || 
  f.toLowerCase().includes('mmlu') && f.endsWith('.md') ||
  f.toLowerCase().includes('source')
);

readmeFiles.forEach(file => {
  console.log('\nFound:', file);
  const content = fs.readFileSync(file, 'utf-8');
  if (content.includes('MMLU') || content.includes('source') || content.includes('dataset')) {
    console.log('Relevant excerpt:', content.substring(0, 300) + '...');
  }
});

// Check superagent directories
console.log('\n\nChecking parent superagent directories for MMLU data...');
const parentDirs = fs.readdirSync('..').filter(d => d.includes('mmlu') || d.includes('superagent'));
parentDirs.forEach(dir => {
  console.log('Found directory:', dir);
});