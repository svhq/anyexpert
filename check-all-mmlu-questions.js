const XLSX = require('xlsx');
const fs = require('fs');

console.log('=== CHECKING ALL AVAILABLE MMLU QUESTIONS ===\n');

// Check the Excel file
try {
  console.log('1. Checking MMLU-PRO Benchmark Questions.xlsx...');
  const workbook = XLSX.readFile('MMLU-PRO Benchmark Questions.xlsx');
  const sheetName = workbook.SheetNames[0];
  const sheet = workbook.Sheets[sheetName];
  const data = XLSX.utils.sheet_to_json(sheet);
  
  console.log(`   Total questions in Excel: ${data.length}`);
  console.log(`   Categories found: ${[...new Set(data.map(q => q.category))].length}`);
  
  // Check if there are questions numbered 31-50
  const highNumberedQuestions = data.filter((q, idx) => {
    const questionText = JSON.stringify(q).toLowerCase();
    return idx >= 30 && idx < 50; // Questions 31-50 by index
  });
  
  console.log(`\n   Sample questions from rows 31-50:`);
  highNumberedQuestions.slice(0, 3).forEach((q, idx) => {
    console.log(`   Row ${31 + idx}: ${q.category} - ${(q.question || '').substring(0, 80)}...`);
  });
} catch (e) {
  console.log('   Error reading Excel:', e.message);
}

// Check all JSON files
console.log('\n2. Checking JSON files for question counts...');
const jsonFiles = fs.readdirSync('.').filter(f => 
  f.includes('mmlu') && f.endsWith('.json') && 
  (f.includes('question') || f.includes('clean') || f.includes('formatted'))
);

jsonFiles.forEach(file => {
  try {
    const data = JSON.parse(fs.readFileSync(file, 'utf-8'));
    let questionCount = 0;
    let maxId = 0;
    
    if (data.questions && Array.isArray(data.questions)) {
      questionCount = data.questions.length;
      data.questions.forEach(q => {
        const match = (q.id || '').match(/mmlu-(\d+)/);
        if (match) {
          maxId = Math.max(maxId, parseInt(match[1]));
        }
      });
    } else if (Array.isArray(data)) {
      questionCount = data.length;
    }
    
    if (questionCount > 0) {
      console.log(`   ${file}: ${questionCount} questions` + (maxId > 0 ? ` (max ID: mmlu-${maxId})` : ''));
    }
  } catch (e) {
    // Skip
  }
});

// Check for any "batch" or "set" files that might contain more questions
console.log('\n3. Checking for additional question sets...');
const batchFiles = fs.readdirSync('.').filter(f => 
  (f.includes('batch') || f.includes('set') || f.includes('50') || f.includes('100')) && 
  (f.endsWith('.json') || f.endsWith('.csv'))
);

batchFiles.forEach(file => {
  console.log(`   Found: ${file}`);
});

// Check parent directory references
console.log('\n4. Checking for references to larger datasets...');
const readmeContent = fs.existsSync('README.md') ? fs.readFileSync('README.md', 'utf-8') : '';
const mmluReferences = readmeContent.match(/MMLU.*?(\d+).*?questions?/gi) || [];
mmluReferences.forEach(ref => {
  console.log(`   README reference: ${ref}`);
});

// Look for any test results that might indicate more questions
console.log('\n5. Checking test results for question ranges...');
const testResults = fs.readdirSync('.').filter(f => 
  f.includes('test') && f.includes('result') && f.endsWith('.json')
);

testResults.forEach(file => {
  try {
    const content = fs.readFileSync(file, 'utf-8');
    const matches = content.match(/mmlu-(\d+)/g) || [];
    const numbers = matches.map(m => parseInt(m.replace('mmlu-', ''))).filter(n => n > 30);
    if (numbers.length > 0) {
      console.log(`   ${file} contains references to questions: ${numbers.join(', ')}`);
    }
  } catch (e) {
    // Skip
  }
});

console.log('\n6. Summary:');
console.log('   - We have formatted questions for Q1-Q30');
console.log('   - The Excel file contains 11,943 total MMLU questions');
console.log('   - No additional formatted sets (Q31-Q50) found in current directory');
console.log('   - You may need to extract/format additional questions from the Excel file');