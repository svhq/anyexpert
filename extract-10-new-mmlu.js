const xlsx = require('xlsx');
const fs = require('fs');
const path = require('path');

// Load previously tested questions to avoid duplicates
function getTestedQuestions() {
  const tested = new Set();
  
  // Original 10 questions
  const original = JSON.parse(fs.readFileSync('./mmlu_test_results_2025-07-30T16-45-40-316Z.json', 'utf8'));
  original.results.forEach(r => {
    tested.add(r.question.toLowerCase().substring(0, 100)); // Use first 100 chars as key
  });
  
  // Batch test 10 questions
  const batch = JSON.parse(fs.readFileSync('./mmlu_batch_test_2025-07-31T02-07-13-672Z.json', 'utf8'));
  batch.batch1.results.forEach(r => tested.add(r.question.toLowerCase().substring(0, 100)));
  batch.batch2.results.forEach(r => tested.add(r.question.toLowerCase().substring(0, 100)));
  
  console.log(`📋 Found ${tested.size} previously tested questions to exclude`);
  return tested;
}

async function extractNewQuestions() {
  try {
    const filePath = path.join(__dirname, 'MMLU-PRO Benchmark Questions.xlsx');
    const workbook = xlsx.readFile(filePath);
    const sheetName = workbook.SheetNames[0];
    const worksheet = workbook.Sheets[sheetName];
    const data = xlsx.utils.sheet_to_json(worksheet);
    
    console.log(`📊 Total questions in Excel: ${data.length}`);
    
    const testedQuestions = getTestedQuestions();
    const selectedQuestions = [];
    const categoryCounts = {};
    
    // Target categories for diversity
    const targetCategories = [
      'math', 'physics', 'chemistry', 'biology', 'medicine', 
      'computer science', 'engineering', 'business', 'economics', 
      'law', 'history', 'geography', 'politics', 'psychology', 
      'sociology', 'philosophy', 'other'
    ];
    
    // Process each row
    for (const row of data) {
      if (selectedQuestions.length >= 10) break;
      
      const category = (row.category || row.Category || row.subject || 'other').toLowerCase();
      const question = row.question || row.Question || '';
      const questionKey = question.toLowerCase().substring(0, 100);
      
      // Skip if already tested
      if (testedQuestions.has(questionKey)) continue;
      
      // Skip if we already have 2 from this category
      if (categoryCounts[category] >= 2) continue;
      
      // Parse options
      let options = [];
      let optionLabels = [];
      
      if (row.options) {
        try {
          const optionsStr = row.options.toString();
          
          // Extract individual options - they seem to be in format: 'option1' 'option2' etc
          const matches = optionsStr.match(/'([^']+)'/g);
          if (matches) {
            options = matches.map(m => m.replace(/'/g, '').trim());
            // Generate labels A, B, C, etc.
            optionLabels = options.map((_, i) => String.fromCharCode(65 + i));
          }
        } catch (e) {
          console.log(`⚠️ Failed to parse options for question ${row.question_id}`);
          continue;
        }
      }
      
      // Validate question data
      if (!question || options.length < 2 || !row.answer) continue;
      
      // Add to selected questions
      selectedQuestions.push({
        question_id: row.question_id || `Q${selectedQuestions.length + 1}`,
        category: category,
        question: question,
        options: options.map((opt, i) => `${optionLabels[i]}) ${opt}`),
        correct_answer: row.answer || row.Answer,
        answer_index: row.answer_index,
        cot_content: row.cot_content || row.explanation || ''
      });
      
      categoryCounts[category] = (categoryCounts[category] || 0) + 1;
      console.log(`✅ Added question from ${category} (${selectedQuestions.length}/10)`);
    }
    
    // Display selected questions
    console.log('\n🎯 Selected 10 New MMLU Questions:\n');
    selectedQuestions.forEach((q, i) => {
      console.log(`${i + 1}. [${q.category.toUpperCase()}] ${q.question.substring(0, 80)}...`);
      console.log(`   Answer: ${q.correct_answer}`);
    });
    
    // Save to JSON
    const output = {
      timestamp: new Date().toISOString(),
      total_questions: selectedQuestions.length,
      questions: selectedQuestions
    };
    
    fs.writeFileSync('mmlu-10-new-questions.json', JSON.stringify(output, null, 2));
    console.log('\n📁 Saved to mmlu-10-new-questions.json');
    
    // Show category distribution
    console.log('\n📊 Category Distribution:');
    Object.entries(categoryCounts).forEach(([cat, count]) => {
      console.log(`   ${cat}: ${count} questions`);
    });
    
    return selectedQuestions;
    
  } catch (error) {
    console.error('❌ Error:', error.message);
    return null;
  }
}

if (require.main === module) {
  extractNewQuestions();
}

module.exports = { extractNewQuestions };