const xlsx = require('xlsx');
const path = require('path');

async function readMMLUQuestions() {
  try {
    // Read the Excel file
    const filePath = path.join(__dirname, 'MMLU-PRO Benchmark Questions.xlsx');
    const workbook = xlsx.readFile(filePath);
    
    // Get the first sheet name
    const sheetName = workbook.SheetNames[0];
    const worksheet = workbook.Sheets[sheetName];
    
    // Convert sheet to JSON
    const data = xlsx.utils.sheet_to_json(worksheet);
    
    console.log(`📊 Found ${data.length} MMLU Pro questions`);
    console.log('\n🔍 Sample data structure:');
    console.log(JSON.stringify(data[0], null, 2));
    
    // Get unique categories
    const categories = [...new Set(data.map(row => row.category || row.Category || row.subject || row.Subject))];
    console.log(`\n📚 Available categories: ${categories.join(', ')}`);
    
    // Select 3 random questions from different categories
    const selectedQuestions = [];
    const usedCategories = new Set();
    
    // Shuffle the data
    const shuffled = data.sort(() => Math.random() - 0.5);
    
    for (const row of shuffled) {
      const category = row.category || row.Category || row.subject || row.Subject;
      
      if (!usedCategories.has(category) && selectedQuestions.length < 3) {
        // Parse options from string format
        let options = [];
        if (row.options) {
          try {
            // The options appear to be in a string format with single quotes and \r\n separators
            const optionsStr = row.options.toString();
            
            // Try different parsing approaches
            if (optionsStr.includes("'") && optionsStr.includes('\r\n')) {
              // Split by \r\n and clean up quotes
              options = optionsStr.split('\r\n')
                .map(opt => opt.trim().replace(/^'|'$/g, ''))
                .filter(opt => opt.length > 0);
            } else if (optionsStr.includes('[') && optionsStr.includes(']')) {
              // Try JSON-like parsing
              const matches = optionsStr.match(/\[(.*)\]/);
              if (matches) {
                const cleaned = matches[1].replace(/'/g, '"');
                const parsed = JSON.parse(`[${cleaned}]`);
                options = parsed;
              }
            }
          } catch (e) {
            console.log('Failed to parse options for question:', row.question_id, ':', row.options);
          }
        }
        
        // Only include questions with valid options and complete data
        if (options.length > 0 && row.question && row.answer) {
          selectedQuestions.push({
            question_id: row.question_id,
            category: category,
            question: row.question || row.Question || row.prompt || row.Prompt,
            options: options,
            correct_answer: row.answer || row.Answer || row.correct || row.Correct,
            answer_index: row.answer_index,
            explanation: row.explanation || row.Explanation || ''
          });
          usedCategories.add(category);
        }
      }
    }
    
    console.log('\n🎯 Selected 3 Random Questions from Different Categories:\n');
    
    selectedQuestions.forEach((q, i) => {
      console.log(`\n=== QUESTION ${i + 1} ===`);
      console.log(`📂 Category: ${q.category}`);
      console.log(`❓ Question: ${q.question}`);
      console.log('\n📝 Options:');
      q.options.forEach((option, idx) => {
        const letter = String.fromCharCode(65 + idx); // A, B, C, D
        console.log(`   ${letter}) ${option}`);
      });
      console.log(`\n✅ Correct Answer: ${q.correct_answer}`);
      if (q.explanation) {
        console.log(`💡 Explanation: ${q.explanation}`);
      }
      console.log('\n' + '─'.repeat(80));
    });
    
    return selectedQuestions;
    
  } catch (error) {
    console.error('❌ Error reading MMLU Excel file:', error.message);
    
    if (error.message.includes('Cannot resolve module')) {
      console.log('\n💡 Installing xlsx package...');
      const { execSync } = require('child_process');
      try {
        execSync('npm install xlsx', { stdio: 'inherit' });
        console.log('✅ xlsx package installed. Please run the script again.');
      } catch (installError) {
        console.error('❌ Failed to install xlsx package:', installError.message);
      }
    }
    
    return null;
  }
}

// Run the function if called directly
if (require.main === module) {
  readMMLUQuestions();
}

module.exports = { readMMLUQuestions };