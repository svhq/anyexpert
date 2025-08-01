const fs = require('fs');
const csv = require('csv-parser');

// Function to clean text by removing extra whitespace and quotes
function cleanText(text) {
  if (!text) return '';
  return text.trim()
    .replace(/\s+/g, ' ')
    .replace(/^["']|["']$/g, '')
    .replace(/\\"/g, '"');
}

// Function to check if a question is clean (no obvious errors)
function isCleanQuestion(row) {
  // Check if essential fields exist
  if (!row['Question'] || !row['Correct Answer'] || 
      !row['Incorrect Answer 1'] || !row['Incorrect Answer 2'] || 
      !row['Incorrect Answer 3']) {
    return false;
  }
  
  // Check if question is not too short or malformed
  if (row['Question'].length < 50) return false;
  
  // Check if answers are distinct and not empty
  const answers = [
    row['Correct Answer'],
    row['Incorrect Answer 1'],
    row['Incorrect Answer 2'],
    row['Incorrect Answer 3']
  ];
  
  const cleanAnswers = answers.map(a => cleanText(a));
  if (cleanAnswers.some(a => a.length < 2)) return false;
  if (new Set(cleanAnswers).size !== 4) return false; // Check all answers are unique
  
  return true;
}

// Collect questions from all GPQA files
async function extractGPQAQuestions() {
  const questions = [];
  const files = [
    'GPQA dataset/gpqa_main.csv',
    'GPQA dataset/gpqa_diamond.csv',
    'GPQA dataset/gpqa_experts.csv',
    'GPQA dataset/gpqa_extended.csv'
  ];
  
  for (const file of files) {
    await new Promise((resolve, reject) => {
      const fileQuestions = [];
      
      fs.createReadStream(file)
        .pipe(csv())
        .on('data', (row) => {
          if (isCleanQuestion(row)) {
            fileQuestions.push({
              question: cleanText(row['Question']),
              correct_answer: cleanText(row['Correct Answer']),
              incorrect_answers: [
                cleanText(row['Incorrect Answer 1']),
                cleanText(row['Incorrect Answer 2']),
                cleanText(row['Incorrect Answer 3'])
              ],
              explanation: cleanText(row['Explanation'] || ''),
              subdomain: row['Subdomain'] || '',
              difficulty: row["Writer's Difficulty Estimate"] || '',
              domain: row['High-level domain'] || '',
              source_file: file
            });
          }
        })
        .on('end', () => {
          console.log(`Found ${fileQuestions.length} clean questions in ${file}`);
          questions.push(...fileQuestions);
          resolve();
        })
        .on('error', reject);
    });
  }
  
  return questions;
}

// Select diverse questions
function selectDiverseQuestions(questions, count = 10) {
  // Group by domain/subdomain
  const byDomain = {};
  questions.forEach(q => {
    const key = q.domain || 'Unknown';
    if (!byDomain[key]) byDomain[key] = [];
    byDomain[key].push(q);
  });
  
  // Select questions ensuring diversity
  const selected = [];
  const domains = Object.keys(byDomain);
  
  // Try to get at least one from each domain
  domains.forEach(domain => {
    if (selected.length < count && byDomain[domain].length > 0) {
      selected.push(byDomain[domain][0]);
    }
  });
  
  // Fill remaining slots
  while (selected.length < count) {
    for (const domain of domains) {
      if (selected.length >= count) break;
      const domainQuestions = byDomain[domain];
      const unselected = domainQuestions.filter(q => !selected.includes(q));
      if (unselected.length > 0) {
        selected.push(unselected[0]);
      }
    }
  }
  
  return selected.slice(0, count);
}

// Format questions for testing
function formatForTesting(questions) {
  return questions.map((q, index) => {
    // Shuffle answers and assign letters
    const allAnswers = [q.correct_answer, ...q.incorrect_answers];
    const shuffled = allAnswers.sort(() => Math.random() - 0.5);
    
    const options = {};
    const letters = ['A', 'B', 'C', 'D'];
    shuffled.forEach((answer, i) => {
      options[letters[i]] = answer;
    });
    
    // Find correct letter
    const correctLetter = letters[shuffled.indexOf(q.correct_answer)];
    
    return {
      id: `gpqa-${index + 1}`,
      category: q.subdomain || q.domain || 'Science',
      question: q.question,
      options,
      correct_answer: correctLetter,
      metadata: {
        domain: q.domain,
        subdomain: q.subdomain,
        difficulty: q.difficulty,
        source_file: q.source_file,
        explanation: q.explanation
      }
    };
  });
}

// Main execution
async function main() {
  console.log('Extracting GPQA questions...\n');
  
  try {
    // Install csv-parser if needed
    const { execSync } = require('child_process');
    try {
      require.resolve('csv-parser');
    } catch (e) {
      console.log('Installing csv-parser...');
      execSync('npm install csv-parser', { stdio: 'inherit' });
    }
    
    const allQuestions = await extractGPQAQuestions();
    console.log(`\nTotal clean questions found: ${allQuestions.length}`);
    
    const selected = selectDiverseQuestions(allQuestions, 10);
    console.log(`\nSelected ${selected.length} diverse questions`);
    
    const formatted = formatForTesting(selected);
    
    // Save to file
    const output = {
      test_name: "GPQA 10 Question Sample",
      created_date: new Date().toISOString(),
      total_questions: formatted.length,
      questions: formatted
    };
    
    fs.writeFileSync('gpqa-10-questions.json', JSON.stringify(output, null, 2));
    
    console.log('\nQuestions saved to gpqa-10-questions.json');
    console.log('\nDomain distribution:');
    const domainCounts = {};
    formatted.forEach(q => {
      const domain = q.metadata.domain || 'Unknown';
      domainCounts[domain] = (domainCounts[domain] || 0) + 1;
    });
    Object.entries(domainCounts).forEach(([domain, count]) => {
      console.log(`  ${domain}: ${count}`);
    });
    
  } catch (error) {
    console.error('Error:', error);
  }
}

main();