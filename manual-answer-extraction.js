const workflowEngine = require('./src/workflow-engine');
const gpqaData = require('./gpqa-10-questions-clean.json');

async function manuallyExtractAnswers() {
  console.log('🔍 MANUAL ANSWER EXTRACTION FOR GPQA QUESTIONS 1-5\n');
  
  const questions = gpqaData.questions.slice(0, 5);
  const results = [];
  
  for (let i = 0; i < questions.length; i++) {
    const q = questions[i];
    console.log(`\n${'='.repeat(60)}`);
    console.log(`QUESTION ${i + 1}: ${q.category}`);
    console.log(`Expected Answer: ${q.correct_answer}`);
    console.log(`${'='.repeat(60)}`);
    
    try {
      const response = await workflowEngine.answer(q.question);
      const content = response.content || '';
      
      console.log(`\n📝 FULL RESPONSE (${content.length} chars):`);
      console.log(content);
      
      // Manual analysis
      console.log(`\n🔍 MANUAL ANALYSIS:`);
      console.log(`Looking for mentions of options A, B, C, D...`);
      
      // Check each option
      ['A', 'B', 'C', 'D'].forEach(opt => {
        const regex1 = new RegExp(`\\b${opt}\\b.*(?:correct|answer|solution|right)`, 'i');
        const regex2 = new RegExp(`(?:correct|answer|solution|option).*\\b${opt}\\b`, 'i');
        const regex3 = new RegExp(`option\\s*${opt}`, 'i');
        
        if (regex1.test(content) || regex2.test(content) || regex3.test(content)) {
          console.log(`   ✓ Found mention of option ${opt}`);
          
          // Get context
          const lines = content.split('\n');
          lines.forEach((line, idx) => {
            if (line.toLowerCase().includes(opt.toLowerCase()) && 
                (line.toLowerCase().includes('answer') || 
                 line.toLowerCase().includes('correct') || 
                 line.toLowerCase().includes('option') ||
                 line.toLowerCase().includes('solution'))) {
              console.log(`     Line ${idx}: ${line.trim()}`);
            }
          });
        }
      });
      
      // Look for specific content based on the question
      if (i === 0) { // Q1 - Molecular Biology
        console.log(`\n🧬 Q1 Analysis (antisense therapy):`);
        if (content.toLowerCase().includes('r-loop')) {
          console.log(`   ✓ Mentions R-loops`);
        }
        if (content.toLowerCase().includes('not involved') || content.toLowerCase().includes('not required')) {
          console.log(`   ✓ Discusses what's "not involved"`);
        }
      }
      
      if (i === 1) { // Q2 - Physics  
        console.log(`\n⚛️ Q2 Analysis (quantum uncertainty):`);
        if (content.includes('10^-4') || content.includes('1e-04') || content.includes('10⁻⁴')) {
          console.log(`   ✓ Mentions 10^-4 eV (option D)`);
        }
        if (content.toLowerCase().includes('option d') || content.toLowerCase().includes('answer d')) {
          console.log(`   ✓ Explicitly mentions option D`);
        }
      }
      
      if (i === 2) { // Q3 - Organic Chemistry
        console.log(`\n🧪 Q3 Analysis (carbon counting):`);
        ['11', '10', '12', '14'].forEach(num => {
          if (content.includes(`${num} carbon`) || content.includes(`C${num}`) || content.includes(`${num}C`)) {
            console.log(`   ✓ Mentions ${num} carbons`);
          }
        });
      }
      
      if (i === 3) { // Q4 - Molecular Biology
        console.log(`\n🔬 Q4 Analysis (ChIP-seq methods):`);
        if (content.toLowerCase().includes('chip-seq')) {
          console.log(`   ✓ Mentions ChIP-seq`);
        }
        if (content.toLowerCase().includes('chromosome conformation')) {
          console.log(`   ✓ Mentions chromosome conformation capture`);
        }
        if (content.toLowerCase().includes('qrt-pcr') || content.toLowerCase().includes('qrt pcr')) {
          console.log(`   ✓ Mentions qRT-PCR`);
        }
      }
      
      if (i === 4) { // Q5 - Quantum Mechanics
        console.log(`\n🌊 Q5 Analysis (expectation value):`);
        ['-0.7', '1.65', '0.85', '-1.4'].forEach(val => {
          if (content.includes(val)) {
            console.log(`   ✓ Contains value ${val}`);
          }
        });
      }
      
      console.log(`\n📊 MANUAL VERDICT: [Please analyze the above content]`);
      
      // Pause between questions
      if (i < questions.length - 1) {
        console.log(`\n⏸️ Press Enter to continue to next question...`);
        // Small delay instead of interactive pause
        await new Promise(resolve => setTimeout(resolve, 2000));
      }
      
    } catch (error) {
      console.error(`❌ Error processing question ${i + 1}:`, error.message);
    }
  }
}

manuallyExtractAnswers().then(() => {
  console.log('\n🏁 Manual extraction completed!');
  process.exit(0);
}).catch(error => {
  console.error('💥 Error:', error);
  process.exit(1);
});