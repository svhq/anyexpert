const workflowEngine = require('./src/workflow-engine');

async function checkQ2Details() {
  console.log('🔍 Checking Q2 Physics calculation details...\n');
  
  const question = "Two quantum states with energies E1 and E2 have a lifetime of 10^-9 sec and 10^-8 sec, respectively. We want to clearly distinguish these two energy levels. Which one of the following options could be their energy difference so that they can be clearly resolved? A) 10^-9 eV, B) 10^-8 eV, C) 10^-11 eV, D) 10^-4 eV";
  
  try {
    const response = await workflowEngine.answer(question);
    
    console.log('Math executed?', response.mathExecuted);
    console.log('Tools used?', response.executionResults ? 'Yes' : 'No');
    
    if (response.executionResults) {
      console.log('\n🔧 Code execution output:');
      response.executionResults.forEach(exec => {
        console.log(exec.result.stdout);
      });
    }
    
    // Look for specific answer patterns
    const content = response.content;
    
    console.log('\n🎯 Looking for answer patterns:');
    
    if (content.includes('10^-4') || content.includes('1e-04')) {
      console.log('✅ Found 10^-4 eV mentioned');
    }
    
    if (content.toLowerCase().includes('option d') || content.toLowerCase().includes('answer d')) {
      console.log('✅ Explicitly mentions option D');
    }
    
    if (content.includes('6.58') && content.includes('10^-7')) {
      console.log('✅ Calculated energy uncertainty ~6.58 × 10^-7 eV');
    }
    
    // Search for conclusion
    const lines = content.split('\n');
    lines.forEach((line, i) => {
      if (line.toLowerCase().includes('conclusion') || 
          line.toLowerCase().includes('therefore') ||
          line.toLowerCase().includes('answer') ||
          (line.includes('10^-4') && line.toLowerCase().includes('option'))) {
        console.log(`Conclusion line ${i}: ${line.trim()}`);
      }
    });
    
    return content.includes('10^-4');
    
  } catch (error) {
    console.error('Error:', error.message);
    return false;
  }
}

checkQ2Details();