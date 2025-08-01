const workflowEngine = require('./src/workflow-engine');

async function quickTest() {
  console.log('Testing Q2 to see actual response...');
  
  const question = "Two quantum states with energies E1 and E2 have a lifetime of 10^-9 sec and 10^-8 sec, respectively. We want to clearly distinguish these two energy levels. Which one of the following options could be their energy difference so that they can be clearly resolved? A) 10^-9 eV, B) 10^-8 eV, C) 10^-11 eV, D) 10^-4 eV";
  
  try {
    const response = await workflowEngine.answer(question);
    
    console.log('Response type:', typeof response);
    console.log('Response keys:', Object.keys(response));
    console.log('Content length:', response.content ? response.content.length : 'No content');
    console.log('First 500 chars:', response.content ? response.content.substring(0, 500) : 'No content');
    
    // Check for "D" in response
    if (response.content && response.content.includes('10^-4') || response.content.includes('D')) {
      const lines = response.content.split('\n');
      lines.forEach((line, i) => {
        if (line.includes('10^-4') || line.includes(' D') || line.toLowerCase().includes('answer')) {
          console.log(`Line ${i}: ${line}`);
        }
      });
    }
    
  } catch (error) {
    console.error('Error:', error.message);
  }
}

quickTest();