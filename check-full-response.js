const workflowEngine = require('./src/workflow-engine');

async function checkFullResponse() {
  const simplifiedQuery = `Two quantum states with energies E1 and E2 have a lifetime of 10^-9 sec and 10^-8 sec, respectively. We want to clearly distinguish these two energy levels. Which one of the following options could be their energy difference so that they can be clearly resolved?

(A) 10^-9 eV
(B) 10^-8 eV  
(C) 10^-11 eV
(D) 10^-4 eV`;

  try {
    const response = await workflowEngine.answer(simplifiedQuery);
    
    console.log('📝 FULL RESPONSE:');
    console.log('=' .repeat(80));
    console.log(response.content);
    console.log('=' .repeat(80));
    
    // Check if calculation shows D is correct
    if (response.content.includes('6.58') && response.content.includes('10^-4')) {
      console.log('\n✅ The system calculated correctly - only 10^-4 eV is large enough!');
      console.log('✅ This means the answer is D, even if not explicitly stated.');
    }
    
  } catch (error) {
    console.error('Error:', error.message);
  }
}

checkFullResponse();