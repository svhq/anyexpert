const GPQATestingFramework = require('./gpqa-testing-framework');

async function demoFramework() {
  console.log('🎯 DEMONSTRATING GPQA TESTING FRAMEWORK');
  console.log('=====================================');
  console.log('This framework:');
  console.log('1. ✅ Sends unconstrained queries to the model');
  console.log('2. ✅ Logs complete responses without any parsing');
  console.log('3. ✅ Allows manual review of each response');
  console.log('4. ✅ Records manual assessments for final scoring\n');
  
  const framework = new GPQATestingFramework();
  
  // Test just one question (Q2 - Physics)
  console.log('Testing Q2 (Physics) as demonstration...\n');
  
  try {
    const result = await framework.testQuestion(1); // Index 1 = Q2
    
    console.log('\n🎉 DEMO COMPLETE!');
    console.log('\n📋 WHAT HAPPENED:');
    console.log('✅ Model received the full question with options');
    console.log('✅ Model gave completely unconstrained response');
    console.log('✅ Full response was logged to JSON file');
    console.log('✅ No automatic answer extraction was attempted');
    console.log('\n📖 MANUAL REVIEW:');
    console.log('Now YOU can read the response and manually determine:');
    console.log('- What is the final answer?');
    console.log('- Is it correct?');
    console.log('- How confident are you?');
    
    // Display the response for manual review
    framework.displayResponseForManualReview(result);
    
  } catch (error) {
    console.error('Demo failed:', error);
  }
}

demoFramework();