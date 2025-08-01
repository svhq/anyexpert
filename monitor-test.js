// Monitor GLM-4.5 test progress
const fs = require('fs');

function checkProgress() {
  const resultsFile = 'glm-4.5-fixed-results.json';
  
  if (!fs.existsSync(resultsFile)) {
    console.log('❌ Results file not found');
    return;
  }
  
  try {
    const data = JSON.parse(fs.readFileSync(resultsFile, 'utf8'));
    
    console.log(`📊 GLM-4.5 Test Progress: ${data.progress || `${data.results.length}/29`}`);
    
    if (data.summary) {
      console.log('\n🏁 TEST COMPLETED!');
      console.log(`✅ Correct: ${data.summary.correct}/${data.summary.total} (${data.summary.accuracy}%)`);
      console.log(`❌ Incorrect: ${data.summary.incorrect}`);
      console.log(`⚠️ Errors: ${data.summary.errors}`);
      
      // Compare with other models
      console.log('\n📈 MODEL COMPARISON:');
      console.log('o4-mini-high: 29/29 (100%)');
      console.log('GLM-4.5-air: 24/29 (82.8%)');
      console.log(`z-ai/glm-4.5: ${data.summary.correct}/${data.summary.total} (${data.summary.accuracy}%)`);
      
    } else if (data.current_accuracy) {
      console.log(`📊 Current accuracy: ${data.current_accuracy}%`);
      console.log(`⏱️ Last updated: ${new Date(data.timestamp).toLocaleTimeString()}`);
      
      // Show recent questions
      const recent = data.results.slice(-3);
      console.log('\n📝 Recent results:');
      recent.forEach(r => {
        const status = r.is_correct ? '✅' : '❌';
        console.log(`  ${status} ${r.category}: ${r.selected_answer || 'null'} (correct: ${r.correct_answer})`);
      });
    }
    
    return data.summary ? true : false;
    
  } catch (error) {
    console.log('❌ Error reading results:', error.message);
  }
}

// Check once, or keep monitoring
const completed = checkProgress();

if (!completed && process.argv.includes('--watch')) {
  console.log('\n👀 Monitoring test progress... (Ctrl+C to stop)');
  const interval = setInterval(() => {
    console.clear();
    const done = checkProgress();
    if (done) {
      clearInterval(interval);
      console.log('\n✨ Test monitoring complete!');
    }
  }, 10000); // Check every 10 seconds
}