// Check if GLM-4.5 used search tools during MMLU testing
const fs = require('fs');

console.log('🔍 CHECKING GLM-4.5 SEARCH TOOL USAGE\n');

// Load the retest results 
const retestResults = JSON.parse(fs.readFileSync('./glm-4.5-complete-final.json', 'utf8'));

console.log('📊 SEARCH USAGE ANALYSIS:\n');

let searchUsed = 0;
let noSearchNeeded = 0;

retestResults.retest_results.forEach(result => {
  console.log(`${result.id} (${result.category}):`);
  
  // Check if search was performed (if the field exists)
  if (result.hasOwnProperty('search_performed')) {
    if (result.search_performed) {
      console.log('  🔍 SEARCH USED');
      searchUsed++;
    } else {
      console.log('  📚 NO SEARCH - Direct answer');
      noSearchNeeded++;
    }
  } else {
    console.log('  ❓ Search info not available');
  }
  
  // Look for search indicators in the response
  if (result.full_response) {
    const response = result.full_response.toLowerCase();
    
    // Look for search-related phrases
    const searchIndicators = [
      'searched for',
      'found information',
      'according to sources',
      'research shows',
      'studies indicate',
      'recent data',
      'current information'
    ];
    
    const hasSearchIndicators = searchIndicators.some(indicator => 
      response.includes(indicator)
    );
    
    if (hasSearchIndicators) {
      console.log('  🔍 Response suggests search was used');
    }
    
    // Check for expert personas (our workflow creates these)
    const expertPatterns = [
      /as (dr\.|professor|expert)/i,
      /i'm (dr\.|professor)/i,
      /\*\*(dr\.|professor).*?\*\*/i
    ];
    
    const hasExpertPersona = expertPatterns.some(pattern => 
      response.match(pattern)
    );
    
    if (hasExpertPersona) {
      console.log('  👨‍🎓 Expert persona used (workflow feature)');
    }
  }
  
  console.log('');
});

console.log('\n' + '='.repeat(60));
console.log('📈 SUMMARY:');
console.log(`Total questions analyzed: ${retestResults.retest_results.length}`);
console.log(`Questions with search used: ${searchUsed}`);
console.log(`Questions without search: ${noSearchNeeded}`);
console.log(`Search usage rate: ${Math.round(searchUsed / retestResults.retest_results.length * 100)}%`);

// Also check the original results from the first test
console.log('\n🔍 CHECKING ORIGINAL TEST RESULTS:\n');

const originalResults = JSON.parse(fs.readFileSync('./glm-4.5-fixed-final.json', 'utf8'));

let originalSearchUsed = 0;
originalResults.results.forEach(result => {
  if (result.hasOwnProperty('search_performed') && result.search_performed) {
    originalSearchUsed++;
  }
});

console.log(`Original test search usage: ${originalSearchUsed}/${originalResults.results.length}`);

// Final assessment
console.log('\n' + '='.repeat(60));
console.log('🎯 SEARCH TOOL USAGE ASSESSMENT:');

if (searchUsed === 0 && originalSearchUsed === 0) {
  console.log('❌ NO SEARCH TOOLS USED');
  console.log('GLM-4.5 relied entirely on its training data');
  console.log('This may explain lower performance on knowledge-intensive questions');
} else {
  console.log('✅ SEARCH TOOLS WERE USED');
  console.log(`Search was used on ${searchUsed + originalSearchUsed} questions`);
}

console.log('\n💡 IMPLICATIONS:');
console.log('- MMLU questions test factual knowledge from training data');
console.log('- Search tools could help with recent/specific information');
console.log('- GLM-4.5 performance reflects pure model capability without search augmentation');