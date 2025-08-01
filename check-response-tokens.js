const fs = require('fs');

// Read the GPQA response
const data = JSON.parse(fs.readFileSync('gpqa-q9-response-2025-08-01T14-40-19-561Z.json', 'utf8'));
const answer = data.answer;

console.log('RESPONSE ANALYSIS:');
console.log('='.repeat(60));
console.log('Total characters:', answer.length);

// Rough token estimate (1 token ≈ 4 characters on average)
const estimatedTokens = Math.round(answer.length / 4);
console.log('Estimated tokens (÷4):', estimatedTokens);
console.log('Estimated tokens (÷3):', Math.round(answer.length / 3));

// Check if 16,384 tokens limit was hit
console.log('\nToken limit analysis:');
console.log('Math agent max_tokens setting: 16,384');
console.log('Estimated tokens used:', estimatedTokens);
console.log('Over limit?', estimatedTokens > 16384 ? 'YES' : 'NO');

// Check where the response might have been cut
console.log('\nResponse ending:');
console.log('Last 200 characters:', answer.substring(answer.length - 200));
console.log('Ends with complete sentence?', /[.!?]\s*$/.test(answer));
console.log('Ends mid-word?', /[a-zA-Z]$/.test(answer));

// Look for patterns indicating truncation
console.log('\nTruncation indicators:');
console.log('Contains "[truncated]"?', answer.includes('[truncated]'));
console.log('Contains "..."?', answer.includes('...'));
console.log('Ends with incomplete code block?', /```[^`]*$/.test(answer));

// Check metadata for token usage
if (data.metadata && data.metadata.tokenUsage) {
  console.log('\nToken usage from metadata:');
  console.log(JSON.stringify(data.metadata.tokenUsage, null, 2));
} else {
  console.log('\nNo token usage metadata found');
}

// Save just the ending for analysis
fs.writeFileSync('gpqa-response-ending.txt', answer.substring(answer.length - 5000));
console.log('\nSaved last 5000 characters to gpqa-response-ending.txt');