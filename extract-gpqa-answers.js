const fs = require('fs');

// Read the GPQA response
const data = JSON.parse(fs.readFileSync('gpqa-q9-response-2025-08-01T14-40-19-561Z.json', 'utf8'));
const answer = data.answer;

console.log('EXTRACTING KEY ANSWERS FROM GPQA Q9 RESPONSE\n');
console.log('='.repeat(80));

// Search for specific answers more carefully
console.log('Response structure analysis:');
console.log('- Total length:', answer.length, 'characters');
console.log('- Contains "negative ρ":', answer.includes('negative ρ'));
console.log('- Contains "Hammett":', answer.includes('Hammett'));
console.log('- Contains actual numbered answers:', /\b[1-5]\.\s/g.test(answer));

// Look for sections that directly answer the questions
const sections = answer.split(/### \d\.|##+ \d\./);
console.log('\nFound', sections.length, 'sections');

// Try to find specific answers by looking for question keywords
console.log('\n' + '='.repeat(80));
console.log('SEARCHING FOR SPECIFIC ANSWERS:\n');

// Look for Hammett/negative rho discussion
const hammetIndex = answer.indexOf('Hammett');
if (hammetIndex > -1) {
  console.log('3. NEGATIVE ρ VALUE (Hammett discussion):');
  console.log('Found at position:', hammetIndex);
  const excerpt = answer.substring(hammetIndex - 100, hammetIndex + 500);
  console.log(excerpt);
  console.log('\n' + '-'.repeat(40) + '\n');
}

// Look for where it actually tries to explain things
const explanationPhrases = [
  'negative ρ value indicates',
  'negative rho indicates',
  'fractional order suggests',
  'Et3N acts as',
  'Et3N role is',
  'rate-determining step is'
];

explanationPhrases.forEach(phrase => {
  const index = answer.toLowerCase().indexOf(phrase.toLowerCase());
  if (index > -1) {
    console.log(`Found "${phrase}" at position ${index}:`);
    console.log(answer.substring(index, index + 300));
    console.log('\n' + '-'.repeat(40) + '\n');
  }
});

// Check if the response was cut off or completed
console.log('Response ending analysis:');
const last500 = answer.substring(answer.length - 500);
console.log('Last 500 characters:', last500);
console.log('\nEnds with complete sentence?', /[.!?]\s*$/.test(answer));
console.log('Contains "Rate-Determining Step" at end?', last500.includes('Rate-Determining Step'));

// Count how many times certain phrases repeat (indicating a loop)
const loopPhrases = [
  'fractional orders suggest',
  'Rate-Determining Step:',
  'A more appropriate explanation'
];

loopPhrases.forEach(phrase => {
  const matches = answer.match(new RegExp(phrase, 'g'));
  if (matches) {
    console.log(`\n"${phrase}" appears ${matches.length} times`);
  }
});

// Save a truncated version focusing on unique content
const truncated = answer.substring(0, 10000) + '\n\n[... Response continues but appears to loop ...]';
fs.writeFileSync('gpqa-q9-truncated.txt', truncated);
console.log('\nSaved first 10,000 characters to gpqa-q9-truncated.txt');