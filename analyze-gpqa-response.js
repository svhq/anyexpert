const fs = require('fs');

// Read the latest GPQA response
const filename = 'gpqa-q9-response-2025-08-01T14-40-19-561Z.json';
const data = JSON.parse(fs.readFileSync(filename, 'utf8'));

console.log('GPQA Q9 RESPONSE ANALYSIS\n');
console.log('='.repeat(80));

// Basic info
console.log('Response length:', data.answer.length, 'characters');
console.log('Agent used:', data.agent);
console.log('');

// Try to find the actual answers to the 5 questions
const answer = data.answer;

// Look for answers to each question
console.log('SEARCHING FOR ANSWERS TO THE 5 QUESTIONS:\n');

// Question 1: Rate-determining step
const rdsMatch = answer.match(/(?:rate[- ]?determining step|RDS)[^.]*?(?:is|likely|involves?)([^.]+\.)/i);
if (rdsMatch) {
  console.log('1. Rate-determining step:');
  console.log('   Found around character', answer.indexOf(rdsMatch[0]));
  console.log('   ', rdsMatch[0].substring(0, 200));
}

// Question 2: Mechanism
const mechMatch = answer.match(/(?:proposed mechanism|mechanism:)([^]{200,500})/i);
if (mechMatch) {
  console.log('\n2. Proposed mechanism:');
  console.log('   Found around character', answer.indexOf(mechMatch[0]));
  console.log('   [Mechanism details found but truncated for brevity]');
}

// Question 3: Negative ρ value
const rhoMatch = answer.match(/negative\s*(?:ρ|rho)[^.]*?(?:indicates?|suggests?|means?)([^.]+\.)/i);
if (rhoMatch) {
  console.log('\n3. Negative ρ value explanation:');
  console.log('   Found around character', answer.indexOf(rhoMatch[0]));
  console.log('   ', rhoMatch[0].substring(0, 200));
}

// Question 4: Fractional order in alkene
const fractionalMatch = answer.match(/(?:order.*alkene.*less.*than.*1|fractional.*order.*alkene)([^]{100,300})/i);
if (fractionalMatch) {
  console.log('\n4. Fractional alkene order:');
  console.log('   Found around character', answer.indexOf(fractionalMatch[0]));
  console.log('   ', fractionalMatch[0].substring(0, 200));
}

// Question 5: Et3N role
const et3nMatch = answer.match(/(?:Et3N|triethylamine).*?(?:role|function|acts?)([^]{100,300})/i);
if (et3nMatch) {
  console.log('\n5. Et3N role:');
  console.log('   Found around character', answer.indexOf(et3nMatch[0]));
  console.log('   ', et3nMatch[0].substring(0, 200));
}

// Look for a summary or conclusion
console.log('\n' + '='.repeat(80));
console.log('LOOKING FOR SUMMARY/CONCLUSION:\n');

// Try to find conclusion section
const conclusionMatch = answer.match(/(?:summary|conclusion|in summary|to summarize)([^]{500,2000})/i);
if (conclusionMatch) {
  console.log('Found conclusion around character', answer.indexOf(conclusionMatch[0]));
  console.log('\nConclusion excerpt:');
  console.log(conclusionMatch[0].substring(0, 1000));
}

// Check for final answer format
const finalAnswerMatch = answer.match(/(?:final answer|answer:)([^]{200,1000})/i);
if (finalAnswerMatch) {
  console.log('\nFound "Final Answer" section:');
  console.log(finalAnswerMatch[0].substring(0, 500));
}

// Look at the very end of the response
console.log('\n' + '='.repeat(80));
console.log('LAST 2000 CHARACTERS OF RESPONSE:\n');
console.log(answer.substring(answer.length - 2000));