const fs = require('fs');
const path = require('path');

// Define the queries
const queries = [
  {
    id: 1,
    category: "Medical/Health",
    question: "I've been experiencing persistent fatigue for 3 months despite sleeping 8 hours nightly. I'm 35, exercise regularly, and eat well. My doctor ran basic blood tests that came back normal. What other causes should I investigate, and what specific tests should I request?"
  },
  {
    id: 2,
    category: "Technology/Programming", 
    question: "I'm building a React app that needs to handle real-time updates for 10,000+ concurrent users. Should I use WebSockets, Server-Sent Events, or long polling? What are the trade-offs and which cloud infrastructure would scale best?"
  },
  {
    id: 3,
    category: "Finance/Investment",
    question: "I'm 28 with $50k saved and want to start investing. My goal is retirement by 55. Should I prioritize maxing out 401k, opening a Roth IRA, or investing in index funds? How should I balance these with current inflation rates?"
  },
  {
    id: 4,
    category: "Science/Environment",
    question: "My garden soil pH is 7.8 and I want to grow blueberries which need acidic soil. What's the safest way to lower pH without harming existing plants nearby? How long will it take to see results?"
  },
  {
    id: 5,
    category: "Legal/Business",
    question: "I'm a freelance designer wanting to protect my work. Should I register copyrights for each project, use watermarks, or rely on automatic copyright? What's the most cost-effective approach for someone doing 20+ projects yearly?"
  },
  {
    id: 6,
    category: "Coding/Software",
    question: "I have a Python list of 1 million integers and need to find all pairs that sum to a target value. My current nested loop solution takes 30+ seconds. What's the most efficient algorithm and implementation? Show me the code with time complexity analysis."
  },
  {
    id: 7,
    category: "Education/Career",
    question: "I'm a 40-year-old marketing manager considering a career switch to data science. I have basic Excel skills but no programming experience. Is this realistic? What's the most efficient learning path and timeline to become job-ready?"
  },
  {
    id: 8,
    category: "Home/DIY",
    question: "My bathroom has black mold growing on the ceiling despite using the exhaust fan. The area is about 2x3 feet. Is this dangerous? Can I remove it myself safely or do I need professionals? What's the long-term solution to prevent regrowth?"
  }
];

// Get list of grok files
const grokFiles = fs.readdirSync(__dirname)
  .filter(f => f.includes('real-world-q') && f.includes('x-ai-grok-3-mini') && f.endsWith('.txt'))
  .sort((a, b) => {
    const numA = parseInt(a.match(/q(\d+)/)[1]);
    const numB = parseInt(b.match(/q(\d+)/)[1]);
    return numA - numB;
  });

// Get list of gemini files  
const geminiFiles = fs.readdirSync(__dirname)
  .filter(f => f.includes('real-world-q') && (f.includes('gemini-flash-lite') || f.includes('google-gemini')) && f.endsWith('.txt'))
  .sort((a, b) => {
    const numA = parseInt(a.match(/q(\d+)/)[1]);
    const numB = parseInt(b.match(/q(\d+)/)[1]);
    return numA - numB;
  });

// Compile Grok responses
let grokContent = `# All Grok-3-mini Expert Responses

This document contains all 8 expert query responses from the x-ai/grok-3-mini model.

---
`;

// Use the first 8 grok files (in case there are duplicates)
const uniqueGrokFiles = [];
for (let i = 1; i <= 8; i++) {
  const file = grokFiles.find(f => f.includes(`q${i}-`));
  if (file) uniqueGrokFiles.push(file);
}

uniqueGrokFiles.forEach((file, index) => {
  const query = queries[index];
  const content = fs.readFileSync(path.join(__dirname, file), 'utf8');
  
  grokContent += `
## Query ${query.id}: ${query.category}

**Question:** ${query.question}

### Full Response:

${content}

---
`;
});

// Compile Gemini responses
let geminiContent = `# All Gemini-2.5-flash-lite Expert Responses

This document contains all 8 expert query responses from the google/gemini-2.5-flash-lite model.

---
`;

// Use the first 8 gemini files
const uniqueGeminiFiles = [];
for (let i = 1; i <= 8; i++) {
  const file = geminiFiles.find(f => f.includes(`q${i}-`));
  if (file) uniqueGeminiFiles.push(file);
}

uniqueGeminiFiles.forEach((file, index) => {
  const query = queries[index];
  const content = fs.readFileSync(path.join(__dirname, file), 'utf8');
  
  geminiContent += `
## Query ${query.id}: ${query.category}

**Question:** ${query.question}

### Full Response:

${content}

---
`;
});

// Write the compiled files
fs.writeFileSync('ALL-GROK-3-MINI-RESPONSES-COMPLETE.md', grokContent);
fs.writeFileSync('ALL-GEMINI-FLASH-LITE-RESPONSES-COMPLETE.md', geminiContent);

console.log('✓ Created ALL-GROK-3-MINI-RESPONSES-COMPLETE.md');
console.log('✓ Created ALL-GEMINI-FLASH-LITE-RESPONSES-COMPLETE.md');
console.log(`\nCompiled ${uniqueGrokFiles.length} Grok responses`);
console.log(`Compiled ${uniqueGeminiFiles.length} Gemini responses`);