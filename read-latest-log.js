const fs = require('fs');
const path = require('path');
const readline = require('readline');

// Path to latest log
const latestLogPath = path.join(__dirname, 'realtime-ui', 'session-logs', 'latest-session.jsonl');

console.log('📖 Reading latest session log...\n');

if (!fs.existsSync(latestLogPath)) {
  console.error('❌ No log file found. Run a question through the UI first.');
  process.exit(1);
}

// Get the actual file path
let actualPath = latestLogPath;
try {
  const firstLine = fs.readFileSync(latestLogPath, 'utf8').split('\n')[0];
  const redirect = JSON.parse(firstLine);
  if (redirect.type === 'redirect' && redirect.actualFile) {
    actualPath = redirect.actualFile;
    console.log(`📁 Log file: ${path.basename(actualPath)}\n`);
  }
} catch (e) {
  // Use latest-session.jsonl directly
}

// Read line by line
const fileStream = fs.createReadStream(latestLogPath);
const rl = readline.createInterface({
  input: fileStream,
  crlfDelay: Infinity
});

let lineCount = 0;
const summary = {
  totalLines: 0,
  types: {},
  errors: [],
  agents: new Set(),
  apiCalls: 0,
  codeExecutions: 0,
  searchQueries: 0,
  duration: 0
};

rl.on('line', (line) => {
  lineCount++;
  try {
    const entry = JSON.parse(line);
    summary.totalLines++;
    
    // Count by type
    summary.types[entry.type] = (summary.types[entry.type] || 0) + 1;
    
    // Track errors
    if (entry.type === 'error' || entry.type === 'workflow_error' || entry.level === 'error') {
      summary.errors.push({
        line: lineCount,
        message: entry.message || entry.error,
        timestamp: entry.timestamp
      });
    }
    
    // Track agents
    if (entry.type === 'agent_start') {
      summary.agents.add(entry.agent);
    }
    
    // Track metrics
    if (entry.type === 'openrouter_call') summary.apiCalls++;
    if (entry.type === 'code_execution_request') summary.codeExecutions++;
    if (entry.type === 'search_request') summary.searchQueries++;
    
    // Track duration
    if (entry.type === 'question_complete') {
      summary.duration = entry.duration;
    }
    
    // Print interesting entries
    if (['error', 'workflow_error', 'search_decision', 'query_classification', 'agent_start', 'code_execution_error'].includes(entry.type)) {
      console.log(`[${entry.timestamp}] ${entry.type.toUpperCase()}:`);
      console.log(JSON.stringify(entry, null, 2));
      console.log('-'.repeat(80));
    }
    
  } catch (e) {
    console.error(`Error parsing line ${lineCount}: ${e.message}`);
  }
});

rl.on('close', () => {
  console.log('\n' + '='.repeat(80));
  console.log('📊 SESSION SUMMARY');
  console.log('='.repeat(80));
  console.log(`Total log entries: ${summary.totalLines}`);
  console.log(`Duration: ${(summary.duration / 1000).toFixed(2)}s`);
  console.log(`\nAPI Calls: ${summary.apiCalls}`);
  console.log(`Code Executions: ${summary.codeExecutions}`);
  console.log(`Search Queries: ${summary.searchQueries}`);
  console.log(`\nAgents Used: ${Array.from(summary.agents).join(', ') || 'None'}`);
  console.log(`\nErrors Found: ${summary.errors.length}`);
  
  if (summary.errors.length > 0) {
    console.log('\n❌ ERRORS:');
    summary.errors.forEach(err => {
      console.log(`  Line ${err.line}: ${err.message}`);
    });
  }
  
  console.log('\n📈 Event Types:');
  Object.entries(summary.types)
    .sort((a, b) => b[1] - a[1])
    .slice(0, 10)
    .forEach(([type, count]) => {
      console.log(`  ${type}: ${count}`);
    });
});

// Also provide a function to get specific entries
function findEntriesOfType(type) {
  const entries = [];
  const data = fs.readFileSync(latestLogPath, 'utf8');
  const lines = data.split('\n');
  
  lines.forEach((line, index) => {
    if (line.trim()) {
      try {
        const entry = JSON.parse(line);
        if (entry.type === type) {
          entries.push({ lineNumber: index + 1, ...entry });
        }
      } catch (e) {
        // Skip invalid lines
      }
    }
  });
  
  return entries;
}

// Export for use in other scripts
module.exports = { findEntriesOfType };