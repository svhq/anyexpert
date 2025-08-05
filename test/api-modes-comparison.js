require('dotenv').config();
const fs = require('fs').promises;
const path = require('path');

/**
 * API Modes Comparison Testing Framework
 * Tests all three API modes (none, search, full) with identical queries
 * Compares outputs for consistency and validates functionality
 */

// Test queries covering different scenarios
const TEST_QUERIES = [
  {
    id: 'math-basic',
    query: 'What is 15% of 250?',
    expectedTools: { none: [], search: [], full: ['code'] },
    category: 'mathematics'
  },
  {
    id: 'math-complex',
    query: 'Calculate the compound interest on $5000 at 4% annual rate for 10 years',
    expectedTools: { none: [], search: [], full: ['code'] },
    category: 'mathematics'
  },
  {
    id: 'current-events',
    query: 'What are the latest developments in AI technology in 2025?',
    expectedTools: { none: [], search: ['search_web'], full: ['search_web'] },
    category: 'current-events'
  },
  {
    id: 'factual-knowledge',
    query: 'What is the capital of Australia and when was it established?',
    expectedTools: { none: [], search: [], full: [] },
    category: 'knowledge'
  },
  {
    id: 'data-analysis',
    query: 'Analyze this dataset: temperatures [23, 25, 27, 24, 26, 28, 30] and find mean, median, std dev',
    expectedTools: { none: [], search: [], full: ['code'] },
    category: 'data-analysis'
  }
];

class ApiModesTester {
  constructor() {
    this.originalEnv = { ...process.env };
    this.results = [];
  }

  /**
   * Run comprehensive testing across all API modes
   */
  async runFullTest() {
    console.log('🚀 Starting API Modes Comparison Test\n');
    console.log(`Testing ${TEST_QUERIES.length} queries across 3 API modes\n`);

    const testResults = {
      timestamp: new Date().toISOString(),
      totalQueries: TEST_QUERIES.length,
      modes: ['none', 'search', 'full'],
      results: [],
      summary: {
        none: { passed: 0, failed: 0, errors: [] },
        search: { passed: 0, failed: 0, errors: [] },
        full: { passed: 0, failed: 0, errors: [] }
      }
    };

    // Test each query across all modes
    for (const testQuery of TEST_QUERIES) {
      console.log(`\n📝 Testing Query: "${testQuery.query}"`);
      console.log(`   Category: ${testQuery.category}`);
      console.log(`   ID: ${testQuery.id}\n`);

      const queryResults = {
        queryId: testQuery.id,
        query: testQuery.query,
        category: testQuery.category,
        modes: {}
      };

      // Test each API mode
      for (const mode of ['none', 'search', 'full']) {
        console.log(`   Testing mode: ${mode}`);
        
        try {
          const result = await this.testApiMode(mode, testQuery);
          queryResults.modes[mode] = result;
          
          if (result.success) {
            testResults.summary[mode].passed++;
            console.log(`   ✅ ${mode} mode: SUCCESS`);
          } else {
            testResults.summary[mode].failed++;
            testResults.summary[mode].errors.push(result.error);
            console.log(`   ❌ ${mode} mode: FAILED - ${result.error}`);
          }
        } catch (error) {
          console.log(`   💥 ${mode} mode: ERROR - ${error.message}`);
          queryResults.modes[mode] = {
            success: false,
            error: error.message,
            duration: 0
          };
          testResults.summary[mode].failed++;
          testResults.summary[mode].errors.push(error.message);
        }
      }

      testResults.results.push(queryResults);
    }

    // Generate comparison report
    await this.generateComparisonReport(testResults);
    
    console.log('\n📊 Test Summary:');
    for (const mode of ['none', 'search', 'full']) {
      const summary = testResults.summary[mode];
      console.log(`   ${mode.toUpperCase()}: ${summary.passed}/${TEST_QUERIES.length} passed`);
    }

    return testResults;
  }

  /**
   * Test a specific API mode with a query
   */
  async testApiMode(mode, testQuery) {
    const startTime = Date.now();
    
    // Set environment for the specific mode
    process.env.API_MODE = mode;
    process.env.USE_MODULAR_SYSTEM = 'true';
    
    // Clear require cache to ensure fresh configuration loading
    delete require.cache[require.resolve('../config/tool-config')];
    delete require.cache[require.resolve('../src/unified-agent-modular')];
    delete require.cache[require.resolve('../src/workflow-engine-modular')];
    
    try {
      const ModularWorkflowEngine = require('../src/workflow-engine-modular');
      const engine = new ModularWorkflowEngine({});
      
      // Validate configuration first
      const validation = await engine.validateConfiguration();
      if (!validation.valid) {
        return {
          success: false,
          error: `Configuration invalid: ${validation.issues.join(', ')}`,
          duration: Date.now() - startTime
        };
      }
      
      // Execute the query
      const response = await engine.answer(testQuery.query, { userId: 'api-test' });
      
      // Validate response structure
      if (!response.content || typeof response.content !== 'string') {
        return {
          success: false,
          error: 'Response missing content',
          duration: Date.now() - startTime
        };
      }
      
      // Validate expected tools were used
      const expectedTools = testQuery.expectedTools[mode] || [];
      const actualTools = response.metadata?.toolsUsed || [];
      
      const toolsMatch = this.validateToolUsage(expectedTools, actualTools);
      
      return {
        success: toolsMatch.valid,
        error: toolsMatch.valid ? null : toolsMatch.error,
        duration: Date.now() - startTime,
        response: {
          content: response.content.substring(0, 500) + '...', // Truncate for storage
          toolsUsed: actualTools,
          apiMode: response.apiMode,
          capabilities: response.capabilities
        },
        toolValidation: toolsMatch
      };
      
    } catch (error) {
      return {
        success: false,
        error: error.message,
        duration: Date.now() - startTime
      };
    } finally {
      // Restore original environment
      Object.assign(process.env, this.originalEnv);
    }
  }

  /**
   * Validate that expected tools were used correctly
   */
  validateToolUsage(expectedTools, actualTools) {
    // Filter out "reason" as it's not a real tool
    const realTools = actualTools.filter(tool => tool !== 'reason');
    
    // For "none" mode, no real tools should be used
    if (expectedTools.length === 0) {
      if (realTools.length > 0) {
        return {
          valid: false,
          error: `Expected no tools but got: ${realTools.join(', ')}`
        };
      }
      return { valid: true };
    }
    
    // Check if expected tools were used
    for (const expectedTool of expectedTools) {
      if (!realTools.includes(expectedTool)) {
        return {
          valid: false,
          error: `Expected tool '${expectedTool}' was not used`
        };
      }
    }
    
    return { valid: true };
  }

  /**
   * Generate comprehensive comparison report
   */
  async generateComparisonReport(testResults) {
    const reportPath = path.join(__dirname, '..', 'api-modes-test-results.json');
    await fs.writeFile(reportPath, JSON.stringify(testResults, null, 2));
    
    // Generate markdown report
    const markdownReport = this.generateMarkdownReport(testResults);
    const markdownPath = path.join(__dirname, '..', 'API-MODES-TEST-REPORT.md');
    await fs.writeFile(markdownPath, markdownReport);
    
    console.log(`\n📄 Detailed results saved to: ${reportPath}`);
    console.log(`📄 Summary report saved to: ${markdownPath}`);
  }

  /**
   * Generate markdown summary report
   */
  generateMarkdownReport(testResults) {
    let markdown = `# API Modes Comparison Test Report\n\n`;
    markdown += `**Test Date:** ${testResults.timestamp}\n`;
    markdown += `**Total Queries:** ${testResults.totalQueries}\n`;
    markdown += `**API Modes:** ${testResults.modes.join(', ')}\n\n`;

    // Summary table
    markdown += `## Summary\n\n`;
    markdown += `| API Mode | Passed | Failed | Success Rate |\n`;
    markdown += `|----------|--------|--------|-------------|\n`;
    
    for (const mode of testResults.modes) {
      const summary = testResults.summary[mode];
      const total = summary.passed + summary.failed;
      const successRate = total > 0 ? ((summary.passed / total) * 100).toFixed(1) : '0.0';
      markdown += `| ${mode} | ${summary.passed} | ${summary.failed} | ${successRate}% |\n`;
    }

    // Detailed results
    markdown += `\n## Detailed Results\n\n`;
    
    for (const result of testResults.results) {
      markdown += `### ${result.queryId}: ${result.category}\n\n`;
      markdown += `**Query:** "${result.query}"\n\n`;
      
      for (const mode of testResults.modes) {
        const modeResult = result.modes[mode];
        const status = modeResult.success ? '✅ PASSED' : '❌ FAILED';
        markdown += `- **${mode.toUpperCase()}:** ${status}`;
        if (!modeResult.success) {
          markdown += ` - ${modeResult.error}`;
        }
        markdown += ` (${modeResult.duration}ms)\n`;
      }
      markdown += `\n`;
    }

    // Error analysis
    markdown += `## Error Analysis\n\n`;
    for (const mode of testResults.modes) {
      const errors = testResults.summary[mode].errors;
      if (errors.length > 0) {
        markdown += `### ${mode.toUpperCase()} Mode Errors\n\n`;
        errors.forEach((error, i) => {
          markdown += `${i + 1}. ${error}\n`;
        });
        markdown += `\n`;
      }
    }

    return markdown;
  }

  /**
   * Test a single query with side-by-side comparison
   */
  async testSingleQuery(query, modes = ['none', 'search', 'full']) {
    console.log(`🔍 Testing single query across modes: ${modes.join(', ')}\n`);
    console.log(`Query: "${query}"\n`);

    const results = {};
    
    for (const mode of modes) {
      console.log(`Testing ${mode} mode...`);
      const testQuery = { query, expectedTools: { [mode]: [] } };
      results[mode] = await this.testApiMode(mode, testQuery);
      console.log(`${mode}: ${results[mode].success ? 'SUCCESS' : 'FAILED'}\n`);
    }

    return results;
  }
}

// Export for use as module
module.exports = ApiModesTester;

// Run tests if called directly
if (require.main === module) {
  const tester = new ApiModesTester();
  
  if (process.argv[2] === '--single') {
    const query = process.argv[3] || 'What is 2 + 2?';
    tester.testSingleQuery(query).then(results => {
      console.log('Results:', JSON.stringify(results, null, 2));
    }).catch(console.error);
  } else {
    tester.runFullTest().then(results => {
      process.exit(results.summary.none.failed + results.summary.search.failed + results.summary.full.failed > 0 ? 1 : 0);
    }).catch(error => {
      console.error('Test framework error:', error);
      process.exit(1);
    });
  }
}