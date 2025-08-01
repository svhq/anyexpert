const fs = require('fs');
const path = require('path');
const { ModelLoader } = require('../../models/model-loader');
const workflowEngine = require('../../src/workflow-engine');
const logger = require('../../src/utils/logger');

/**
 * Model Test Runner - Runs standardized tests across different models
 */
class ModelTestRunner {
  constructor(modelKey) {
    this.modelLoader = new ModelLoader();
    this.modelKey = modelKey;
    this.model = this.modelLoader.getModel(modelKey);
    this.testQueries = this.loadTestQueries();
    this.results = {
      model: modelKey,
      modelInfo: this.model,
      timestamp: new Date().toISOString(),
      tests: [],
      summary: {}
    };
  }

  /**
   * Load test queries from JSON file
   */
  loadTestQueries() {
    const queriesPath = path.join(__dirname, 'test-queries.json');
    const data = fs.readFileSync(queriesPath, 'utf8');
    return JSON.parse(data);
  }

  /**
   * Override environment to use specific model
   */
  setModel() {
    process.env.OPENROUTER_MODEL = this.model.id;
    console.log(`\n🔄 Switched to model: ${this.model.name} (${this.model.id})`);
  }

  /**
   * Run a single test query
   */
  async runQuery(queryInfo, category) {
    console.log(`\n📝 Testing: ${queryInfo.id}`);
    console.log(`   Query: ${queryInfo.query}`);
    
    const startTime = Date.now();
    const result = {
      id: queryInfo.id,
      category: category,
      query: queryInfo.query,
      timestamp: new Date().toISOString()
    };

    try {
      // Call the workflow engine
      const response = await workflowEngine.answer(queryInfo.query, {});
      const duration = Date.now() - startTime;
      
      // Extract key information
      result.response = {
        content: response.content,
        duration_ms: duration,
        searchPerformed: response.searchPerformed || false,
        codeExecuted: response.metadata?.codeExecuted || false,
        expertPersona: this.extractExpertPersona(response.content),
        responseLength: response.content.length
      };

      // Validate response if expected answer exists
      if (queryInfo.expectedAnswer) {
        result.validation = this.validateResponse(
          response.content, 
          queryInfo.expectedAnswer, 
          queryInfo.validateFunction
        );
      }

      // Estimate token usage (rough estimate)
      const estimatedTokens = {
        input: Math.ceil((queryInfo.query.length + 1000) / 4), // query + system prompt
        output: Math.ceil(response.content.length / 4)
      };
      
      result.cost = this.modelLoader.calculateCost(
        this.modelKey,
        estimatedTokens.input,
        estimatedTokens.output
      );

      console.log(`   ✅ Response received (${duration}ms)`);
      console.log(`   💰 Estimated cost: $${result.cost.totalCost}`);
      
      if (result.validation) {
        console.log(`   🎯 Validation: ${result.validation.passed ? 'PASSED' : 'FAILED'}`);
      }

    } catch (error) {
      result.error = error.message;
      result.response = { error: true, duration_ms: Date.now() - startTime };
      console.log(`   ❌ Error: ${error.message}`);
    }

    return result;
  }

  /**
   * Extract expert persona from response
   */
  extractExpertPersona(content) {
    const patterns = [
      /I am (Dr\.|Professor) ([^,\n]+)/i,
      /\*\*(Dr\.|Professor) ([^*]+)\*\*/,
      /^(Dr\.|Professor) ([^,\n]+)/m
    ];

    for (const pattern of patterns) {
      const match = content.match(pattern);
      if (match) {
        return match[0];
      }
    }
    return null;
  }

  /**
   * Validate response against expected answer
   */
  validateResponse(response, expected, validateFunction) {
    const validators = {
      exact_match: (r, e) => r.includes(e),
      contains: (r, e) => r.toLowerCase().includes(e.toLowerCase()),
      contains_concept: (r, e) => {
        const concepts = e.split(' ');
        return concepts.some(c => r.toLowerCase().includes(c.toLowerCase()));
      },
      has_search_results: (r) => r.includes('[1]') || r.includes('Source:'),
      contains_last_prime: (r, e) => r.includes(e),
      has_statistics: (r) => r.includes('mean') && r.includes('median'),
      contains_logic: (r, e) => r.toLowerCase().includes(e.toLowerCase()),
      contains_months: (r, e) => r.includes(e) || r.includes('month'),
      contains_speed: (r, e) => r.includes(e),
      has_expert_persona: (r) => this.extractExpertPersona(r) !== null
    };

    const validator = validators[validateFunction] || validators.contains;
    const passed = validator(response, expected);

    return {
      passed,
      expected,
      validateFunction,
      details: passed ? 'Response matches expected pattern' : 'Response does not match expected pattern'
    };
  }

  /**
   * Run all tests for a specific profile
   */
  async runTestProfile(profileName = 'quick') {
    const profile = this.testQueries.testProfiles[profileName];
    if (!profile) {
      throw new Error(`Test profile '${profileName}' not found`);
    }

    console.log(`\n🧪 Running '${profileName}' test profile with ${profile.length} queries`);
    this.setModel();

    const allQueries = [];
    
    // Collect all queries from profile
    Object.entries(this.testQueries.categories).forEach(([catKey, category]) => {
      category.queries.forEach(query => {
        if (profile.includes(query.id)) {
          allQueries.push({ ...query, category: catKey });
        }
      });
    });

    // Run each query
    for (const queryInfo of allQueries) {
      const result = await this.runQuery(queryInfo, queryInfo.category);
      this.results.tests.push(result);
      
      // Small delay between queries
      await new Promise(resolve => setTimeout(resolve, 1000));
    }

    this.generateSummary();
    return this.results;
  }

  /**
   * Run tests for specific category
   */
  async runCategory(categoryName) {
    const category = this.testQueries.categories[categoryName];
    if (!category) {
      throw new Error(`Category '${categoryName}' not found`);
    }

    console.log(`\n🧪 Running ${category.name} tests (${category.queries.length} queries)`);
    this.setModel();

    for (const queryInfo of category.queries) {
      const result = await this.runQuery(queryInfo, categoryName);
      this.results.tests.push(result);
      
      await new Promise(resolve => setTimeout(resolve, 1000));
    }

    this.generateSummary();
    return this.results;
  }

  /**
   * Generate summary statistics
   */
  generateSummary() {
    const tests = this.results.tests;
    const successful = tests.filter(t => !t.error);
    const validated = tests.filter(t => t.validation);
    const passed = validated.filter(t => t.validation.passed);

    this.results.summary = {
      totalTests: tests.length,
      successfulTests: successful.length,
      failedTests: tests.length - successful.length,
      validatedTests: validated.length,
      passedValidation: passed.length,
      failedValidation: validated.length - passed.length,
      
      averageDuration: successful.length > 0 
        ? Math.round(successful.reduce((sum, t) => sum + t.response.duration_ms, 0) / successful.length)
        : 0,
      
      totalCost: tests.reduce((sum, t) => sum + (t.cost?.totalCost || 0), 0),
      
      toolUsage: {
        search: successful.filter(t => t.response.searchPerformed).length,
        code: successful.filter(t => t.response.codeExecuted).length,
        expertPersona: successful.filter(t => t.response.expertPersona).length
      },
      
      categorySummary: this.getCategorySummary()
    };
  }

  /**
   * Get summary by category
   */
  getCategorySummary() {
    const summary = {};
    
    this.results.tests.forEach(test => {
      if (!summary[test.category]) {
        summary[test.category] = {
          total: 0,
          successful: 0,
          passed: 0,
          avgDuration: 0,
          durations: []
        };
      }
      
      summary[test.category].total++;
      if (!test.error) {
        summary[test.category].successful++;
        summary[test.category].durations.push(test.response.duration_ms);
      }
      if (test.validation?.passed) {
        summary[test.category].passed++;
      }
    });

    // Calculate average durations
    Object.values(summary).forEach(cat => {
      if (cat.durations.length > 0) {
        cat.avgDuration = Math.round(
          cat.durations.reduce((a, b) => a + b, 0) / cat.durations.length
        );
      }
      delete cat.durations;
    });

    return summary;
  }

  /**
   * Save results to file
   */
  saveResults() {
    const timestamp = new Date().toISOString().replace(/[:.]/g, '-');
    const filename = `${this.modelKey}_${timestamp}.json`;
    const filepath = path.join(__dirname, 'results', filename);
    
    fs.writeFileSync(filepath, JSON.stringify(this.results, null, 2));
    console.log(`\n💾 Results saved to: ${filename}`);
    
    return filepath;
  }

  /**
   * Print summary to console
   */
  printSummary() {
    const s = this.results.summary;
    
    console.log('\n' + '='.repeat(60));
    console.log(`📊 TEST SUMMARY - ${this.model.name}`);
    console.log('='.repeat(60));
    
    console.log(`\n✅ Success Rate: ${s.successfulTests}/${s.totalTests} (${Math.round(s.successfulTests/s.totalTests*100)}%)`);
    console.log(`🎯 Validation: ${s.passedValidation}/${s.validatedTests} passed`);
    console.log(`⏱️  Avg Duration: ${s.averageDuration}ms`);
    console.log(`💰 Total Cost: $${s.totalCost.toFixed(6)}`);
    
    console.log('\n🛠️  Tool Usage:');
    console.log(`   Search: ${s.toolUsage.search} queries`);
    console.log(`   Code: ${s.toolUsage.code} queries`);
    console.log(`   Expert: ${s.toolUsage.expertPersona} queries`);
    
    console.log('\n📈 Category Performance:');
    Object.entries(s.categorySummary).forEach(([cat, data]) => {
      console.log(`   ${cat}: ${data.successful}/${data.total} (${data.avgDuration}ms avg)`);
    });
  }
}

// Export class and convenience functions
module.exports = {
  ModelTestRunner,
  
  // Quick test function
  async testModel(modelKey, profile = 'quick') {
    const runner = new ModelTestRunner(modelKey);
    const results = await runner.runTestProfile(profile);
    runner.printSummary();
    runner.saveResults();
    return results;
  }
};

// If run directly, test the model specified in command line
if (require.main === module) {
  const modelKey = process.argv[2] || 'o4-mini-high';
  const profile = process.argv[3] || 'quick';
  
  console.log(`🚀 Testing model: ${modelKey} with profile: ${profile}`);
  
  const runner = new ModelTestRunner(modelKey);
  runner.runTestProfile(profile)
    .then(() => {
      runner.printSummary();
      runner.saveResults();
      process.exit(0);
    })
    .catch(err => {
      console.error('❌ Test failed:', err.message);
      process.exit(1);
    });
}