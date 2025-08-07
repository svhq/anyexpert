require('dotenv').config();
const agent = require('./src/unified-agent-modular');

const stressQueries = [
  {
    name: "Healthcare Crisis Analysis",
    query: "My 68-year-old mother has been diagnosed with stage 2 pancreatic cancer. She also has type 2 diabetes and mild kidney disease. What are the latest treatment options available in 2025, including immunotherapy and targeted therapies? What clinical trials might she qualify for, and what are the survival rates with current treatments? Also explain the interaction concerns between cancer treatments and her existing conditions."
  },
  {
    name: "Complex Financial Planning",
    query: "I'm 35 years old with $150,000 in savings, earning $120,000/year. I want to retire by 55 with $3 million. I'm considering investing in a mix of index funds, real estate, and cryptocurrency. Given current market conditions and inflation rates, create a detailed 20-year investment strategy. Calculate the required monthly investments, expected returns, and tax implications. Also factor in having 2 kids and their college expenses."
  },
  {
    name: "Advanced Technical Architecture",
    query: "Design a scalable microservices architecture for a real-time trading platform that needs to handle 1 million concurrent users, process 100,000 transactions per second, and maintain 99.99% uptime. Include database choices (SQL vs NoSQL), message queue systems, caching strategies, load balancing, disaster recovery, and estimated AWS costs. Also explain how to implement circuit breakers and rate limiting."
  },
  {
    name: "Climate Change Impact Analysis",
    query: "Analyze the specific impacts of climate change on coastal cities in Southeast Asia by 2050. Include sea level rise projections, economic costs, population displacement estimates, and adaptation strategies. Compare the effectiveness and costs of different mitigation approaches like sea walls, floating cities, and managed retreat. What are the latest technologies being developed to address these challenges?"
  },
  {
    name: "Biotech Investment Research",
    query: "Evaluate the investment potential of CRISPR gene editing companies for the next 5 years. Analyze the competitive landscape, patent situations, FDA approval pipelines, and potential market sizes for different applications (cancer, genetic diseases, agriculture). Which companies have the strongest portfolios and why? Include recent clinical trial results and breakthrough developments from 2024-2025."
  },
  {
    name: "Legal & Regulatory Compliance",
    query: "Our AI startup is launching a health diagnostic app that uses machine learning to analyze medical images and provide preliminary diagnoses. What are the regulatory requirements for FDA approval, HIPAA compliance, and data privacy laws across US, EU (including AI Act), and Asia? What's our liability exposure, and what insurance do we need? Include recent 2025 regulatory changes."
  },
  {
    name: "Educational Career Pivot",
    query: "I'm a 42-year-old marketing manager wanting to transition into data science. Compare bootcamps vs master's degrees vs self-learning paths. What specific skills and certifications are most valuable in 2025? Create a 12-month learning roadmap with resources, costs, and time commitments. What's the realistic salary expectation and job market outlook? Include success stories and failure rates."
  },
  {
    name: "Quantum Computing Applications",
    query: "Explain the current state of quantum computing in 2025 and its practical applications beyond research labs. How close are we to quantum advantage in drug discovery, cryptography, and financial modeling? Which companies are leading, what are the main technical challenges remaining, and when can businesses realistically start implementing quantum solutions? Include recent breakthroughs and setbacks."
  }
];

async function runStressTests() {
  console.log("=".repeat(70));
  console.log("STRESS TEST - COMPLEX REAL-WORLD QUERIES");
  console.log("=".repeat(70));
  
  const results = [];
  
  for (let i = 0; i < stressQueries.length; i++) {
    const query = stressQueries[i];
    console.log(`\n[${i+1}/${stressQueries.length}] ${query.name}`);
    console.log("-".repeat(70));
    console.log(`Query: "${query.query.substring(0, 100)}..."`);
    
    const startTime = Date.now();
    
    try {
      const result = await agent.process(query.query);
      const duration = Date.now() - startTime;
      
      const summary = {
        query: query.name,
        success: true,
        duration: duration,
        steps: result.metadata.totalSteps,
        tools: result.metadata.toolsUsed,
        confidence: result.metadata.finalConfidence,
        answerLength: result.content.length,
        searchPerformed: result.searchPerformed || false,
        codeExecuted: result.metadata.toolsUsed.includes('code'),
        tokensUsed: result.metadata.tokensUsed || 0
      };
      
      results.push(summary);
      
      console.log(`✓ Completed in ${(duration/1000).toFixed(1)}s`);
      console.log(`  Steps: ${summary.steps}`);
      console.log(`  Tools: ${summary.tools.join(', ')}`);
      console.log(`  Answer length: ${summary.answerLength} chars`);
      console.log(`  Search: ${summary.searchPerformed ? 'Yes' : 'No'}`);
      console.log(`  Code: ${summary.codeExecuted ? 'Yes' : 'No'}`);
      
      // Show first 200 chars of answer
      console.log(`  Preview: ${result.content.substring(0, 200)}...`);
      
    } catch (error) {
      console.log(`✗ Error: ${error.message}`);
      results.push({
        query: query.name,
        success: false,
        error: error.message
      });
    }
    
    // Brief pause between queries to avoid overwhelming the system
    if (i < stressQueries.length - 1) {
      console.log("\nWaiting 2 seconds before next query...");
      await new Promise(resolve => setTimeout(resolve, 2000));
    }
  }
  
  // Summary statistics
  console.log("\n" + "=".repeat(70));
  console.log("STRESS TEST SUMMARY");
  console.log("=".repeat(70));
  
  const successful = results.filter(r => r.success);
  const failed = results.filter(r => !r.success);
  
  console.log(`\nTotal queries: ${results.length}`);
  console.log(`Successful: ${successful.length}`);
  console.log(`Failed: ${failed.length}`);
  
  if (successful.length > 0) {
    const avgDuration = successful.reduce((a, b) => a + b.duration, 0) / successful.length;
    const avgSteps = successful.reduce((a, b) => a + b.steps, 0) / successful.length;
    const searchCount = successful.filter(r => r.searchPerformed).length;
    const codeCount = successful.filter(r => r.codeExecuted).length;
    
    console.log(`\nAverage response time: ${(avgDuration/1000).toFixed(1)}s`);
    console.log(`Average steps: ${avgSteps.toFixed(1)}`);
    console.log(`Used web search: ${searchCount}/${successful.length}`);
    console.log(`Used code execution: ${codeCount}/${successful.length}`);
    
    // Tool usage breakdown
    const toolUsage = {};
    successful.forEach(r => {
      r.tools.forEach(tool => {
        toolUsage[tool] = (toolUsage[tool] || 0) + 1;
      });
    });
    
    console.log("\nTool usage:");
    Object.entries(toolUsage).forEach(([tool, count]) => {
      console.log(`  ${tool}: ${count} times`);
    });
  }
  
  // Save results to file
  const fs = require('fs');
  const timestamp = new Date().toISOString().replace(/[:.]/g, '-');
  const filename = `stress-test-results-${timestamp}.json`;
  
  fs.writeFileSync(filename, JSON.stringify({
    timestamp: new Date().toISOString(),
    queries: stressQueries.map(q => q.name),
    results: results,
    summary: {
      total: results.length,
      successful: successful.length,
      failed: failed.length,
      avgResponseTime: successful.length > 0 ? 
        (successful.reduce((a, b) => a + b.duration, 0) / successful.length / 1000).toFixed(1) + 's' : 'N/A'
    }
  }, null, 2));
  
  console.log(`\nResults saved to: ${filename}`);
  console.log("\n" + "=".repeat(70));
  console.log("STRESS TEST COMPLETE");
  console.log("=".repeat(70));
}

// Run tests
console.log("Starting stress test with complex real-world queries...\n");
runStressTests().catch(error => {
  console.error("Fatal error:", error);
  process.exit(1);
});