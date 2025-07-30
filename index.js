require('dotenv').config();
const workflowEngine = require('./src/workflow-engine');
const logger = require('./src/utils/logger');

/**
 * Main entry point for Ask Any Expert system
 */
async function main() {
  logger.info('Ask Any Expert system starting...');
  
  // Validate environment
  const requiredEnvVars = ['OPENROUTER_API_KEY', 'SERPER_API_KEY'];
  const missingVars = requiredEnvVars.filter(v => !process.env[v]);
  
  if (missingVars.length > 0) {
    logger.error(`Missing required environment variables: ${missingVars.join(', ')}`);
    process.exit(1);
  }

  logger.info({
    model: process.env.OPENROUTER_MODEL || 'z-ai/glm-4.5-air:free',
    maxRounds: process.env.MAX_ROUNDS || 3,
    confidenceThreshold: process.env.CONF_THRESHOLD || 0.8,
    environment: process.env.NODE_ENV || 'development'
  }, 'Configuration loaded');

  // Example usage - in production this would be an API server
  if (process.argv[2]) {
    const query = process.argv.slice(2).join(' ');
    logger.info({ query }, 'Processing command line query');
    
    try {
      const result = await workflowEngine.answer(query);
      console.log('\n--- Response ---');
      console.log(result.content);
      
      if (result.sources && result.sources.length > 0) {
        console.log('\n--- Sources ---');
        result.sources.forEach(source => {
          console.log(`[${source.number}] ${source.title} - ${source.url}`);
        });
      }
      
      logger.info({
        searchPerformed: result.searchPerformed,
        rounds: result.metadata?.rounds,
        confidence: result.metadata?.finalConfidence,
        duration: result.metadata?.totalDuration
      }, 'Query completed');
      
    } catch (error) {
      logger.error({ error: error.message }, 'Query failed');
      console.error('Error:', error.message);
    }
  } else {
    console.log('Ask Any Expert - Web-Augmented Expert System');
    console.log('\nUsage: node index.js <your question>');
    console.log('Example: node index.js What are the latest AI breakthroughs?');
  }
}

// Run if called directly
if (require.main === module) {
  main().catch(error => {
    logger.error({ error: error.message, stack: error.stack }, 'Fatal error');
    process.exit(1);
  });
}

module.exports = { main };