const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const e2bManager = require('./src/e2b-manager');

async function checkE2B() {
  console.log('=== E2B FUNCTIONALITY CHECK ===\n');
  
  try {
    // Test 1: Simple calculation
    console.log('Test 1: Basic Python...');
    const test1 = await e2bManager.executeWithFallback('print("Hello from E2B!")\nprint(10 * 20)', {
      requestId: 'check-1'
    });
    console.log('Success:', test1.success);
    console.log('Output:', test1.stdout);
    console.log('Errors:', test1.stderr || 'None');
    
    // Test 2: Scientific library
    console.log('\nTest 2: NumPy calculation...');
    const test2 = await e2bManager.executeWithFallback(
      'import numpy as np\narr = np.array([1, 2, 3, 4, 5])\nprint(f"Mean: {arr.mean()}")\nprint(f"Sum: {arr.sum()}")',
      { requestId: 'check-2' }
    );
    console.log('Success:', test2.success);
    console.log('Output:', test2.stdout);
    
    // Show metrics
    console.log('\n=== E2B Metrics ===');
    const metrics = e2bManager.getMetrics();
    console.log('Total executions:', metrics.totalExecutions);
    console.log('Failures:', metrics.failures);
    console.log('Average execution time:', metrics.avgExecutionTime, 'ms');
    
  } catch (error) {
    console.error('\nError occurred:', error.message);
    console.error('Stack:', error.stack);
  } finally {
    console.log('\nShutting down E2B manager...');
    await e2bManager.shutdown();
    console.log('Done.');
  }
}

checkE2B();