const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const { Sandbox } = require('@e2b/code-interpreter');

console.log('=== CLEARING ALL E2B SANDBOXES AND TESTING FRESH ===\n');

async function clearAndTest() {
  try {
    // Step 1: List and kill all existing sandboxes
    console.log('Step 1: Clearing all existing sandboxes...');
    const existingSandboxes = await Sandbox.list();
    console.log(`Found ${existingSandboxes.length} existing sandboxes`);
    
    for (const sandboxInfo of existingSandboxes) {
      try {
        console.log(`  Killing sandbox ${sandboxInfo.sandboxId}...`);
        const sandbox = await Sandbox.connect(sandboxInfo.sandboxId);
        await sandbox.kill();
      } catch (error) {
        console.log(`  Failed to kill ${sandboxInfo.sandboxId}: ${error.message}`);
      }
    }
    
    console.log('\nAll sandboxes cleared.\n');
    
    // Step 2: Wait a moment
    console.log('Step 2: Waiting 3 seconds for cleanup...');
    await new Promise(resolve => setTimeout(resolve, 3000));
    
    // Step 3: Create a fresh sandbox
    console.log('\nStep 3: Creating fresh sandbox...');
    const sandbox = await Sandbox.create({
      template: 'prod-all',
      apiKey: process.env.E2B_API_KEY
    });
    
    console.log(`Created fresh sandbox: ${sandbox.sandboxId}`);
    
    // Step 4: Test numpy_financial
    console.log('\nStep 4: Testing numpy_financial...\n');
    
    const testCode = `
import os
import sys
print(f"Python: {sys.version}")
print(f"BUILD_TAG: {os.environ.get('BUILD_TAG', 'missing')}")

# Test numpy_financial
try:
    import numpy_financial as npf
    print("✅ numpy_financial is installed!")
    print(f"Version: {npf.__version__}")
    
    # Quick test
    payment = npf.pmt(0.05/12, 30*12, -300000)
    print(f"\\nTest calculation - Monthly payment: $" + f"{payment:.2f}")
except ImportError as e:
    print("❌ numpy_financial is NOT installed")
    print(f"Error: {e}")
    
    # Try auto-install
    print("\\nAttempting auto-install...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "numpy-financial", "--quiet"])
    import numpy_financial as npf
    print("✅ Auto-install successful!")
`;
    
    const result = await sandbox.runCode(testCode, { 
      language: 'python',
      timeoutMs: 30000 
    });
    
    if (result.error) {
      console.log('Execution error:', result.error);
    }
    
    if (result.logs?.stdout) {
      console.log('Output:');
      console.log(result.logs.stdout.join(''));
    }
    
    if (result.logs?.stderr) {
      console.log('\nErrors:');
      console.log(result.logs.stderr.join(''));
    }
    
    // Cleanup
    console.log('\nCleaning up...');
    await sandbox.kill();
    
  } catch (error) {
    console.error('Error:', error.message);
  }
}

clearAndTest().catch(console.error);