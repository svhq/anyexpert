const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const { Sandbox } = require('@e2b/code-interpreter');

console.log('=== KILLING ALL SANDBOXES AND TESTING FRESH ===\n');

async function killAllAndTestFresh() {
  try {
    // Step 1: List and kill ALL existing sandboxes
    console.log('Step 1: Finding ALL existing sandboxes...');
    const existingSandboxes = await Sandbox.list();
    console.log(`Found ${existingSandboxes.length} sandboxes\n`);
    
    if (existingSandboxes.length > 0) {
      console.log('Killing each sandbox:');
      for (const sandboxInfo of existingSandboxes) {
        try {
          console.log(`  Killing ${sandboxInfo.sandboxId}...`);
          const sandbox = await Sandbox.connect(sandboxInfo.sandboxId);
          await sandbox.kill();
          console.log(`  ✅ Killed ${sandboxInfo.sandboxId}`);
        } catch (error) {
          console.log(`  ❌ Failed to kill ${sandboxInfo.sandboxId}: ${error.message}`);
        }
      }
      console.log('\nAll sandboxes killed.');
      
      // Wait for cleanup
      console.log('\nWaiting 5 seconds for complete cleanup...');
      await new Promise(resolve => setTimeout(resolve, 5000));
    }
    
    // Step 2: Create brand new sandbox
    console.log('\nStep 2: Creating BRAND NEW sandbox with prod-all...');
    const newSandbox = await Sandbox.create({
      template: 'prod-all',
      apiKey: process.env.E2B_API_KEY,
      timeout: 60000
    });
    
    console.log(`\n✅ Created NEW sandbox: ${newSandbox.sandboxId}`);
    console.log('This should be a completely fresh sandbox!\n');
    
    // Step 3: Test BUILD_TAG and numpy-financial
    console.log('Step 3: Testing BUILD_TAG and numpy-financial...\n');
    
    const testCode = `
import os
import sys

print("=== NEW SANDBOX INFO ===")
print(f"Sandbox ID: {os.environ.get('E2B_SANDBOX_ID', 'NOT SET')}")
print(f"Template ID: {os.environ.get('E2B_TEMPLATE_ID', 'NOT SET')}")
print(f"BUILD_TAG: {os.environ.get('BUILD_TAG', 'NOT SET')}")
print(f"Python: {sys.version}")

print("\\n=== NUMPY-FINANCIAL CHECK ===")
try:
    import numpy_financial as npf
    print("✅ numpy-financial is INSTALLED!")
    print(f"   Version: {npf.__version__}")
    print(f"   Location: {npf.__file__}")
    
    # Test calculation
    payment = npf.pmt(0.05/12, 30*12, -300000)
    print(f"   Test: Monthly payment = $" + f"{payment:.2f}")
except ImportError as e:
    print("❌ numpy-financial is NOT installed")
    print(f"   Error: {e}")
    print("\\n   This means the template hasn't been updated yet")
    print("   Expected BUILD_TAG: npf-20250805-003130")
`;
    
    const result = await newSandbox.runCode(testCode, {
      language: 'python',
      timeoutMs: 30000
    });
    
    if (result.logs?.stdout) {
      console.log(result.logs.stdout.join(''));
    }
    
    if (result.logs?.stderr) {
      const stderr = result.logs.stderr.join('');
      if (stderr.trim()) {
        console.log('\nStderr:', stderr);
      }
    }
    
    // Cleanup
    console.log('\nCleaning up...');
    await newSandbox.kill();
    console.log('✅ Done');
    
  } catch (error) {
    console.error('Error:', error.message);
  }
}

killAllAndTestFresh().catch(console.error);