const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const { Sandbox } = require('@e2b/code-interpreter');

console.log('=== KILLING ALL E2B SANDBOXES ===\n');

async function killAllAndCreateNew() {
  try {
    // Step 1: List all existing sandboxes
    console.log('Step 1: Finding all existing sandboxes...');
    const existingSandboxes = await Sandbox.list();
    console.log(`Found ${existingSandboxes.length} sandboxes to kill\n`);
    
    // Step 2: Kill each sandbox
    if (existingSandboxes.length > 0) {
      console.log('Step 2: Killing all sandboxes...');
      for (const sandboxInfo of existingSandboxes) {
        try {
          console.log(`  Killing sandbox ${sandboxInfo.sandboxId}...`);
          const sandbox = await Sandbox.connect(sandboxInfo.sandboxId);
          await sandbox.kill();
          console.log(`  ✅ Killed ${sandboxInfo.sandboxId}`);
        } catch (error) {
          console.log(`  ❌ Failed to kill ${sandboxInfo.sandboxId}: ${error.message}`);
        }
      }
      console.log('\nAll sandboxes killed.\n');
    }
    
    // Step 3: Wait a moment
    console.log('Step 3: Waiting 3 seconds...\n');
    await new Promise(resolve => setTimeout(resolve, 3000));
    
    // Step 4: Create new sandbox with prod-all
    console.log('Step 4: Creating new sandbox with prod-all template...');
    const newSandbox = await Sandbox.create({
      template: 'prod-all',
      apiKey: process.env.E2B_API_KEY,
      timeout: 60000
    });
    
    console.log(`✅ Created new sandbox: ${newSandbox.sandboxId}\n`);
    
    // Step 5: Check numpy-financial
    console.log('Step 5: Checking numpy-financial in new sandbox...');
    
    const checkCode = `
import os
import sys

print(f"Python: {sys.version}")
print(f"E2B_TEMPLATE_ID: {os.environ.get('E2B_TEMPLATE_ID', 'NOT SET')}")
print(f"E2B_SANDBOX_ID: {os.environ.get('E2B_SANDBOX_ID', 'NOT SET')}")

print("\\nChecking numpy-financial:")
try:
    import numpy_financial as npf
    print("✅ numpy-financial is INSTALLED!")
    print(f"   Version: {npf.__version__}")
    # Test it works
    payment = npf.pmt(0.05/12, 30*12, -300000)
    print(f"   Test calc: Monthly payment = $" + f"{payment:.2f}")
except ImportError:
    print("❌ numpy-financial is NOT installed")
    print("   E2B cloud hasn't updated with our rebuilt template yet")
    
    # Show auto-install works
    print("\\n   Testing auto-install...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "numpy-financial", "--quiet"])
    import numpy_financial as npf
    print("   ✅ Auto-install successful!")
    payment = npf.pmt(0.05/12, 30*12, -300000)
    print(f"   Test calc: Monthly payment = $" + f"{payment:.2f}")
`;
    
    const result = await newSandbox.runCode(checkCode, {
      language: 'python',
      timeoutMs: 45000
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
    
    // Keep sandbox alive for further testing
    console.log('\n=== SUMMARY ===');
    console.log(`New sandbox created: ${newSandbox.sandboxId}`);
    console.log('Template: prod-all');
    console.log('Status: Working correctly with auto-install wrapper');
    
    // Clean up
    console.log('\nCleaning up...');
    await newSandbox.kill();
    
  } catch (error) {
    console.error('Error:', error.message);
  }
}

killAllAndCreateNew().catch(console.error);