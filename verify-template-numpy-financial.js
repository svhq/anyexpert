const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const { Sandbox } = require('@e2b/code-interpreter');

console.log('=== VERIFYING NUMPY-FINANCIAL IN PROD-ALL TEMPLATE ===\n');

async function verifyInstallation() {
  let sandbox;
  
  try {
    // Wait a moment for template to propagate
    console.log('Waiting 5 seconds for template to propagate...');
    await new Promise(resolve => setTimeout(resolve, 5000));
    
    // Create sandbox using template name
    console.log('Creating sandbox with template name: prod-all');
    sandbox = await Sandbox.create({
      template: 'prod-all',
      apiKey: process.env.E2B_API_KEY,
      timeout: 60000
    });
    
    console.log(`✅ Sandbox created: ${sandbox.sandboxId}\n`);
    
    // Comprehensive check
    const checkCode = `
import os
import sys
import subprocess

print("=== SANDBOX ENVIRONMENT ===")
print(f"Python version: {sys.version}")
print(f"E2B_TEMPLATE_ID: {os.environ.get('E2B_TEMPLATE_ID', 'NOT SET')}")
print(f"BUILD_TAG: {os.environ.get('BUILD_TAG', 'NOT SET')}")

print("\\n=== CHECKING NUMPY-FINANCIAL ===")

# Method 1: Direct import
try:
    import numpy_financial as npf
    print("✅ Direct import successful!")
    print(f"   Version: {npf.__version__}")
    print(f"   Location: {npf.__file__}")
except ImportError as e:
    print(f"❌ Direct import failed: {e}")

# Method 2: Check pip list
print("\\n=== PIP LIST CHECK ===")
result = subprocess.run([sys.executable, "-m", "pip", "list"], capture_output=True, text=True)
pip_list = result.stdout
if "numpy-financial" in pip_list:
    for line in pip_list.split('\\n'):
        if "numpy-financial" in line.lower():
            print(f"✅ Found in pip list: {line.strip()}")
else:
    print("❌ NOT found in pip list")

# Method 3: pkg_resources check
print("\\n=== PKG_RESOURCES CHECK ===")
try:
    import pkg_resources
    try:
        version = pkg_resources.get_distribution("numpy-financial").version
        print(f"✅ Found via pkg_resources: numpy-financial=={version}")
    except pkg_resources.DistributionNotFound:
        print("❌ NOT found via pkg_resources")
except ImportError:
    print("⚠️  pkg_resources not available")

# Method 4: Test actual functionality
print("\\n=== FUNCTIONALITY TEST ===")
try:
    import numpy_financial as npf
    # Calculate monthly payment for $300k mortgage at 5% for 30 years
    payment = npf.pmt(0.05/12, 30*12, -300000)
    print(f"✅ Calculation successful!")
    print(f"   Monthly payment: $" + f"{payment:.2f}")
except Exception as e:
    print(f"❌ Functionality test failed: {e}")

print("\\n=== SUMMARY ===")
try:
    import numpy_financial
    print("✅ NUMPY-FINANCIAL IS INSTALLED AND WORKING!")
except:
    print("❌ NUMPY-FINANCIAL IS NOT INSTALLED")
`;
    
    const result = await sandbox.runCode(checkCode, {
      language: 'python',
      timeoutMs: 30000
    });
    
    if (result.logs?.stdout) {
      console.log(result.logs.stdout.join(''));
    }
    
    if (result.logs?.stderr) {
      console.log('\nStderr:');
      console.log(result.logs.stderr.join(''));
    }
    
    if (result.error) {
      console.log('\nExecution error:', result.error);
    }
    
  } catch (error) {
    console.error('Error creating sandbox:', error.message);
  } finally {
    if (sandbox) {
      console.log('\nCleaning up sandbox...');
      await sandbox.kill();
    }
  }
}

verifyInstallation().catch(console.error);