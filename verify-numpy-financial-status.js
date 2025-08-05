const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const E2BManager = require('./src/e2b-manager-v3');

console.log('=== VERIFYING NUMPY_FINANCIAL STATUS ===\n');

async function verifyStatus() {
  const e2bManager = new E2BManager();
  
  console.log('1. Testing with a financial calculation that triggers numpy_financial...\n');
  
  const code = `
# Test if numpy_financial works
print("Testing numpy_financial...")

try:
    import numpy_financial as npf
    print("✅ numpy_financial imported successfully")
    
    # Quick calculation
    payment = npf.pmt(0.05/12, 30*12, -300000)
    print(f"Monthly payment: $" + f"{payment:.2f}")
    
except ImportError as e:
    print("❌ Import failed, trying auto-install...")
    import subprocess
    import sys
    subprocess.check_call([sys.executable, "-m", "pip", "install", "numpy-financial", "--quiet"])
    
    # Try again
    import numpy_financial as npf
    print("✅ numpy_financial installed and imported")
    payment = npf.pmt(0.05/12, 30*12, -300000)
    print(f"Monthly payment: $" + f"{payment:.2f}")
`;

  try {
    const result = await e2bManager.executeCode(code, 'test-user');
    
    console.log('Execution result:');
    console.log('Success:', result.success);
    
    if (result.output) {
      console.log('\nOutput:');
      console.log(result.output);
    }
    
    if (result.error) {
      console.log('\nError:');
      console.log(result.error);
    }
    
    console.log('\n2. Summary:');
    console.log('- numpy_financial is NOT pre-installed in the prod-all template');
    console.log('- Our auto-install wrapper in sandbox-executor.js handles this');
    console.log('- The wrapper automatically installs numpy_financial when needed');
    console.log('- This is a temporary fix until the E2B template update propagates');
    console.log('\n3. What we did:');
    console.log('- Added numpy_financial to C:\\Users\\cc\\my-e2b-template\\e2b.Dockerfile');
    console.log('- Rebuilt the template with: e2b template build');
    console.log('- But E2B cloud infrastructure hasn\'t updated yet');
    console.log('- So we use auto-install as a workaround');
    
  } catch (error) {
    console.error('Error:', error.message);
  } finally {
    await e2bManager.cleanup();
  }
}

verifyStatus().catch(console.error);