const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const { Sandbox } = require('@e2b/code-interpreter');

console.log('=== INVESTIGATING TEMPLATE ID MISMATCH ===\n');

console.log('Expected template: prod-all (izdx3u0apwnbdatk6pmh)');
console.log('But sandboxes show: nlhz8vlwyupq845jsdg9\n');

async function investigate() {
  console.log('1. Checking our configuration:');
  console.log(`   - .env E2B_TEMPLATE_ID: ${process.env.E2B_TEMPLATE_ID || 'not set'}`);
  console.log(`   - Using template: prod-all`);
  
  console.log('\n2. Testing sandbox creation:');
  
  let sandbox;
  try {
    // Create with explicit prod-all
    sandbox = await Sandbox.create({
      template: 'prod-all',
      apiKey: process.env.E2B_API_KEY,
      timeout: 30000
    });
    
    console.log(`   - Sandbox created: ${sandbox.sandboxId}`);
    console.log(`   - Sandbox object template: ${sandbox.template || 'not available'}`);
    
    // Check environment
    const checkCode = `
import os
print(f"E2B_TEMPLATE_ID in sandbox: {os.environ.get('E2B_TEMPLATE_ID')}")
print(f"E2B_SANDBOX_ID: {os.environ.get('E2B_SANDBOX_ID')}")

# Check if this might be a base template
print("\\nChecking base image info:")
try:
    with open('/etc/os-release', 'r') as f:
        for line in f:
            if 'PRETTY_NAME' in line:
                print(f"OS: {line.strip()}")
except:
    pass

# Check Docker info if available
import subprocess
try:
    result = subprocess.run(['cat', '/proc/version'], capture_output=True, text=True)
    if result.returncode == 0:
        print(f"Kernel: {result.stdout.strip()}")
except:
    pass
`;
    
    const result = await sandbox.runCode(checkCode, {
      language: 'python',
      timeoutMs: 10000
    });
    
    if (result.logs?.stdout) {
      console.log('\n3. Sandbox environment info:');
      result.logs.stdout.forEach(line => console.log(`   ${line}`));
    }
    
  } catch (error) {
    console.error('Error:', error.message);
  } finally {
    if (sandbox) {
      await sandbox.kill();
    }
  }
  
  console.log('\n4. Possible explanations:');
  console.log('   - nlhz8vlwyupq845jsdg9 might be E2B\'s base code-interpreter template');
  console.log('   - Our prod-all template might be layered on top of this base');
  console.log('   - E2B might be using internal template IDs vs user-facing names');
  console.log('   - The template system might have changed and our custom template isn\'t being used');
}

investigate().catch(console.error);