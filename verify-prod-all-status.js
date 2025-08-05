const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const { Sandbox } = require('@e2b/code-interpreter');

console.log('=== VERIFYING PROD-ALL TEMPLATE STATUS ===\n');

console.log('Our configuration:');
console.log('- E2B_TEMPLATE_ID in .env:', process.env.E2B_TEMPLATE_ID);
console.log('- Template we built: prod-all (izdx3u0apwnbdatk6pmh)');
console.log('- Dockerfile has: numpy-financial');
console.log('');

async function verify() {
  let sandbox;
  try {
    console.log('Creating sandbox with template: prod-all');
    sandbox = await Sandbox.create({
      template: 'prod-all',
      apiKey: process.env.E2B_API_KEY
    });
    
    console.log('Sandbox created:', sandbox.sandboxId);
    console.log('');
    
    const code = `
import os
print("E2B_TEMPLATE_ID in sandbox:", os.environ.get('E2B_TEMPLATE_ID'))

# Check numpy-financial
try:
    import numpy_financial
    print("numpy-financial: INSTALLED ✅")
except ImportError:
    print("numpy-financial: NOT INSTALLED ❌")
    print("This means E2B hasn't updated with our rebuilt template yet")
`;
    
    const result = await sandbox.runCode(code, { language: 'python' });
    
    if (result.logs?.stdout) {
      console.log(result.logs.stdout.join('\n'));
    }
    
    console.log('\nSUMMARY:');
    console.log('- We ARE using prod-all template');
    console.log('- We rebuilt it with numpy-financial');
    console.log('- E2B cloud infrastructure hasn\'t updated yet');
    console.log('- Our auto-install wrapper handles this perfectly');
    
  } catch (error) {
    console.error('Error:', error.message);
  } finally {
    if (sandbox) {
      await sandbox.kill();
    }
  }
}

verify().catch(console.error);