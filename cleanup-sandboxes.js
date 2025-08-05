const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const { Sandbox } = require('@e2b/code-interpreter');

async function cleanupAllSandboxes() {
  console.log('=== E2B SANDBOX CLEANUP ===\n');
  
  try {
    // List all running sandboxes
    console.log('Fetching list of running sandboxes...');
    const sandboxes = await Sandbox.list();
    
    console.log(`Found ${sandboxes.length} running sandboxes\n`);
    
    if (sandboxes.length === 0) {
      console.log('No sandboxes to clean up!');
      return;
    }
    
    // Kill each sandbox
    let killed = 0;
    let failed = 0;
    
    for (const sandboxInfo of sandboxes) {
      try {
        console.log(`Killing sandbox ${sandboxInfo.sandboxId}...`);
        
        // Connect to the sandbox and kill it
        const sandbox = await Sandbox.connect(sandboxInfo.sandboxId);
        await sandbox.kill();
        
        killed++;
        console.log(`✅ Successfully killed sandbox ${sandboxInfo.sandboxId}`);
      } catch (error) {
        failed++;
        console.log(`❌ Failed to kill sandbox ${sandboxInfo.sandboxId}: ${error.message}`);
      }
    }
    
    console.log('\n=== CLEANUP SUMMARY ===');
    console.log(`Total sandboxes: ${sandboxes.length}`);
    console.log(`Successfully killed: ${killed}`);
    console.log(`Failed to kill: ${failed}`);
    
    // Verify cleanup
    console.log('\nVerifying cleanup...');
    const remainingSandboxes = await Sandbox.list();
    console.log(`Remaining sandboxes: ${remainingSandboxes.length}`);
    
    if (remainingSandboxes.length > 0) {
      console.log('\n⚠️  Some sandboxes are still running:');
      remainingSandboxes.forEach(s => {
        console.log(`  - ${s.sandboxId}`);
      });
    } else {
      console.log('\n✅ All sandboxes successfully cleaned up!');
    }
    
  } catch (error) {
    console.error('Error during cleanup:', error);
  }
}

cleanupAllSandboxes().catch(console.error);