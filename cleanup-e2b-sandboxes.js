const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const { Sandbox } = require('@e2b/code-interpreter');

async function cleanupSandboxes() {
  console.log('=== E2B SANDBOX CLEANUP ===\n');
  console.log('API Key:', process.env.E2B_API_KEY ? '✓ Set' : '✗ Missing');
  
  try {
    console.log('Fetching all sandboxes...\n');
    
    // List all running sandboxes
    const sandboxes = await Sandbox.list();
    console.log(`Found ${sandboxes.length} sandboxes\n`);
    
    if (sandboxes.length === 0) {
      console.log('No sandboxes to clean up!');
      return;
    }
    
    // Display sandbox information
    console.log('Sandbox Details:');
    sandboxes.forEach((sandbox, i) => {
      console.log(`${i + 1}. ID: ${sandbox.sandboxId}`);
      console.log(`   Started: ${new Date(sandbox.startedAt).toLocaleString()}`);
      console.log(`   Metadata: ${JSON.stringify(sandbox.metadata || {})}`);
      console.log();
    });
    
    // Ask for confirmation
    console.log(`\nPreparing to close ${sandboxes.length} sandboxes...`);
    console.log('Starting cleanup in 3 seconds (Ctrl+C to cancel)...\n');
    
    await new Promise(resolve => setTimeout(resolve, 3000));
    
    // Close each sandbox
    let successCount = 0;
    let failCount = 0;
    
    for (const sandboxInfo of sandboxes) {
      process.stdout.write(`Closing sandbox ${sandboxInfo.sandboxId}... `);
      try {
        // Connect to the sandbox and kill it
        const sandbox = await Sandbox.connect(sandboxInfo.sandboxId);
        await sandbox.kill();
        console.log('✓ Success');
        successCount++;
      } catch (error) {
        console.log(`✗ Failed: ${error.message}`);
        failCount++;
      }
    }
    
    console.log('\n=== CLEANUP SUMMARY ===');
    console.log(`✓ Successfully closed: ${successCount}`);
    console.log(`✗ Failed to close: ${failCount}`);
    console.log(`Total processed: ${sandboxes.length}`);
    
    // Verify cleanup
    console.log('\nVerifying cleanup...');
    const remainingSandboxes = await Sandbox.list();
    console.log(`Remaining sandboxes: ${remainingSandboxes.length}`);
    
    if (remainingSandboxes.length === 0) {
      console.log('\n✅ All sandboxes successfully cleaned up!');
    } else {
      console.log('\n⚠️  Some sandboxes still remain. You may need to wait or try again.');
    }
    
  } catch (error) {
    console.error('\n❌ Cleanup error:', error.message);
    if (error.message.includes('401')) {
      console.error('Invalid API key. Please check your E2B_API_KEY.');
    }
  }
}

// Run cleanup
cleanupSandboxes().catch(console.error);