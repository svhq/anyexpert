// E2B Template Configuration
// This file centralizes our E2B template configuration

const BUILD_TAG = 'npf-20250805-003130'; // Update this when rebuilding template

module.exports = {
  // Template identification
  TEMPLATE_NAME: 'prod-all',
  TEMPLATE_SLUG: 'izdx3u0apwnbdatk6pmh',
  // Note: We don't have the tpl_ ID yet, using name for now
  TEMPLATE_ID: process.env.E2B_TEMPLATE_ID || 'prod-all',
  
  // Build tracking
  EXPECTED_BUILD_TAG: BUILD_TAG,
  
  // Template features
  HAS_NUMPY_FINANCIAL: true,
  
  // Logging
  logTemplateInfo() {
    console.log('=== E2B Template Configuration ===');
    console.log(`Template: ${this.TEMPLATE_ID}`);
    console.log(`Expected BUILD_TAG: ${this.EXPECTED_BUILD_TAG}`);
    console.log(`Has numpy_financial: ${this.HAS_NUMPY_FINANCIAL}`);
    console.log('================================\n');
  }
};