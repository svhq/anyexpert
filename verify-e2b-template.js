const fetch = require('node-fetch');
const colors = {
  green: '\x1b[32m',
  red: '\x1b[31m',
  yellow: '\x1b[33m',
  blue: '\x1b[34m',
  reset: '\x1b[0m'
};

async function verifyTemplate() {
  console.log(`${colors.blue}Verifying E2B Custom Template Setup...${colors.reset}\n`);
  
  // Check health endpoint
  try {
    const healthRes = await fetch('http://localhost:3001/health');
    const health = await healthRes.json();
    
    console.log('🏥 Health Check:');
    console.log(`   Status: ${health.status === 'healthy' ? colors.green + '✓' : colors.red + '✗'} ${health.status}${colors.reset}`);
    console.log(`   E2B Key: ${health.e2bKey === 'configured' ? colors.green + '✓' : colors.red + '✗'} ${health.e2bKey}${colors.reset}`);
    console.log(`   Template ID: ${health.templateId === 'u8w4jzut60mgigs4ofp9' ? colors.green + '✓' : colors.red + '✗'} ${health.templateId || 'none'}${colors.reset}`);
    
    if (health.templateId !== 'u8w4jzut60mgigs4ofp9') {
      console.log(`\n${colors.yellow}⚠️  Template not loaded! The service needs to be restarted.${colors.reset}`);
      console.log(`${colors.yellow}   Run: cd microservices/run-code && node server-e2b-custom.js${colors.reset}`);
      return;
    }
  } catch (error) {
    console.log(`${colors.red}❌ Cannot connect to E2B service at http://localhost:3001${colors.reset}`);
    return;
  }
  
  // Test template functionality
  console.log('\n📋 Testing Template Features:');
  
  const tests = [
    {
      name: 'Template Info',
      code: `
import os
with open('/home/user/template_info.txt', 'r') as f:
    print(f.read().strip())
`,
      expected: 'Custom E2B Template'
    },
    {
      name: 'NumPy Version',
      code: 'import numpy; print(numpy.__version__)',
      expected: '1.'
    },
    {
      name: 'RDKit Import',
      code: `
from rdkit import Chem
from rdkit.Chem import Descriptors
mol = Chem.MolFromSmiles('CCO')
print(f"Ethanol MW: {Descriptors.MolWt(mol):.2f}")
`,
      expected: '46.07'
    },
    {
      name: 'RDChiral Import',
      code: `
import rdchiral
print("RDChiral imported successfully")
`,
      expected: 'successfully'
    }
  ];
  
  let allPassed = true;
  
  for (const test of tests) {
    process.stdout.write(`   ${test.name}: `);
    
    try {
      const response = await fetch('http://localhost:3001/run_code', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          language: 'python',
          source: test.code,
          timeout: 10000
        })
      });
      
      const result = await response.json();
      
      if (result.exitCode === 0 && result.stdout && result.stdout.includes(test.expected)) {
        console.log(`${colors.green}✓ PASS${colors.reset}`);
      } else {
        console.log(`${colors.red}✗ FAIL${colors.reset}`);
        if (result.stdout) console.log(`     Got: ${result.stdout.trim()}`);
        if (result.stderr) console.log(`     Error: ${result.stderr.split('\\n')[0]}`);
        allPassed = false;
      }
    } catch (error) {
      console.log(`${colors.red}✗ ERROR: ${error.message}${colors.reset}`);
      allPassed = false;
    }
  }
  
  console.log('\n' + '─'.repeat(50));
  if (allPassed) {
    console.log(`${colors.green}✅ All tests passed! Custom template is working correctly.${colors.reset}`);
  } else {
    console.log(`${colors.red}❌ Some tests failed. Check the configuration.${colors.reset}`);
  }
}

verifyTemplate();