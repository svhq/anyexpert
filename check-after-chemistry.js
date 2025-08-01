const fetch = require('node-fetch');

async function checkAfterChemistry() {
  // First run chemistry code to trigger installations
  console.log('Running chemistry code to trigger library installation...\n');
  
  const chemCode = `
# This will trigger library installation
from rdkit import Chem
import numpy as np
import rdchiral

mol = Chem.MolFromSmiles('CCO')
print(f"Ethanol has {mol.GetNumAtoms()} atoms")
print(f"NumPy version: {np.__version__}")
`;

  try {
    let response = await fetch('http://localhost:3001/run_code', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        language: 'python',
        source: chemCode,
        timeout: 60000  // Allow time for installation
      })
    });

    let result = await response.json();
    console.log('Chemistry code result:', result.stdout || 'No output');
    console.log('Exit code:', result.exitCode);
    
    // Now check what's installed
    console.log('\n' + '='.repeat(50));
    console.log('Checking installed libraries after chemistry code...\n');
    
    const checkCode = `
import subprocess
result = subprocess.run([sys.executable, "-m", "pip", "list"], capture_output=True, text=True)
lines = result.stdout.strip().split('\\n')[2:]  # Skip headers

print("Installed packages after chemistry code:")
print("-" * 40)
for line in sorted(lines):
    if line.strip():
        print(line)
`;

    response = await fetch('http://localhost:3001/run_code', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        language: 'python',
        source: checkCode,
        timeout: 10000
      })
    });

    result = await response.json();
    
    if (result.stdout) {
      console.log(result.stdout);
    }
    
  } catch (error) {
    console.error('Error:', error.message);
  }
}

checkAfterChemistry();