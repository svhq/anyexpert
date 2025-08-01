const fetch = require('node-fetch');

async function smokeTest() {
  console.log('Smoke testing RDKit in warm sandbox...\n');
  
  const response = await fetch('http://localhost:3001/run_code', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      language: 'python',
      source: `from rdkit import Chem
mol = Chem.MolFromSmiles('CCO')
print(f"RDKit OK, atoms: {mol.GetNumAtoms()}")`
    })
  });
  
  const result = await response.json();
  console.log('Exit code:', result.exitCode);
  console.log('STDOUT:');
  console.log(result.stdout);
  if (result.stderr) {
    console.log('STDERR:');
    console.log(result.stderr);
  }
}

smokeTest().catch(console.error);