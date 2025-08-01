const fetch = require('node-fetch');

async function debug() {
  const response = await fetch('http://localhost:3001/run_code', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      language: 'python',
      source: `
try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    mol = Chem.MolFromSmiles("CC")
    complexity = rdMolDescriptors.CalcBertzCT(mol)
    print("Success! Complexity:", complexity)
except Exception as e:
    print("Error:", str(e))
    import traceback
    traceback.print_exc()
`
    })
  });
  
  const result = await response.json();
  console.log('Result:', result);
}

debug();