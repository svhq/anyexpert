const fetch = require('node-fetch');

async function verifyAspirin() {
  console.log('Verifying Aspirin molecular properties...\n');
  
  const code = `
from rdkit import Chem
from rdkit.Chem import Descriptors

# Aspirin SMILES
smiles = "CC(=O)Oc1ccccc1C(=O)O"
mol = Chem.MolFromSmiles(smiles)

# Get molecular formula
formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
print(f"Molecular Formula: {formula}")

# Count atoms
print(f"Carbon atoms: {sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')}")
print(f"Hydrogen atoms: {sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'H')}")
print(f"Oxygen atoms: {sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O')}")

# Calculate properties
mw = Descriptors.MolWt(mol)
logp = Descriptors.MolLogP(mol)

print(f"\\nMolecular Weight: {mw:.3f} g/mol")
print(f"LogP: {logp:.3f}")

# Additional properties
print(f"\\nRotatable bonds: {Descriptors.NumRotatableBonds(mol)}")
print(f"H-bond donors: {Descriptors.NumHDonors(mol)}")
print(f"H-bond acceptors: {Descriptors.NumHAcceptors(mol)}")
print(f"TPSA: {Descriptors.TPSA(mol):.2f}")
`;

  const response = await fetch('http://localhost:3001/run_code', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      language: 'python',
      source: code
    })
  });
  
  const result = await response.json();
  console.log(result.stdout);
}

verifyAspirin().catch(console.error);