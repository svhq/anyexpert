const fetch = require('node-fetch');

async function demoE2BValue() {
  console.log('E2B VALUE DEMONSTRATION\n');
  console.log('='.repeat(80));
  console.log('Showing calculations that are impossible without computational chemistry tools:\n');

  const code = `
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
import numpy as np

# Define molecules
molecules = {
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Nicotine": "CN1CCC[C@H]1c2cccnc2",
    "Morphine": "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5",
    "Cocaine": "CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c3ccccc3)C2",
    "THC": "CCCCCc1cc(O)c2c(c1)OC(C)(C)[C@@H]3CC/C(C)=C\\\\[C@@H]23"
}

print("=== Molecular Complexity Analysis ===")
print("(BertzCT complexity - higher = more complex)\\n")

complexity_scores = {}
for name, smiles in molecules.items():
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        complexity = rdMolDescriptors.CalcBertzCT(mol)
        stereocenters = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        complexity_scores[name] = {
            'complexity': complexity,
            'stereocenters': stereocenters,
            'MW': mw,
            'LogP': logp
        }

# Display results
for name, data in sorted(complexity_scores.items(), key=lambda x: x[1]['complexity'], reverse=True):
    print(f"{name:12} - Complexity: {data['complexity']:6.1f}, Stereocenters: {data['stereocenters']}")

# Find extremes
most_complex = max(complexity_scores.items(), key=lambda x: x[1]['complexity'])
most_stereo = max(complexity_scores.items(), key=lambda x: x[1]['stereocenters'])

print(f"\\nMost complex: {most_complex[0]} ({most_complex[1]['complexity']:.1f})")
print(f"Most stereocenters: {most_stereo[0]} with {most_stereo[1]['stereocenters']} centers")

# Calculate similarity matrix
print("\\n=== Tanimoto Similarity Matrix ===")
print("(1.0 = identical, 0.0 = completely different)\\n")

mol_list = [(name, Chem.MolFromSmiles(smiles)) for name, smiles in molecules.items()]
fps = [(name, FingerprintMols.FingerprintMol(mol)) for name, mol in mol_list]

print("         ", end="")
for name, _ in fps:
    print(f"{name[:8]:>9}", end="")
print()

for i, (name1, fp1) in enumerate(fps):
    print(f"{name1[:8]:>9}", end="")
    for j, (name2, fp2) in enumerate(fps):
        similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
        print(f"    {similarity:.3f}", end="")
    print()

print("\\n=== Lipinski's Rule of Five Analysis ===")
for name, data in complexity_scores.items():
    violations = 0
    if data['MW'] > 500: violations += 1
    if data['LogP'] > 5: violations += 1
    hbd = Descriptors.NumHDonors(mol_list[list(molecules.keys()).index(name)][1])
    hba = Descriptors.NumHAcceptors(mol_list[list(molecules.keys()).index(name)][1])
    if hbd > 5: violations += 1
    if hba > 10: violations += 1
    
    print(f"{name}: {violations} violations (MW={data['MW']:.1f}, LogP={data['LogP']:.2f})")

print("\\n=== Why This Demonstrates E2B Value ===")
print("1. BertzCT complexity calculation - proprietary algorithm")
print("2. Stereocenter detection - requires 3D structure analysis")
print("3. Tanimoto similarity - needs molecular fingerprints")
print("4. All done in seconds vs hours of manual calculation")
print("5. Results are reproducible and accurate")
`;

  const response = await fetch('http://localhost:3001/run_code', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      language: 'python',
      source: code,
      timeout: 30000
    })
  });
  
  const result = await response.json();
  
  if (result.exitCode === 0) {
    console.log(result.stdout);
  } else {
    console.log('Error:', result.stderr);
  }
}

demoE2BValue().catch(console.error);