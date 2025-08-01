const fetch = require('node-fetch');

async function demoE2BValue() {
  console.log('E2B VALUE DEMONSTRATION - Real Chemistry Calculations\n');
  console.log('='.repeat(80));
  console.log('Showing calculations that are impossible without computational tools:\n');

  const code = `
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs

# Real-world drug molecules
molecules = {
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Ibuprofen": "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O",
    "Acetaminophen": "CC(=O)Nc1ccc(O)cc1",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Morphine": "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5"
}

print("=== Drug Property Analysis ===\\n")

# Calculate properties for each drug
drug_data = {}
for name, smiles in molecules.items():
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        # Calculate various descriptors
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        tpsa = Descriptors.TPSA(mol)
        rotatable = Descriptors.NumRotatableBonds(mol)
        rings = Descriptors.RingCount(mol)
        
        drug_data[name] = {
            'MW': mw,
            'LogP': logp,
            'HBD': hbd,
            'HBA': hba,
            'TPSA': tpsa,
            'Rotatable': rotatable,
            'Rings': rings
        }
        
        print(f"{name}:")
        print(f"  Molecular Weight: {mw:.1f} g/mol")
        print(f"  LogP: {logp:.2f}")
        print(f"  H-bond donors/acceptors: {hbd}/{hba}")
        print(f"  Polar Surface Area: {tpsa:.1f} Ų")
        print(f"  Rotatable bonds: {rotatable}")
        print()

# Similarity analysis
print("\\n=== Drug Similarity Matrix ===")
print("(Using Tanimoto coefficient with molecular fingerprints)\\n")

# Calculate fingerprints
fps = {}
for name, smiles in molecules.items():
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        fps[name] = FingerprintMols.FingerprintMol(mol)

# Print similarity matrix
names = list(fps.keys())
print("           ", end="")
for name in names:
    print(f"{name[:10]:>11}", end="")
print()

for i, name1 in enumerate(names):
    print(f"{name1[:10]:>11}", end="")
    for j, name2 in enumerate(names):
        sim = DataStructs.TanimotoSimilarity(fps[name1], fps[name2])
        print(f"      {sim:.3f}", end="")
    print()

# Lipinski Rule of Five analysis
print("\\n=== Lipinski's Rule of Five Compliance ===")
print("(For oral bioavailability prediction)\\n")

for name, data in drug_data.items():
    violations = 0
    issues = []
    
    if data['MW'] > 500:
        violations += 1
        issues.append(f"MW > 500 ({data['MW']:.1f})")
    if data['LogP'] > 5:
        violations += 1
        issues.append(f"LogP > 5 ({data['LogP']:.2f})")
    if data['HBD'] > 5:
        violations += 1
        issues.append(f"HBD > 5 ({data['HBD']})")
    if data['HBA'] > 10:
        violations += 1
        issues.append(f"HBA > 10 ({data['HBA']})")
    
    status = "PASS" if violations == 0 else f"FAIL ({violations} violations)"
    print(f"{name}: {status}")
    if issues:
        for issue in issues:
            print(f"  - {issue}")

# Most similar drugs
print("\\n=== Key Findings ===")
max_sim = 0
most_similar = ("", "")
for i in range(len(names)):
    for j in range(i+1, len(names)):
        sim = DataStructs.TanimotoSimilarity(fps[names[i]], fps[names[j]])
        if sim > max_sim:
            max_sim = sim
            most_similar = (names[i], names[j])

print(f"Most similar drugs: {most_similar[0]} and {most_similar[1]} (similarity: {max_sim:.3f})")

# BBB permeability prediction (simple model)
print("\\n=== Blood-Brain Barrier Permeability Prediction ===")
print("(Simple model: MW < 450 and LogP 1-4)\\n")

for name, data in drug_data.items():
    bbb_likely = data['MW'] < 450 and 1 < data['LogP'] < 4
    print(f"{name}: {'Likely' if bbb_likely else 'Unlikely'} to cross BBB")

print("\\n" + "="*60)
print("These calculations demonstrate E2B's value:")
print("- Molecular descriptor calculations (MW, LogP, TPSA, etc)")
print("- Fingerprint-based similarity analysis")
print("- Drug-likeness assessment (Lipinski's Rule)")
print("- ADME property predictions")
print("All computed instantly using RDKit's validated algorithms!")
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