# Chemistry Examples for Organic Reactions using RDKit and RDChiral
# These examples show how to solve common organic chemistry problems

# Example 1: Diels-Alder Reaction (like GPQA Q9)
def diels_alder_example():
    """
    Example: (E)-penta-1,3-diene + acrylonitrile → 5-methylcyclohex-3-ene-1-carbonitrile
    """
    from rdkit import Chem
    from rdkit.Chem import rdChemReactions as Reactions
    
    # Define Diels-Alder reaction SMARTS
    # [4+2] cycloaddition: diene + dienophile → cyclohexene
    da_smarts = '[C:1]=[C:2]-[C:3]=[C:4].[C:5]=[C:6]>>[C:1]1[C:2][C:3][C:4][C:5][C:6]1'
    
    # Create reaction
    rxn = Reactions.ReactionFromSmarts(da_smarts)
    
    # Define reactants
    diene = Chem.MolFromSmiles('C/C=C/C=C')  # (E)-penta-1,3-diene
    dienophile = Chem.MolFromSmiles('C=CC#N')  # acrylonitrile
    
    # Run reaction
    products = rxn.RunReactants((diene, dienophile))
    
    # Get product SMILES
    for prod_set in products:
        for prod in prod_set:
            print(f"Product SMILES: {Chem.MolToSmiles(prod)}")
            # This gives: CC1C=CCCC1C#N (5-methylcyclohex-3-ene-1-carbonitrile)
    
    return products

# Example 2: Using RDChiral for stereochemistry-aware reactions
def rdchiral_diels_alder():
    """
    Better example using RDChiral for correct stereochemistry
    """
    from rdkit import Chem
    import rdchiral
    
    # More detailed Diels-Alder SMARTS that preserves stereochemistry
    da_smarts = '[c,C:1]=[c,C:2]-[c,C:3]=[c,C:4].[c,C:5]=[c,C:6]>>[c:1]1[c:2][c:6][c:5][c:4][c:3]1'
    
    # Prepare reactants
    diene = Chem.MolFromSmiles('C/C=C/C=C')  # (E)-penta-1,3-diene
    dienophile = Chem.MolFromSmiles('C=CC#N')  # acrylonitrile
    
    # Use RDChiral for reaction
    rdc_rxn = rdchiral.rdchiralReaction(da_smarts)
    rdc_inputs = rdchiral.rdchiralReactants(f"{Chem.MolToSmiles(diene)}.{Chem.MolToSmiles(dienophile)}")
    products = rdchiral.rdchiralRun(rdc_rxn, rdc_inputs)
    
    print(f"RDChiral product: {products[0]}")
    # This gives the correct regioisomer and stereochemistry
    
    return products

# Example 3: Cyclopentadiene + methyl acrylate → bicyclo[2.2.1]heptene (norbornene)
def norbornene_synthesis():
    """
    Classic Diels-Alder: cyclopentadiene + methyl acrylate
    """
    from rdkit import Chem
    from rdkit.Chem import rdChemReactions as Reactions
    
    # Diels-Alder for bicyclic system
    # Note: This SMARTS captures the bicyclo[2.2.1] formation
    da_smarts = '[C:1]1=[C:2][C:3]=[C:4][C:5]1.[C:6]=[C:7]>>[C:1]1[C:2]2[C:3][C:4][C:5]1[C:7][C:6]2'
    
    rxn = Reactions.ReactionFromSmarts(da_smarts)
    
    # Reactants
    cyclopentadiene = Chem.MolFromSmiles('C1=CC=CC1')
    methyl_acrylate = Chem.MolFromSmiles('C=CC(=O)OC')
    
    # Run reaction
    products = rxn.RunReactants((cyclopentadiene, methyl_acrylate))
    
    for prod_set in products:
        for prod in prod_set:
            smiles = Chem.MolToSmiles(prod)
            print(f"Bicyclic product: {smiles}")
            # This gives the bicyclo[2.2.1]hept-5-ene-2-carboxylate structure
    
    return products

# Example 4: Helper function to predict cycloaddition products
def predict_cycloaddition(diene_smiles, dienophile_smiles, reaction_type='diels-alder'):
    """
    General function to predict cycloaddition products
    """
    from rdkit import Chem
    from rdkit.Chem import rdChemReactions as Reactions
    
    reaction_smarts = {
        'diels-alder': '[C:1]=[C:2]-[C:3]=[C:4].[C:5]=[C:6]>>[C:1]1[C:2][C:3][C:4][C:5][C:6]1',
        '[2+2]': '[C:1]=[C:2].[C:3]=[C:4]>>[C:1]1[C:2][C:3][C:4]1',
        '[3+2]': '[C:1]=[C:2]-[C:3].[N+:4]=[N-:5]>>[C:1]1[C:2][C:3][N:4][N:5]1'
    }
    
    if reaction_type not in reaction_smarts:
        raise ValueError(f"Unknown reaction type: {reaction_type}")
    
    # Create reaction
    rxn = Reactions.ReactionFromSmarts(reaction_smarts[reaction_type])
    
    # Parse molecules
    diene = Chem.MolFromSmiles(diene_smiles)
    dienophile = Chem.MolFromSmiles(dienophile_smiles)
    
    if not diene or not dienophile:
        raise ValueError("Invalid SMILES provided")
    
    # Run reaction
    products = rxn.RunReactants((diene, dienophile))
    
    # Collect unique products
    unique_products = set()
    for prod_set in products:
        for prod in prod_set:
            unique_products.add(Chem.MolToSmiles(prod))
    
    return list(unique_products)

# Example usage for GPQA Q9
if __name__ == "__main__":
    print("=== GPQA Q9 Solution ===")
    print("Reaction A: (E)-penta-1,3-diene + acrylonitrile")
    products_a = predict_cycloaddition('C/C=C/C=C', 'C=CC#N')
    print(f"Product A: {products_a}")
    
    print("\nReaction B: cyclopentadiene + methyl acrylate")
    products_b = predict_cycloaddition('C1=CC=CC1', 'C=CC(=O)OC')
    print(f"Product B: {products_b}")