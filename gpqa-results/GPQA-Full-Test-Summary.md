# GPQA Full Test Results Summary

## Overall Results
- **Total Questions**: 10
- **Correct Answers**: 2/10 (20%)
- **Average Response Time**: ~21.3 seconds
- **Testing Date**: August 1, 2025

## Detailed Results by Question

### Question 1 - Molecular Biology ✅
- **Topic**: Morpholino therapy for exon skipping
- **Expected**: D (R-loops)
- **Model Answer**: D
- **Status**: CORRECT
- **Response Time**: 3.5s

### Question 2 - Physics (general) ❌
- **Topic**: Quantum energy level resolution
- **Expected**: D (10^-4 eV)
- **Model Answer**: Not found
- **Status**: INCORRECT
- **Response Time**: 5.8s

### Question 3 - Organic Chemistry ❌
- **Topic**: Carbon atom count in multi-step synthesis
- **Expected**: A (11)
- **Model Answer**: C (12)
- **Status**: INCORRECT
- **Response Time**: 21.5s

### Question 4 - Molecular Biology ❌
- **Topic**: Methods for investigating chromatin and gene expression
- **Expected**: B (ChIP-seq, chromosome conformation capture, and qRT-PCR)
- **Model Answer**: Not found
- **Status**: INCORRECT
- **Response Time**: 7.6s

### Question 5 - Quantum Mechanics ❌
- **Topic**: Spin state expectation values
- **Expected**: A (-0.7)
- **Model Answer**: Not found
- **Status**: INCORRECT
- **Response Time**: 7.9s

### Question 6 - Chemistry (general) ❌
- **Topic**: Cassie-Baxter contact angles
- **Expected**: B (124°)
- **Model Answer**: D (129°)
- **Status**: INCORRECT
- **Response Time**: 27.9s

### Question 7 - Genetics ❌
- **Topic**: Three-point testcross mapping
- **Expected**: D (Double crossover event)
- **Model Answer**: A (Recombinant interference)
- **Status**: INCORRECT
- **Response Time**: 7.2s

### Question 8 - Electromagnetism and Photonics ❌
- **Topic**: Maxwell's equations with magnetic monopoles
- **Expected**: C (Circulation of electric field and divergence of magnetic field)
- **Model Answer**: Not found
- **Status**: INCORRECT
- **Response Time**: 5.9s

### Question 9 - Organic Chemistry ❌
- **Topic**: Cycloaddition reactions (Diels-Alder)
- **Expected**: A (5-methylcyclohex-3-ene-1-carbonitrile, endo product)
- **Model Answer**: C (5-methylcyclohex-3-ene-1-carbonitrile, exo product)
- **Status**: INCORRECT (wrong stereochemistry for product B)
- **Response Time**: 122.2s
- **Note**: Very long response (149,314 characters)

### Question 10 - Genetics ✅
- **Topic**: Gene interactions and epistasis
- **Expected**: D (G2 transcription factor, G1 and G3 redundancy, G1 epistatic to G3)
- **Model Answer**: D
- **Status**: CORRECT
- **Response Time**: 22.1s

## Analysis by Category
- **Molecular Biology**: 1/2 (50%)
- **Physics**: 0/1 (0%)
- **Organic Chemistry**: 0/2 (0%)
- **Quantum Mechanics**: 0/1 (0%)
- **Chemistry (general)**: 0/1 (0%)
- **Genetics**: 1/2 (50%)
- **Electromagnetism**: 0/1 (0%)

## Key Observations

1. **Answer Extraction Issues**: 4/10 questions had no clear answer extracted, suggesting the model may not be consistently formatting final answers in a recognizable way.

2. **Chemistry Performance**: Chemistry questions performed poorly (0/3 correct), despite having RDKit libraries available. The model seems to struggle with:
   - Stereochemistry (Q9)
   - Physical chemistry calculations (Q6)
   - Carbon counting in synthesis (Q3)

3. **Response Time Variability**: Response times ranged from 3.5s to 122.2s, with chemistry questions taking significantly longer.

4. **Model Routing**: Questions were routed to different agents:
   - math-complex: Q1, Q3, Q9, Q10
   - unknown: Q2, Q4, Q5, Q6, Q7, Q8

5. **Long Responses**: Q9 generated an extremely long response (149,314 characters), suggesting the model may be getting stuck in loops or over-explaining.

## Recommendations

1. **Answer Format**: Need to enforce clearer answer selection format (e.g., "Final Answer: (X)")
2. **Chemistry Tools**: Despite having chemistry libraries, the model needs better prompting for chemical reasoning
3. **Response Length**: Consider adding constraints to prevent extremely long responses
4. **Agent Routing**: Review why many questions are routing to "unknown" agent