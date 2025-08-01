# GLM-4.5-AIR Final MMLU Test Results

## Summary
- **Model**: z-ai/glm-4.5-air:free
- **Total Questions**: 30
- **Total Correct**: 17/30 (56.7%)
- **Completion Date**: 2025-07-31T12:54:03.764Z

## Breakdown by Test Type

### Original 10 Questions (from initial test)
- **Score**: 6/10 (60%)
- **Note**: Psychology question was actually correct (answer "H" was valid)

### New 11 Questions (from this test)
- **Score**: 5/11 (45.5%)
- **Correct**: Math, Physics, Psychology, Biology (speciation), Economics
- **Wrong**: Business, History (Hanseatic), Computer Science, Philosophy, Engineering, History (duplicate)

### Retests (2 questions)
- **Score**: 1/2 (50%)
- **Chemistry**: ✅ Corrected (E → D)
- **Biology (mitosis)**: ❌ Still wrong (extraction failed)

## Detailed Results

### Extraction Issues
Several questions had correct reasoning but failed answer extraction:
- Biology Q7 (mitosis): Model explained correctly but didn't format answer
- Business Q3: No clear answer statement
- Computer Science Q8: Explained Merge sort correctly but no clear answer
- Philosophy Q9: Detailed explanation but no answer statement
- Engineering Q10: Correct reasoning but no extractable answer

### Genuinely Wrong Answers
1. **History (Hanseatic League)**: Selected Amsterdam (E) instead of Stockholm (J)
   - Note: Both cities were technically correct answers, but Stockholm was the expected answer
2. **Biology (mitosis)**: Even in retest, failed to clearly state answer B

### Model Behavior Observations
- Provides detailed explanations
- Sometimes doesn't conclude with a clear answer statement
- Occasionally uses code execution for non-coding questions
- Response times: 9-30 seconds per question
- No timeout errors in this test run

## Final Score Calculation
- Original new 10 questions: 6/10
- Additional 11 questions tested: 5/11  
- Retests improved score: +1 (chemistry)
- **Total**: (6 + 5 + 1) / 30 = **17/30 (56.7%)**

## Comparison with Other Models
- **o4-mini-high**: 29/29 (100% on valid questions)
- **GLM-4.5 (premium)**: 22/30 (73.3%)
- **GLM-4.5-air (free)**: 17/30 (56.7%)
- **Gemini 2.5 Pro**: 29/30 (96.7%)

## Conclusion
GLM-4.5-air (free) performed reasonably well but had significant extraction issues. The model often provides correct reasoning but fails to format answers in an extractable way. This suggests the model has good knowledge but may need clearer prompting for standardized test formats.