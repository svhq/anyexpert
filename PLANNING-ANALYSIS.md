# Planning Addition Analysis

## Date: 2025-08-03

### Summary of Changes
Added "First, make a plan" instruction to the planning assistant prompt in `unified-agent.js` along with a comprehensive list of available tools and libraries.

### Test Results

#### Overall Performance
- **Questions Tested**: Q11, Q14, Q15, Q16
- **Correct**: 0/4 (0%)
- **Previous Performance**: Q11 was wrong before, Q14 timed out, Q15-Q16 untested

#### Detailed Results

| Question | Category | Previous | With Planning | Tools Used | Duration |
|----------|----------|----------|---------------|------------|----------|
| Q11 | Molecular Biology | ❌ A | ❌ B | reason → synthesize | 16.4s |
| Q14 | Quantum Mechanics | Timeout | ❌ B | reason×6 (maxed out) | 523s |
| Q15 | Astrophysics | Untested | ❌ C | synthesize only | 11.1s |
| Q16 | Particle Physics | Untested | ❌ B | search×2 → reason → synthesize | 24.7s |

### Key Observations

#### 1. Planning Did Not Improve Accuracy
Despite adding planning instructions, all 4 questions were answered incorrectly. The planning prompt did not lead to better tool selection or more accurate answers.

#### 2. Tool Usage Patterns

**Q11 (Molecular Biology)**:
- Used: reason → synthesize
- Should have used: search (for specialized biology methods)
- Issue: Model didn't recognize need for domain-specific knowledge despite mentioning biopython in tools

**Q14 (Quantum Mechanics)**:
- Used: reason×6 (hit max steps)
- Should have used: code (for calculations)
- Issue: Got stuck in reasoning loop without using computational tools

**Q15 (Astrophysics)**:
- Used: synthesize only
- Should have used: search or code
- Issue: Jumped straight to answer without gathering information

**Q16 (Particle Physics)**:
- Used: search×2 → reason → synthesize
- Good: Actually used search for specialized physics knowledge
- Bad: Still got wrong answer (B instead of D)

#### 3. Planning Behavior
Looking at the logs, the model did receive the planning prompt with tool listings, but:
- It didn't create comprehensive upfront plans
- Tool selection remained reactive rather than strategic
- The biological tools (biopython, etc.) were never considered for Q11

#### 4. JSON Parsing Issues
Multiple instances of "Failed to parse JSON response" suggest the model sometimes struggled with the structured format requirement.

### Conclusions

1. **Planning Alone Is Insufficient**: Simply adding "make a plan" didn't improve performance. The model needs stronger guidance on when to use specific tools.

2. **Domain Recognition Problem**: The model fails to recognize when specialized tools (like biopython for molecular biology) would be helpful.

3. **Tool Selection Heuristics Too Weak**: Current heuristics in `fallbackPlan()` are too simple - they only look for keywords like "calculate" or "current".

### Recommendations

1. **Strengthen Domain Detection**: Add explicit domain recognition that maps question types to relevant tools.

2. **Mandatory Tool Exploration**: For complex scientific questions, require at least one search or code execution before synthesizing.

3. **Example-Based Planning**: Add few-shot examples in the planning prompt showing how to match question types to appropriate tools.

4. **Confidence Thresholds**: Lower the initial confidence for pure reasoning to encourage tool use.

5. **Parallel Execution**: The current restriction of parallel execution to first step only is limiting. Consider allowing it at any step.

### Next Steps

The planning addition shows that simply asking the model to plan is not enough. The system needs:
- Better domain recognition
- Stronger tool selection guidance  
- More sophisticated confidence assessment
- Examples of good planning behavior

Without these improvements, the model defaults to reasoning alone, missing opportunities to leverage specialized tools that could provide accurate answers.