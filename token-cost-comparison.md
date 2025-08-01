# Token Cost Comparison: Gemini 2.5 Pro vs o4-mini-high

## Pricing (per million tokens)

### Google Gemini 2.5 Pro
- **Input**: $1.25 per million tokens
- **Output**: $10.00 per million tokens

### OpenAI o4-mini-high
- **Input**: $1.10 per million tokens  
- **Output**: $4.40 per million tokens

## Cost Analysis for MMLU Testing

Based on typical MMLU question characteristics:
- Average question + options: ~300-500 input tokens
- Average model response: ~800-1500 output tokens

### Per Question Cost Estimates

**Gemini 2.5 Pro:**
- Input: 400 tokens × $1.25/1M = $0.0005
- Output: 1200 tokens × $10.00/1M = $0.012
- **Total per question: ~$0.0125**

**o4-mini-high:**
- Input: 400 tokens × $1.10/1M = $0.00044
- Output: 1200 tokens × $4.40/1M = $0.00528
- **Total per question: ~$0.00572**

## Cost for 30 Questions

- **Gemini 2.5 Pro**: 30 × $0.0125 = **$0.375**
- **o4-mini-high**: 30 × $0.00572 = **$0.172**

## Conclusion

**o4-mini-high is significantly cheaper** than Gemini 2.5 Pro:
- **54% cheaper overall** for typical MMLU usage
- Input tokens: 12% cheaper
- Output tokens: 56% cheaper

Given that their performance was similar in our tests:
- Gemini 2.5 Pro: 28/30 (93.3%)
- o4-mini-high: 27/30 (90%)

**o4-mini-high offers better value** for MMLU testing, providing nearly identical accuracy at less than half the cost.