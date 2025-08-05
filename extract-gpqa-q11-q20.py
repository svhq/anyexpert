import pandas as pd
import json
import sys

# Set UTF-8 encoding for console output
sys.stdout.reconfigure(encoding='utf-8')

# Load the GPQA dataset
df = pd.read_csv('GPQA dataset/gpqa_main.csv')

# Extract questions 11-20 (rows 10-19 in 0-indexed pandas)
questions = []

for idx in range(10, 20):
    if idx < len(df):
        row = df.iloc[idx]
        
        # Extract the question and answers
        question_text = row['Question']
        correct_answer = row['Correct Answer']
        incorrect_1 = row['Incorrect Answer 1']
        incorrect_2 = row['Incorrect Answer 2']
        incorrect_3 = row['Incorrect Answer 3']
        explanation = row['Explanation']
        subdomain = row['Subdomain']
        
        # Determine which option (A, B, C, D) is the correct answer
        # This requires careful checking
        options = {
            'A': correct_answer,
            'B': incorrect_1,
            'C': incorrect_2,
            'D': incorrect_3
        }
        
        # For GPQA, the correct answer is always in position A initially
        # But we should verify this
        print(f"\n{'='*80}")
        print(f"Question {idx + 1} - {subdomain}")
        print(f"{'='*80}")
        print(f"\nQuestion Text:\n{question_text}")
        print(f"\nOptions:")
        print(f"A) {options['A']}")
        print(f"B) {options['B']}")
        print(f"C) {options['C']}")
        print(f"D) {options['D']}")
        print(f"\nCorrect Answer: A")
        print(f"\nExplanation preview: {explanation[:200]}...")
        
        # Create structured data
        question_data = {
            "id": f"gpqa-{idx + 1}",
            "category": subdomain,
            "question": question_text,
            "options": options,
            "correct_answer": "A",  # Based on GPQA structure
            "metadata": {
                "domain": row.get('Domain', ''),
                "subdomain": subdomain,
                "explanation": explanation
            }
        }
        
        questions.append(question_data)

# Save to JSON file
output_data = {
    "test_name": "GPQA Questions 11-20",
    "total_questions": len(questions),
    "questions": questions
}

with open('gpqa-q11-q20-extracted.json', 'w', encoding='utf-8') as f:
    json.dump(output_data, f, indent=2, ensure_ascii=False)

print(f"\n{'='*80}")
print(f"Extraction complete! Saved {len(questions)} questions to gpqa-q11-q20-extracted.json")