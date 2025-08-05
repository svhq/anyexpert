import pandas as pd
import json
import random
import sys

# Set UTF-8 encoding
sys.stdout.reconfigure(encoding='utf-8')

# Load the GPQA dataset
df = pd.read_csv('GPQA dataset/gpqa_main.csv')

# For GPQA, we need to randomize the answer positions
# Let's create a proper extraction with randomized options

questions = []

for idx in range(10, 20):
    if idx < len(df):
        row = df.iloc[idx]
        
        # Extract the question and answers
        question_text = row['Question']
        correct_answer_text = row['Correct Answer']
        incorrect_1 = row['Incorrect Answer 1']
        incorrect_2 = row['Incorrect Answer 2']
        incorrect_3 = row['Incorrect Answer 3']
        explanation = row['Explanation']
        subdomain = row['Subdomain']
        
        # Create a list of all answers with their correctness
        all_answers = [
            {'text': correct_answer_text, 'is_correct': True},
            {'text': incorrect_1, 'is_correct': False},
            {'text': incorrect_2, 'is_correct': False},
            {'text': incorrect_3, 'is_correct': False}
        ]
        
        # Shuffle the answers to randomize positions
        random.seed(42 + idx)  # Use consistent seed for reproducibility
        random.shuffle(all_answers)
        
        # Create options dictionary and find correct answer letter
        options = {}
        correct_letter = None
        letters = ['A', 'B', 'C', 'D']
        
        for i, answer in enumerate(all_answers):
            options[letters[i]] = answer['text']
            if answer['is_correct']:
                correct_letter = letters[i]
        
        # Print for manual verification
        print(f"\n{'='*80}")
        print(f"Question {idx + 1} - {subdomain}")
        print(f"{'='*80}")
        print(f"\nQuestion: {question_text[:200]}...")
        print(f"\nOptions:")
        for letter in letters:
            print(f"{letter}) {options[letter][:100]}...")
        print(f"\nCorrect Answer: {correct_letter}")
        print(f"Correct Answer Text: {correct_answer_text[:100]}...")
        
        # Create structured data
        question_data = {
            "id": f"gpqa-{idx + 1}",
            "category": subdomain,
            "question": question_text.strip(),
            "options": options,
            "correct_answer": correct_letter,
            "metadata": {
                "domain": row.get('Domain', ''),
                "subdomain": subdomain,
                "explanation": explanation
            }
        }
        
        questions.append(question_data)

# Save to JSON file
output_data = {
    "test_name": "GPQA Questions 11-20 (Randomized)",
    "created_date": pd.Timestamp.now().isoformat(),
    "total_questions": len(questions),
    "questions": questions
}

with open('gpqa-q11-q20-randomized.json', 'w', encoding='utf-8') as f:
    json.dump(output_data, f, indent=2, ensure_ascii=False)

print(f"\n{'='*80}")
print(f"Extraction complete! Saved {len(questions)} questions to gpqa-q11-q20-randomized.json")
print("Note: Answer positions have been randomized for proper testing.")