// Retest GLM-4.5-air on the 5 questions it got wrong with clean options
const workflowEngine = require('./src/workflow-engine');
const fs = require('fs');

// The 5 questions GLM-4.5-air got wrong, with cleaned options
const questionsToRetest = [
  {
    id: 'law-q3',
    original_number: 3,
    category: 'law',
    question: "A woman was standing in the aisle of a subway car and put her purse on the seat next to her. A man approached the woman from behind and grabbed the purse off the seat. He then pushed the woman out of the way and ran out of the subway car while carrying the purse. The man was apprehended on the subway platform while in possession of the purse. In a jurisdiction that follows the common law with respect to criminal offenses, of what crime can the man properly be convicted?",
    options: [
      "A) Larceny, because he took the purse without the woman's consent.",
      "B) Burglary, because he entered the subway car with the intention of committing a theft.",
      "C) Robbery, because he used force in leaving with the purse.",
      "D) Robbery, because he used force to take possession of the purse.",
      "E) Theft from a person's vicinity.",
      "F) Assault and larceny.",
      "G) Grand theft.",
      "H) Robbery, because he physically took the purse from the woman's presence."
    ],
    correct_answer: "D",
    glm_previously_selected: "A"
  },
  {
    id: 'law-q4',
    original_number: 4,
    category: 'law',
    question: "A federal grand jury was investigating a corporation whose tanker ship had spilled crude oil into environmentally sensitive waters. The grand jury issued a subpoena requiring the corporation to produce all emails and internal documents regarding the corporation's knowledge of the risks of an oil spill. The corporation has objected, citing its Fifth Amendment privilege against self-incrimination. Can the subpoena be enforced?",
    options: [
      "A) Yes, because the Fifth Amendment privilege only applies to personal testimonies.",
      "B) Yes, because a corporation has no Fifth Amendment privilege.",
      "C) No, because the corporation was not granted use-and-derivative-use immunity.",
      "D) No, because the corporation was not granted transactional immunity.",
      "E) Yes, because the Fifth Amendment privilege does not apply to the compelled production of documents.",
      "F) No, because the documents are protected under attorney-client privilege.",
      "G) No, because a corporation has the same Fifth Amendment rights as an individual.",
      "H) No, because of the corporation's Fourth Amendment rights."
    ],
    correct_answer: "B",
    glm_previously_selected: "H"
  },
  {
    id: 'psych-q6',
    original_number: 6,
    category: 'psychology',
    question: "You receive a phone call from Hermann H., age 28, who says he is 'totally miserable' because of the recent breakup with his girlfriend and that he would like to begin therapy with you. During the first session with Hermann, you find out that his political views are completely repugnant to you, and you feel that you would not enjoy working with him. As an ethical psychologist, you should:",
    options: [
      "A) Discuss your discomfort with Hermann's political views in the first session.",
      "B) Ignore your personal feelings and continue therapy without discussing the difference in political views.",
      "C) Discuss the difference in political views with Hermann only if they become relevant to the psychotherapy process.",
      "D) Suggest Hermann to find a psychologist who shares his political views.",
      "E) Decline Hermann's request for therapy because of the difference in political views.",
      "F) See Hermann in therapy until his current crisis is over and then make a referral if necessary.",
      "G) Tell Hermann outright that his political views are repugnant and continue the therapy.",
      "H) Provide Hermann with appropriate referrals."
    ],
    correct_answer: "H", // Fixed from the corrupted "I"
    glm_previously_selected: "H"
  },
  {
    id: 'bio-q7',
    original_number: 7,
    category: 'biology',
    question: "Which of the following would most likely provide examples of mitotic cell divisions?",
    options: [
      "A) cross section of muscle tissue",
      "B) longitudinal section of a shoot tip",
      "C) longitudinal section of a leaf vein",
      "D) cross section of a fruit",
      "E) cross section of a leaf",
      "F) longitudinal section of a petal",
      "G) longitudinal section of a seed",
      "H) cross section of an anther (site of pollen production in a flower)"
    ],
    correct_answer: "B",
    glm_previously_selected: "H"
  },
  {
    id: 'chem-q10',
    original_number: 10,
    category: 'chemistry',
    question: "Which of the following has an octet of electrons around the central atom?",
    options: [
      "A) BF3",
      "B) BeF2",
      "C) PF5",
      "D) NH4+",
      "E) SF6"
    ],
    correct_answer: "D",
    glm_previously_selected: "E"
  }
];

async function retestQuestion(questionData) {
  console.log('\n' + '█'.repeat(100));
  console.log(`🔍 RETESTING: Question ${questionData.original_number} (${questionData.category.toUpperCase()})`);
  console.log(`📝 Previous GLM answer: ${questionData.glm_previously_selected} | Expected: ${questionData.correct_answer}`);
  console.log('█'.repeat(100));
  
  console.log(`\n📋 QUESTION:`);
  console.log(questionData.question);
  
  console.log(`\n📋 OPTIONS:`);
  questionData.options.forEach(opt => console.log(opt));
  
  const query = `${questionData.question}\n\nOptions:\n${questionData.options.join('\n')}\n\nPlease select the correct answer and explain your reasoning.`;
  
  console.log(`\n🚀 SENDING TO GLM-4.5-AIR...`);
  const startTime = Date.now();
  
  try {
    const response = await workflowEngine.answer(query, {});
    const duration = Date.now() - startTime;
    
    console.log(`\n✅ RESPONSE RECEIVED (${Math.round(duration/1000)}s, ${response.content.length} chars)`);
    
    // Show full response for manual verification
    console.log(`\n📄 FULL MODEL RESPONSE:`);
    console.log('─'.repeat(100));
    console.log(response.content);
    console.log('─'.repeat(100));
    
    // Look for answer patterns
    console.log(`\n🔍 ANSWER EXTRACTION:`);
    
    // Pattern 1: "The answer is X)"
    let match = response.content.match(/(?:answer is|correct answer is|the answer is)[\s:]*([A-H])\)/i);
    if (match) {
      console.log(`   ✅ Found answer via pattern 1: ${match[1]}`);
    }
    
    // Pattern 2: "Option X is correct"
    if (!match) {
      match = response.content.match(/option ([A-H]) is (?:the )?correct/i);
      if (match) {
        console.log(`   ✅ Found answer via pattern 2: ${match[1]}`);
      }
    }
    
    // Pattern 3: Final answer statement
    if (!match) {
      match = response.content.match(/therefore[,\s]*(?:the )?answer is ([A-H])\)/i);
      if (match) {
        console.log(`   ✅ Found answer via pattern 3: ${match[1]}`);
      }
    }
    
    const extractedAnswer = match ? match[1].toUpperCase() : null;
    
    console.log(`\n📊 RESULTS:`);
    console.log(`   - Expected: ${questionData.correct_answer}`);
    console.log(`   - Previously selected: ${questionData.glm_previously_selected}`);
    console.log(`   - Now extracted: ${extractedAnswer || 'NULL'}`);
    console.log(`   - Is correct: ${extractedAnswer === questionData.correct_answer ? '✅ YES' : '❌ NO'}`);
    
    return {
      id: questionData.id,
      category: questionData.category,
      correct_answer: questionData.correct_answer,
      previous_answer: questionData.glm_previously_selected,
      new_extracted_answer: extractedAnswer,
      is_correct: extractedAnswer === questionData.correct_answer,
      response_length: response.content.length,
      duration_ms: duration,
      full_response: response.content,
      timestamp: new Date().toISOString()
    };
    
  } catch (error) {
    console.log(`\n❌ ERROR: ${error.message}`);
    return {
      id: questionData.id,
      error: error.message,
      timestamp: new Date().toISOString()
    };
  }
}

async function runRetest() {
  console.log('🧪 GLM-4.5-AIR RETEST - 5 QUESTIONS WITH CLEAN OPTIONS');
  console.log('Model:', process.env.OPENROUTER_MODEL || 'z-ai/glm-4.5-air:free');
  console.log('Date:', new Date().toISOString());
  console.log(`\n📊 Total questions to retest: ${questionsToRetest.length}`);
  
  const results = [];
  
  for (let i = 0; i < questionsToRetest.length; i++) {
    const question = questionsToRetest[i];
    
    console.log(`\n\n>>> RETESTING ${i + 1}/${questionsToRetest.length} <<<`);
    
    const result = await retestQuestion(question);
    results.push(result);
    
    // Wait between questions
    if (i < questionsToRetest.length - 1) {
      console.log(`\n⏳ Waiting 5 seconds...`);
      await new Promise(resolve => setTimeout(resolve, 5000));
    }
  }
  
  // Summary
  console.log('\n\n' + '█'.repeat(100));
  console.log('📊 RETEST SUMMARY');
  console.log('█'.repeat(100));
  
  const successful = results.filter(r => !r.error);
  const correct = results.filter(r => r.is_correct);
  const correctedFromBefore = results.filter(r => r.is_correct && r.previous_answer !== r.correct_answer);
  
  console.log(`\n📈 RESULTS:`);
  console.log(`- Total questions: ${results.length}`);
  console.log(`- Successful responses: ${successful.length}`);
  console.log(`- Correct answers: ${correct.length}`);
  console.log(`- Errors: ${results.filter(r => r.error).length}`);
  console.log(`- Improved from before: ${correctedFromBefore.length}`);
  
  console.log(`\n📋 DETAILED BREAKDOWN:`);
  results.forEach(result => {
    if (result.error) {
      console.log(`❌ ${result.id}: ERROR`);
    } else {
      const status = result.is_correct ? '✅' : '❌';
      const change = result.new_extracted_answer !== result.previous_answer ? 
        ` (was ${result.previous_answer})` : ' (same as before)';
      console.log(`${status} ${result.id}: ${result.new_extracted_answer} → ${result.correct_answer}${change}`);
    }
  });
  
  // Special note about psych question
  const psychQuestion = results.find(r => r.id === 'psych-q6');
  if (psychQuestion && psychQuestion.is_correct) {
    console.log(`\n✅ NOTE: Psychology question was actually CORRECT before (H was the right answer, not I)`);
  }
  
  // Save results
  const summary = {
    model: 'z-ai/glm-4.5-air:free',
    test_type: 'retest_5_wrong_questions_clean',
    completed_timestamp: new Date().toISOString(),
    summary: {
      total_questions: results.length,
      successful_responses: successful.length,
      correct_answers: correct.length,
      improved_from_before: correctedFromBefore.length
    },
    results: results
  };
  
  fs.writeFileSync('glm45air-retest-results.json', JSON.stringify(summary, null, 2));
  console.log(`\n💾 Results saved to: glm45air-retest-results.json`);
  
  return summary;
}

runRetest().catch(console.error);