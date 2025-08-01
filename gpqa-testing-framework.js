const workflowEngine = require('./src/workflow-engine');
const gpqaData = require('./gpqa-10-questions-clean.json');
const fs = require('fs');
const path = require('path');

class GPQATestingFramework {
  constructor() {
    this.sessionId = `gpqa-session-${new Date().toISOString().replace(/[:.]/g, '-')}`;
    this.logDir = path.join(__dirname, 'gpqa-test-logs');
    this.sessionFile = path.join(this.logDir, `${this.sessionId}.json`);
    
    // Ensure log directory exists
    if (!fs.existsSync(this.logDir)) {
      fs.mkdirSync(this.logDir, { recursive: true });
    }
    
    this.results = {
      sessionId: this.sessionId,
      timestamp: new Date().toISOString(),
      questions: [],
      summary: {}
    };
  }

  createSimplifiedQuery(questionData) {
    return `${questionData.question}

(A) ${questionData.options.A}
(B) ${questionData.options.B}
(C) ${questionData.options.C}
(D) ${questionData.options.D}`;
  }

  async testQuestion(questionIndex) {
    const questionData = gpqaData.questions[questionIndex];
    const questionNum = questionIndex + 1;
    
    console.log(`\n${'='.repeat(70)}`);
    console.log(`🧪 TESTING Q${questionNum}: ${questionData.category}`);
    console.log(`Expected Answer: ${questionData.correct_answer}`);
    console.log(`${'='.repeat(70)}`);
    
    const simplifiedQuery = this.createSimplifiedQuery(questionData);
    
    console.log(`📤 Query being sent:`);
    console.log(`"${simplifiedQuery}"`);
    
    const startTime = Date.now();
    
    try {
      console.log(`\n⏳ Sending to model (no constraints)...`);
      
      const response = await workflowEngine.answer(simplifiedQuery);
      const duration = Date.now() - startTime;
      
      console.log(`✅ Response received (${duration}ms)`);
      console.log(`📝 Response length: ${response.content.length} chars`);
      console.log(`🔧 Tools used: ${response.executionResults ? 'Yes' : 'No'}`);
      
      // Store complete unconstrained response
      const questionResult = {
        questionId: questionData.id,
        questionNumber: questionNum,
        category: questionData.category,
        expectedAnswer: questionData.correct_answer,
        query: simplifiedQuery,
        response: {
          content: response.content,
          executionResults: response.executionResults || null,
          mathExecuted: response.mathExecuted || false,
          metadata: response.metadata || {}
        },
        timing: {
          responseTime: duration,
          timestamp: new Date().toISOString()
        },
        // These will be filled manually
        manualReview: {
          extractedAnswer: null,
          isCorrect: null,
          confidence: null,
          reasoning: null,
          reviewedBy: null,
          reviewTimestamp: null
        }
      };
      
      this.results.questions.push(questionResult);
      
      // Save after each question
      this.saveResults();
      
      console.log(`💾 Full response logged to: ${this.sessionFile}`);
      console.log(`\n📋 MANUAL REVIEW NEEDED:`);
      console.log(`Please read the full response and determine:`);
      console.log(`- What is the final answer (A, B, C, or D)?`);
      console.log(`- Is it correct?`);
      console.log(`- How confident are you in this assessment?`);
      
      return questionResult;
      
    } catch (error) {
      console.error(`❌ Error: ${error.message}`);
      
      const errorResult = {
        questionId: questionData.id,
        questionNumber: questionNum,
        category: questionData.category,
        expectedAnswer: questionData.correct_answer,
        query: simplifiedQuery,
        error: error.message,
        timing: {
          responseTime: Date.now() - startTime,
          timestamp: new Date().toISOString()
        }
      };
      
      this.results.questions.push(errorResult);
      this.saveResults();
      
      return errorResult;
    }
  }

  saveResults() {
    fs.writeFileSync(this.sessionFile, JSON.stringify(this.results, null, 2));
  }

  displayResponseForManualReview(questionResult) {
    console.log(`\n${'='.repeat(80)}`);
    console.log(`📖 MANUAL REVIEW: Q${questionResult.questionNumber}`);
    console.log(`Expected: ${questionResult.expectedAnswer}`);
    console.log(`${'='.repeat(80)}`);
    console.log(`\n📝 FULL UNCONSTRAINED RESPONSE:`);
    console.log(questionResult.response.content);
    console.log(`\n${'='.repeat(80)}`);
    console.log(`❓ MANUAL REVIEW QUESTIONS:`);
    console.log(`1. What is the model's final answer? (A, B, C, or D)`);
    console.log(`2. Is this answer correct? (Yes/No)`);
    console.log(`3. How confident are you? (High/Medium/Low)`);
    console.log(`4. Brief reasoning for your assessment:`);
    console.log(`${'='.repeat(80)}`);
  }

  async runFullTest(startQuestion = 0, endQuestion = 5) {
    console.log(`🎯 GPQA TESTING FRAMEWORK - UNCONSTRAINED RESPONSES`);
    console.log(`Session ID: ${this.sessionId}`);
    console.log(`Testing questions ${startQuestion + 1} to ${endQuestion}`);
    console.log(`Logs will be saved to: ${this.sessionFile}\n`);
    
    const results = [];
    
    for (let i = startQuestion; i < endQuestion && i < gpqaData.questions.length; i++) {
      const result = await this.testQuestion(i);
      results.push(result);
      
      // Brief pause between questions
      if (i < endQuestion - 1) {
        console.log(`\n⏸️ Brief pause before next question...`);
        await new Promise(resolve => setTimeout(resolve, 3000));
      }
    }
    
    console.log(`\n${'='.repeat(70)}`);
    console.log(`🏁 TESTING COMPLETE`);
    console.log(`${'='.repeat(70)}`);
    console.log(`Questions tested: ${results.length}`);
    console.log(`Results saved to: ${this.sessionFile}`);
    console.log(`\n📋 NEXT STEPS:`);
    console.log(`1. Review each response manually`);
    console.log(`2. Use the manual review section to record your findings`);
    console.log(`3. Run the analysis script to get final scores`);
    
    return results;
  }

  static async reviewSession(sessionFile) {
    console.log(`📖 MANUAL REVIEW MODE`);
    console.log(`Session file: ${sessionFile}\n`);
    
    const data = JSON.parse(fs.readFileSync(sessionFile, 'utf8'));
    
    for (const questionResult of data.questions) {
      if (!questionResult.error) {
        console.log(`\n${'='.repeat(80)}`);
        console.log(`📖 QUESTION ${questionResult.questionNumber}: ${questionResult.category}`);
        console.log(`Expected Answer: ${questionResult.expectedAnswer}`);
        console.log(`${'='.repeat(80)}`);
        console.log(`\n📝 FULL RESPONSE:`);
        console.log(questionResult.response.content);
        console.log(`\n${'='.repeat(80)}`);
        console.log(`❓ Please manually determine:`);
        console.log(`- Final answer: A, B, C, or D?`);
        console.log(`- Correct: Yes/No?`);
        console.log(`- Confidence: High/Medium/Low?`);
        console.log(`${'='.repeat(80)}`);
        
        // Wait for user input (in a real scenario, you'd implement interactive input)
        console.log(`⏸️ Review complete for Q${questionResult.questionNumber}\n`);
      }
    }
  }
}

// Example usage
if (require.main === module) {
  const framework = new GPQATestingFramework();
  
  // Test first 5 questions
  framework.runFullTest(0, 5).then(() => {
    console.log('\n🎉 Testing session complete!');
    console.log(`📁 Review the logs at: ${framework.sessionFile}`);
  }).catch(error => {
    console.error('💥 Testing failed:', error);
  });
}

module.exports = GPQATestingFramework;