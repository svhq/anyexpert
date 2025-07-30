const workflowEngine = require('./src/workflow-engine');

console.log('🧮 ASK ANY EXPERT - Math Capabilities Demo\n');

async function demoMath() {
  // Demo 1: Simple arithmetic
  console.log('1️⃣ Simple Arithmetic');
  const q1 = "What is 45 + 78?";
  console.log(`Q: ${q1}`);
  const r1 = await workflowEngine.answer(q1, {});
  console.log(`A: ${r1.content.substring(0, 200)}...`);
  console.log(`Math executed: ${r1.mathExecuted ? 'Yes' : 'No (direct calculation)'}\n`);
  
  // Demo 2: High precision
  console.log('2️⃣ High Precision Calculation');
  const q2 = "Calculate the factorial of 15";
  console.log(`Q: ${q2}`);
  const r2 = await workflowEngine.answer(q2, {});
  console.log(`A: ${r2.content.substring(0, 300)}...`);
  console.log(`Math executed: ${r2.mathExecuted ? 'Yes' : 'No'}`);
  if (r2.executionResults && r2.executionResults[0]?.result?.stdout) {
    console.log(`Computed result: ${r2.executionResults[0].result.stdout}`);
  }
  console.log('');
  
  // Demo 3: Symbolic math
  console.log('3️⃣ Symbolic Mathematics');
  const q3 = "Find the derivative of x^3 - 2x^2 + 5x - 7";
  console.log(`Q: ${q3}`);
  const r3 = await workflowEngine.answer(q3, {});
  console.log(`A: ${r3.content.substring(0, 300)}...`);
  console.log(`Math executed: ${r3.mathExecuted ? 'Yes' : 'No'}\n`);
  
  // Demo 4: Statistics
  console.log('4️⃣ Statistical Analysis');
  const q4 = "Calculate the mean and standard deviation of [12, 15, 18, 20, 22, 25, 28]";
  console.log(`Q: ${q4}`);
  const r4 = await workflowEngine.answer(q4, {});
  console.log(`A: ${r4.content.substring(0, 300)}...`);
  console.log(`Math executed: ${r4.mathExecuted ? 'Yes' : 'No'}\n`);
  
  console.log('✅ Ask Any Expert now handles:');
  console.log('- Simple arithmetic (instant mental math)');
  console.log('- High-precision calculations (factorials, roots)');
  console.log('- Symbolic mathematics (calculus, algebra)');
  console.log('- Statistical analysis');
  console.log('- All with appropriate expert personas!');
}

demoMath().catch(console.error);