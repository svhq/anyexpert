async function debugServiceCall() {
  console.log('Debug: Testing service with simple code first...');
  
  // Test 1: Simple Python without numpy
  const simpleCode = `
print("Hello from Python!")
import sys
print(f"Python version: {sys.version}")
`;

  try {
    console.log('Test 1: Simple Python code');
    const response1 = await fetch('http://localhost:3001/run_code', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        language: 'python',
        source: simpleCode
      })
    });

    const result1 = await response1.json();
    console.log('Simple code result:');
    console.log('Exit code:', result1.exitCode);
    console.log('Stdout:', result1.stdout);
    console.log('Stderr:', result1.stderr);
    
    if (result1.exitCode === 0) {
      console.log('✅ Basic Python execution works');
      
      // Test 2: With numpy
      console.log('\nTest 2: Code with numpy');
      const numpyCode = `
import numpy as np
print(f"NumPy version: {np.__version__}")
print(f"Pi = {np.pi}")
`;
      
      const response2 = await fetch('http://localhost:3001/run_code', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          language: 'python',
          source: numpyCode
        })
      });

      const result2 = await response2.json();
      console.log('NumPy code result:');
      console.log('Exit code:', result2.exitCode);
      console.log('Stdout:', result2.stdout);
      console.log('Stderr:', result2.stderr);
      
    } else {
      console.log('❌ Basic Python execution failed');
    }
    
  } catch (error) {
    console.error('Error:', error.message);
  }
}

debugServiceCall();