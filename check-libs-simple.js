const fetch = require('node-fetch');

async function checkLibs() {
  const code = `
import subprocess
import sys

# Get pip list
result = subprocess.run([sys.executable, "-m", "pip", "list"], capture_output=True, text=True)
print("Installed packages:")
print(result.stdout)
`;

  try {
    const response = await fetch('http://localhost:3001/run_code', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        language: 'python',
        source: code,
        timeout: 10000
      })
    });

    const result = await response.json();
    
    if (result.stdout) {
      console.log(result.stdout);
    }
    
    if (result.exitCode !== 0) {
      console.log('Exit code:', result.exitCode);
      if (result.stderr) console.log('Error:', result.stderr);
    }
    
  } catch (error) {
    console.error('Error:', error.message);
  }
}

checkLibs();