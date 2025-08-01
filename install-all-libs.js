const http = require('http');

// Install all required libraries
const installCode = `
import subprocess
import sys

libs_to_install = [
    'numpy', 'pandas', 'sympy', 'scipy', 'matplotlib', 
    'seaborn', 'scikit-learn', 'statsmodels', 'plotly',
    'networkx', 'beautifulsoup4', 'lxml', 'pillow', 'mpmath', 'requests'
]

print("Installing all required libraries...")
for lib in libs_to_install:
    print(f"Installing {lib}...")
    result = subprocess.run([sys.executable, '-m', 'pip', 'install', lib], 
                          capture_output=True, text=True)
    if result.returncode == 0:
        print(f"✅ {lib} installed")
    else:
        print(f"❌ {lib} failed: {result.stderr[:100]}")

print("\\nVerifying imports...")
# Test imports
import numpy as np
import pandas as pd
import sympy as sp
import scipy
print("✅ Core libraries imported successfully")
print(f"numpy: {np.__version__}")
print(f"pandas: {pd.__version__}")
print(f"sympy: {sp.__version__}")
print(f"scipy: {scipy.__version__}")
`;

const data = JSON.stringify({
  language: 'python',
  source: installCode.trim(),
  timeout: 300000 // 5 minutes
});

const options = {
  hostname: 'localhost',
  port: 3001,
  path: '/run_code',
  method: 'POST',
  headers: {
    'Content-Type': 'application/json',
    'Content-Length': Buffer.byteLength(data)
  }
};

console.log('Installing all libraries in E2B (this may take a few minutes)...');

const req = http.request(options, (res) => {
  let body = '';
  res.on('data', (chunk) => body += chunk);
  res.on('end', () => {
    try {
      const result = JSON.parse(body);
      console.log('Exit Code:', result.exitCode);
      console.log('Output:');
      console.log(result.stdout);
      if (result.stderr) {
        console.log('Errors:');
        console.log(result.stderr);
      }
    } catch (e) {
      console.error('Parse error:', e.message);
      console.log('Raw response:', body);
    }
  });
});

req.on('error', (e) => {
  console.error('Request error:', e.message);
});

req.write(data);
req.end();