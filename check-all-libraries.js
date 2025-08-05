const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const { Sandbox } = require('@e2b/code-interpreter');

console.log('=== CHECKING ALL LIBRARIES IN PROD-ALL TEMPLATE ===\n');

async function checkAllLibraries() {
  let sandbox;
  try {
    console.log('Creating fresh sandbox with prod-all template...');
    sandbox = await Sandbox.create({
      template: 'prod-all',
      apiKey: process.env.E2B_API_KEY,
      timeout: 60000
    });
    
    console.log(`Sandbox created: ${sandbox.sandboxId}\n`);
    
    // Check Python version and pip
    const versionCode = `
import sys
import subprocess
print(f"Python version: {sys.version}")
print(f"\\nPip version:")
subprocess.run([sys.executable, "-m", "pip", "--version"])
`;
    
    console.log('=== PYTHON AND PIP INFO ===');
    const versionResult = await sandbox.runCode(versionCode, { 
      language: 'python',
      timeoutMs: 10000 
    });
    
    if (versionResult.logs?.stdout) {
      console.log(versionResult.logs.stdout.join(''));
    }
    
    // List all installed packages
    const listPackagesCode = `
import subprocess
import json

# Get list of all installed packages with versions
result = subprocess.run(
    [sys.executable, "-m", "pip", "list", "--format=json"],
    capture_output=True,
    text=True
)

packages = json.loads(result.stdout)
print(f"\\nTotal packages installed: {len(packages)}\\n")

# Group packages by category
math_packages = []
data_packages = []
ml_packages = []
bio_packages = []
viz_packages = []
web_packages = []
finance_packages = []
other_packages = []

# Categorize packages
for pkg in packages:
    name = pkg['name'].lower()
    pkg_str = f"{pkg['name']}=={pkg['version']}"
    
    if any(x in name for x in ['numpy', 'scipy', 'sympy', 'mpmath', 'decimal', 'fractions']):
        math_packages.append(pkg_str)
    elif any(x in name for x in ['pandas', 'statistics', 'statsmodels', 'tables', 'numexpr']):
        data_packages.append(pkg_str)
    elif any(x in name for x in ['scikit', 'sklearn', 'joblib', 'threadpoolctl']):
        ml_packages.append(pkg_str)
    elif any(x in name for x in ['bio', 'splice', 'deep', 'rdkit', 'rdchiral']):
        bio_packages.append(pkg_str)
    elif any(x in name for x in ['matplotlib', 'seaborn', 'plotly', 'pillow', 'kaleido']):
        viz_packages.append(pkg_str)
    elif any(x in name for x in ['requests', 'beautifulsoup', 'urllib', 'certifi', 'charset', 'idna']):
        web_packages.append(pkg_str)
    elif any(x in name for x in ['finance', 'yfinance', 'financial']):
        finance_packages.append(pkg_str)
    else:
        other_packages.append(pkg_str)

# Print categorized packages
print("=== MATH & NUMERICAL PACKAGES ===")
for pkg in sorted(math_packages):
    print(f"  {pkg}")

print("\\n=== DATA ANALYSIS PACKAGES ===")
for pkg in sorted(data_packages):
    print(f"  {pkg}")

print("\\n=== MACHINE LEARNING PACKAGES ===")
for pkg in sorted(ml_packages):
    print(f"  {pkg}")

print("\\n=== BIOLOGY & CHEMISTRY PACKAGES ===")
for pkg in sorted(bio_packages):
    print(f"  {pkg}")

print("\\n=== VISUALIZATION PACKAGES ===")
for pkg in sorted(viz_packages):
    print(f"  {pkg}")

print("\\n=== WEB & NETWORKING PACKAGES ===")
for pkg in sorted(web_packages):
    print(f"  {pkg}")

print("\\n=== FINANCE PACKAGES ===")
for pkg in sorted(finance_packages):
    print(f"  {pkg}")

# Check specifically for numpy_financial
print("\\n=== NUMPY_FINANCIAL CHECK ===")
numpy_financial_installed = any(pkg['name'].lower() == 'numpy-financial' for pkg in packages)
if numpy_financial_installed:
    npf_pkg = next(pkg for pkg in packages if pkg['name'].lower() == 'numpy-financial')
    print(f"✅ numpy-financial is installed: {npf_pkg['name']}=={npf_pkg['version']}")
else:
    print("❌ numpy-financial is NOT installed")

# Try importing it
print("\\nTrying to import numpy_financial...")
try:
    import numpy_financial as npf
    print(f"✅ Import successful! Version: {npf.__version__}")
    print(f"Available functions: {', '.join([f for f in dir(npf) if not f.startswith('_')][:10])}...")
except ImportError as e:
    print(f"❌ Import failed: {e}")
`;
    
    console.log('\n=== INSTALLED PACKAGES BY CATEGORY ===');
    const packagesResult = await sandbox.runCode(listPackagesCode, { 
      language: 'python',
      timeoutMs: 30000 
    });
    
    if (packagesResult.logs?.stdout) {
      console.log(packagesResult.logs.stdout.join(''));
    }
    
    if (packagesResult.logs?.stderr) {
      console.log('\nStderr output:');
      console.log(packagesResult.logs.stderr.join(''));
    }
    
    // Check environment variables
    const envCode = `
import os
print("\\n=== RELEVANT ENVIRONMENT VARIABLES ===")
for key, value in os.environ.items():
    if any(x in key.upper() for x in ['BUILD', 'TAG', 'E2B', 'TEMPLATE', 'VERSION']):
        print(f"{key}: {value}")
`;
    
    const envResult = await sandbox.runCode(envCode, { 
      language: 'python',
      timeoutMs: 5000 
    });
    
    if (envResult.logs?.stdout) {
      console.log(envResult.logs.stdout.join(''));
    }
    
  } catch (error) {
    console.error('Error:', error.message);
  } finally {
    if (sandbox) {
      console.log('\n\nCleaning up sandbox...');
      await sandbox.kill();
    }
  }
}

checkAllLibraries().catch(console.error);