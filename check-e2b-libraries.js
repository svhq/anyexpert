const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const e2bManager = require('./src/e2b-manager');

console.log('=== CHECKING E2B INSTALLED LIBRARIES ===\n');
console.log('E2B Template ID:', process.env.E2B_TEMPLATE_ID);
console.log('Date:', new Date().toISOString());

// Libraries to check
const librariesToCheck = {
  'Math Core': ['numpy', 'sympy', 'scipy', 'mpmath'],
  'Bio Helpers': ['splicekit', 'rnalib', 'deeptools', 'pandas', 'biopython'],
  'Chemistry': ['rdkit', 'chempy'],
  'General Scientific': ['matplotlib', 'seaborn', 'scikit-learn'],
  'Additional': ['requests', 'beautifulsoup4', 'lxml']
};

async function checkLibraries() {
  console.log('\nChecking installed Python libraries in E2B environment...\n');
  
  // Python code to check all libraries
  const checkCode = `
import subprocess
import json

# Get all installed packages
result = subprocess.run(['pip', 'list', '--format=json'], capture_output=True, text=True)
installed_packages = json.loads(result.stdout)
installed_dict = {pkg['name'].lower(): pkg['version'] for pkg in installed_packages}

# Libraries to check
libraries_to_check = ${JSON.stringify(librariesToCheck)}

# Check each category
results = {}
for category, libs in libraries_to_check.items():
    results[category] = {}
    for lib in libs:
        lib_lower = lib.lower()
        if lib_lower in installed_dict:
            results[category][lib] = {
                'installed': True,
                'version': installed_dict[lib_lower]
            }
        else:
            results[category][lib] = {
                'installed': False,
                'version': None
            }

# Print summary
print("\\n=== INSTALLED LIBRARIES ===")
for category, libs in results.items():
    print(f"\\n{category}:")
    for lib, info in libs.items():
        if info['installed']:
            print(f"  ✓ {lib} v{info['version']}")
        else:
            print(f"  ✗ {lib} (NOT INSTALLED)")

# Count totals
total_checked = sum(len(libs) for libs in libraries_to_check.values())
total_installed = sum(1 for cat in results.values() for info in cat.values() if info['installed'])

print(f"\\n=== SUMMARY ===")
print(f"Total libraries checked: {total_checked}")
print(f"Installed: {total_installed}")
print(f"Missing: {total_checked - total_installed}")

# Generate pip install command for missing libraries
missing = []
for category, libs in results.items():
    for lib, info in libs.items():
        if not info['installed']:
            missing.append(lib)

if missing:
    print(f"\\n=== INSTALL MISSING LIBRARIES ===")
    print(f"pip install {' '.join(missing)}")

# Output JSON for parsing
print("\\n=== JSON RESULTS ===")
print(json.dumps(results, indent=2))
`;

  try {
    const result = await e2bManager.executeWithFallback(checkCode, {
      language: 'python',
      timeoutMs: 30000
    });
    
    if (result.success) {
      console.log(result.stdout);
      
      // Parse results to create recommendations
      console.log('\n=== RECOMMENDATIONS ===\n');
      console.log('For GPQA improvements, we specifically need:');
      console.log('\n1. MUST HAVE (Math Core):');
      console.log('   - numpy (fundamental for calculations)');
      console.log('   - sympy (symbolic math)');
      console.log('   - scipy (scientific computing)');
      console.log('   - mpmath (arbitrary precision)');
      
      console.log('\n2. HIGHLY RECOMMENDED (Bio Helpers):');
      console.log('   - splicekit (RNA splicing analysis)');
      console.log('   - rnalib (RNA transcriptomics - NEW 2025)');
      console.log('   - deeptools (ChIP-seq visualization)');
      console.log('   - pandas (data manipulation)');
      console.log('   - biopython (general bioinformatics)');
      
      console.log('\n3. OPTIONAL (Chemistry):');
      console.log('   - rdkit (molecular informatics)');
      console.log('   - chempy (chemical computations)');
      
    } else {
      console.error('Failed to check libraries:', result.stderr);
    }
    
  } catch (error) {
    console.error('Error:', error.message);
  } finally {
    await e2bManager.shutdown();
  }
}

checkLibraries().catch(console.error);