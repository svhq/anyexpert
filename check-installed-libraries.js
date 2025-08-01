const fetch = require('node-fetch');

async function checkInstalledLibraries() {
  console.log('Checking installed libraries in E2B environment...\n');
  
  const code = `
import subprocess
import sys

print("=== Python Environment ===")
print(f"Python: {sys.version}")
print(f"Executable: {sys.executable}")

print("\\n=== Installed Packages ===")
result = subprocess.run([sys.executable, "-m", "pip", "list"], capture_output=True, text=True)

# Parse and display packages
if result.returncode == 0:
    lines = result.stdout.strip().split('\\n')
    packages = []
    
    # Skip header lines
    for line in lines[2:]:  # Skip "Package Version" and "------- -------"
        if line.strip():
            parts = line.split()
            if len(parts) >= 2:
                packages.append((parts[0], parts[1]))
    
    # Sort packages alphabetically
    packages.sort()
    
    # Group packages by category
    categories = {
        'Core Scientific': ['numpy', 'scipy', 'pandas', 'matplotlib', 'seaborn', 'scikit-learn'],
        'Math & Symbolic': ['sympy', 'mpmath'],
        'Chemistry': ['rdkit', 'rdchiral'],
        'Data Processing': ['pyarrow', 'polars', 'openpyxl', 'xlrd'],
        'Web & Network': ['requests', 'beautifulsoup4', 'lxml', 'urllib3'],
        'Image Processing': ['pillow', 'opencv-python'],
        'Geospatial': ['geopandas', 'shapely', 'folium'],
        'NLP': ['spacy', 'nltk', 'transformers'],
        'Finance': ['yfinance', 'pandas-datareader'],
        'Statistics': ['statsmodels'],
        'Plotting': ['plotly', 'bokeh'],
        'Other': []
    }
    
    # Categorize packages
    categorized = {cat: [] for cat in categories}
    
    for pkg, version in packages:
        pkg_lower = pkg.lower()
        categorized_flag = False
        
        for category, keywords in categories.items():
            if category != 'Other':
                for keyword in keywords:
                    if keyword.lower() in pkg_lower or pkg_lower in keyword.lower():
                        categorized[category].append(f"{pkg} ({version})")
                        categorized_flag = True
                        break
                if categorized_flag:
                    break
        
        if not categorized_flag:
            categorized['Other'].append(f"{pkg} ({version})")
    
    # Display by category
    for category, pkgs in categorized.items():
        if pkgs:
            print(f"\\n{category}:")
            for pkg in pkgs:
                print(f"  - {pkg}")
    
    print(f"\\nTotal packages: {len(packages)}")
else:
    print("Error getting package list:", result.stderr)

# Check specific important packages
print("\\n=== Testing Key Imports ===")
test_imports = [
    ('numpy', 'np'),
    ('scipy', None),
    ('pandas', 'pd'),
    ('matplotlib.pyplot', 'plt'),
    ('sympy', None),
    ('rdkit.Chem', 'Chem'),
    ('rdchiral', None),
    ('requests', None),
    ('sklearn', None)
]

for module, alias in test_imports:
    try:
        if alias:
            exec(f"import {module} as {alias}")
            exec(f"print(f'✓ {module}: {{str({alias}.__version__)[:10] if hasattr({alias}, \"__version__\") else \"OK\"}}')")
        else:
            exec(f"import {module}")
            exec(f"print(f'✓ {module}: {{str({module}.__version__)[:10] if hasattr({module}, \"__version__\") else \"OK\"}}')")
    except ImportError:
        print(f"✗ {module}: Not installed")
    except Exception as e:
        print(f"✗ {module}: Error - {type(e).__name__}")
`;

  try {
    const response = await fetch('http://localhost:3001/run_code', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        language: 'python',
        source: code,
        timeout: 30000
      })
    });

    const result = await response.json();
    
    if (result.stdout) {
      console.log(result.stdout);
    }
    
    if (result.stderr && result.stderr.trim()) {
      console.log('\nErrors/Warnings:', result.stderr.substring(0, 200) + '...');
    }
    
  } catch (error) {
    console.error('Error:', error.message);
  }
}

checkInstalledLibraries();