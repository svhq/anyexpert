# Library Installation Status

## Current Situation

### E2B Environment
- **Current Mode**: Dynamic installation (no custom template)
- **Base Environment**: Minimal (only pip, setuptools, wheel)
- **Template Status**: Created but not properly built with libraries

### Dynamic Library Installation
When code is executed, the E2B service (`server-e2b-custom.js`) automatically installs:

#### Scientific Libraries (triggered by numpy, scipy, pandas, etc.)
- numpy
- scipy
- sympy
- pandas
- matplotlib

#### Chemistry Libraries (triggered by rdkit, Chem, etc.)
- numpy<2 (for RDKit compatibility)
- rdkit-pypi
- rdchiral

## Custom Template Issue

We created a custom template (`u8w4jzut60mgigs4ofp9`) but it didn't preserve the installed libraries because:

1. The E2B template build process requires specific steps to persist packages
2. Our Dockerfile installed packages but they weren't preserved in the final image
3. The template build process might require additional configuration

## Current Working Solution

The system currently works with dynamic installation:
1. Each new E2B sandbox starts clean
2. Libraries are installed on-demand based on code content
3. Installation adds ~30-60 seconds to first execution
4. Subsequent runs in the same sandbox are fast

## Libraries Available Through Dynamic Installation

### Core Scientific Computing
- **numpy** (<2 for chemistry, latest for others)
- **scipy** - Scientific computing
- **sympy** - Symbolic mathematics
- **pandas** - Data manipulation
- **matplotlib** - Plotting

### Chemistry
- **rdkit-pypi** - Cheminformatics
- **rdchiral** - Chirality analysis

### Additional Libraries (can be added to installation)
- **scikit-learn** - Machine learning
- **requests** - HTTP requests
- **beautifulsoup4** - Web scraping
- **lxml** - XML/HTML processing
- **pillow** - Image processing
- **networkx** - Graph algorithms
- **plotly** - Interactive plots
- **seaborn** - Statistical plotting
- **statsmodels** - Statistical models
- **yfinance** - Financial data

## Math Agent Capabilities

The Math Agent has access to:
1. **Python execution** via E2B
2. **All dynamically installed libraries**
3. **Chemistry tools** (RDKit, RDChiral)
4. **Mathematical tools** (SymPy, NumPy, SciPy)
5. **Data analysis** (Pandas, Matplotlib)

## To Fix Custom Template

To properly create a custom template with pre-installed libraries:

1. Use a different approach in Dockerfile:
```dockerfile
# Install in user space
RUN pip install --user numpy<2 rdkit rdchiral scipy sympy pandas matplotlib

# Or use a requirements file
COPY requirements.txt /tmp/
RUN pip install -r /tmp/requirements.txt
```

2. Follow E2B's specific template guidelines for persisting packages

3. Test the template thoroughly before deployment

## Current Recommendation

Continue using dynamic installation because:
- It works reliably
- Libraries are always up-to-date
- No template maintenance needed
- Only installs what's needed
- Installation time is acceptable (30-60s on first use)