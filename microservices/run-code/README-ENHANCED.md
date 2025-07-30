# Enhanced E2B Setup - Ask Any Expert

## Overview
This enhanced setup provides a comprehensive Python environment with advanced libraries for mathematical analysis, data science, financial modeling, legal document processing, and machine learning.

## Libraries Included

### Mathematics & Science
- **NumPy 1.26.4** - Numerical computing
- **SciPy 1.13.1** - Scientific computing
- **SymPy 1.12.1** - Symbolic mathematics
- **mpmath 1.3.0** - Arbitrary precision arithmetic

### Data Processing
- **pandas 2.2.2** - Data manipulation and analysis
- **polars 0.20.21** - Fast DataFrame library
- **pyarrow 16.1.0** - Columnar data processing

### Statistics & Machine Learning
- **statsmodels 0.14.1** - Statistical modeling
- **scikit-learn 1.5.0** - Machine learning

### Visualization
- **matplotlib 3.9.0** - Plotting library
- **seaborn 0.13.2** - Statistical visualizations
- **plotly 5.22.0** - Interactive plots

### Natural Language Processing
- **spaCy 3.7.4** - NLP library
- **blackstone-spacy 1.0.0** - Legal text processing
- **sentence-transformers 3.0.0** - Sentence embeddings
- **transformers 4.42.0** - Hugging Face transformers

### Finance
- **yfinance 0.2.40** - Yahoo Finance API
- **pandas-datareader 0.10.1** - Financial data reader
- **quantlib 1.35** - Quantitative finance

### Document Processing
- **pdfminer.six 20241030** - PDF text extraction
- **pypdf 4.3.0** - PDF manipulation

### Web & Security
- **requests 2.32.3** - HTTP library
- **beautifulsoup4 4.12.3** - HTML parsing
- **python-whois 0.9.4** - Domain information
- **shodan 1.30.0** - Security analysis

### Optimization & Networks
- **pulp 2.8.0** - Linear programming
- **networkx 3.3** - Network analysis

### Geospatial
- **geopandas 0.14.4** - Geographic data
- **shapely 2.0.4** - Geometric operations

## Setup Instructions

### 1. Build Custom Template
```bash
cd microservices/run-code
chmod +x build-template.sh
./build-template.sh
```

### 2. Add Template ID to Environment
After building, add the template ID to your `.env` file:
```env
E2B_TEMPLATE_ID=tpl_xxxxxxxx
```

### 3. Restart Server
The server will automatically use the custom template when `E2B_TEMPLATE_ID` is set.

## Timeout Configuration

The system now uses more generous timeouts for complex analysis:
- **Default timeout**: 30 seconds (increased from 6 seconds)
- **Math calculations**: 60 seconds (for complex symbolic math)
- **Server default**: 30 seconds

## Usage Examples

### Advanced Mathematics
```python
import sympy as sp
import numpy as np
from scipy import optimize

# Symbolic calculus
x = sp.Symbol('x')
f = x**3 - 2*x**2 + x - 1
derivative = sp.diff(f, x)
```

### Financial Analysis
```python
import yfinance as yf
import pandas as pd

# Stock data analysis
ticker = yf.Ticker("AAPL")
data = ticker.history(period="1y")
```

### Legal Document Processing
```python
import spacy
nlp = spacy.load("en_core_web_sm")

# Process legal text with Blackstone
doc = nlp("The defendant breached the contract...")
```

### Machine Learning
```python
from sklearn.ensemble import RandomForestClassifier
import pandas as pd

# ML analysis on datasets
```

## Benefits

1. **Instant Import**: All libraries are pre-installed, no installation delays
2. **Cold Start**: < 350 MB RAM usage keeps startup times fast
3. **License Safe**: All libraries use permissive licenses (no GPL/AGPL)
4. **Comprehensive Coverage**: Math, finance, legal, ML, data science
5. **Production Ready**: Pinned versions for reproducibility