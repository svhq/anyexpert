# e2b.Dockerfile - Enhanced Python environment for Ask Any Expert
FROM e2bdev/code-interpreter:latest            # base Debian slim + Python 3.11

# System deps for lxml & shapely wheels (small)
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential libxml2-dev libxslt-dev libgeos-dev \
    && rm -rf /var/lib/apt/lists/*

# Python packages (pin for reproducibility)
RUN pip install --no-cache-dir \
      numpy==1.26.4 scipy==1.13.1 mpmath==1.3.0 sympy==1.12.1 \
      pandas==2.2.2 polars==0.20.21 pyarrow==16.1.0 \
      statsmodels==0.14.1 scikit-learn==1.5.0 \
      matplotlib==3.9.0 seaborn==0.13.2 plotly==5.22.0 \
      spacy==3.7.4 blackstone-spacy==1.0.0 sentence-transformers==3.0.0 \
      transformers==4.42.0 \
      pdfminer.six==20241030 pypdf==4.3.0 \
      yfinance==0.2.40 pandas-datareader==0.10.1 quantlib==1.35 \
      pulp==2.8.0 networkx==3.3 \
      python-whois==0.9.4 shodan==1.30.0 \
      requests==2.32.3 beautifulsoup4==4.12.3 lxml==5.2.1 \
      pillow==10.3.0 geopandas==0.14.4 shapely==2.0.4

# Pre-download spaCy model so runtime is instant
RUN python -m spacy download en_core_web_sm