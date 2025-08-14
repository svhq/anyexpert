  # E2B Unified Template for Ask Any Expert
  # Version: 2.0.0

  FROM e2bdev/code-interpreter:latest

  ENV DEBIAN_FRONTEND=noninteractive
  ENV PYTHONUNBUFFERED=1

  # Install system dependencies
  RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential \
      libxml2-dev \
      libxslt-dev \
      libgeos-dev \
      libgdal-dev \
      libproj-dev \
      libspatialindex-dev \
      libxrender1 \
      libxext6 \
      git \
      && rm -rf /var/lib/apt/lists/*

  # Upgrade pip
  RUN python -m pip install --upgrade pip setuptools wheel

  # Install all libraries
  RUN pip install --no-cache-dir \
      numpy==1.26.4 \
      scipy==1.13.1 \
      pandas==2.2.2 \
      sympy==1.12.1 \
      mpmath==1.3.0 \
      numpy-financial==1.0.0 \
      statsmodels==0.14.1 \
      yfinance==0.2.40 \
      pandas-datareader==0.10.0 \
      matplotlib==3.9.0 \
      seaborn==0.13.2 \
      plotly==5.22.0 \
      pillow==10.3.0 \
      scikit-learn==1.5.0 \
      transformers==4.42.0 \
      sentence-transformers==3.0.0 \
      spacy==3.7.4 \
      biopython==1.83 \
      rdkit \
      geopandas==0.14.4 \
      shapely==2.0.4 \
      networkx==3.3 \
      beautifulsoup4==4.12.3 \
      requests==2.32.3 \
      lxml==5.2.1 \
      pdfminer.six==20231228 \
      pypdf==4.3.0 \
      QuantLib \
      python-dateutil \
      openpyxl \
      xlrd \
      polars==0.20.21 \
      pyarrow==16.1.0 \
      pulp==2.8.0

  # Download spaCy model
  RUN python -m spacy download en_core_web_sm

  WORKDIR /home/user

  RUN echo "Ask Any Expert E2B Template Ready!"
