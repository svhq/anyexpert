  # E2B Unified Template for Ask Any Expert
  # Version: 2.0.0
  # This template includes all libraries needed for the Ask Any Expert system

  FROM e2bdev/code-interpreter:latest

  # Set environment variables
  ENV DEBIAN_FRONTEND=noninteractive
  ENV PYTHONUNBUFFERED=1

  # Install system dependencies required by various Python packages
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

  # Upgrade pip and install build tools
  RUN python -m pip install --upgrade pip setuptools wheel

  # Install core scientific computing libraries
  RUN pip install --no-cache-dir \
      numpy==1.26.4 \
      scipy==1.13.1 \
      pandas==2.2.2 \
      sympy==1.12.1 \
      mpmath==1.3.0

  # Install financial and statistical libraries
  RUN pip install --no-cache-dir \
      numpy-financial==1.0.0 \
      statsmodels==0.14.1 \
      yfinance==0.2.40 \
      pandas-datareader==0.10.0

  # Note: QuantLib requires special installation
  RUN pip install --no-cache-dir QuantLib

  # Install visualization libraries
  RUN pip install --no-cache-dir \
      matplotlib==3.9.0 \
      seaborn==0.13.2 \
      plotly==5.22.0 \
      pillow==10.3.0

  # Install machine learning libraries (these are large, install separately)
  RUN pip install --no-cache-dir scikit-learn==1.5.0

  # Install transformers and NLP libraries (very large, separate installation)
  RUN pip install --no-cache-dir \
      transformers==4.42.0 \
      sentence-transformers==3.0.0

  # Install spaCy and download English model
  RUN pip install --no-cache-dir spacy==3.7.4 && \
      python -m spacy download en_core_web_sm

  # Install domain-specific libraries
  RUN pip install --no-cache-dir biopython==1.83

  # Install RDKit for chemistry (requires specific handling)
  RUN pip install --no-cache-dir rdkit

  # Install geospatial libraries
  RUN pip install --no-cache-dir \
      geopandas==0.14.4 \
      shapely==2.0.4

  # Install network analysis
  RUN pip install --no-cache-dir networkx==3.3

  # Install web scraping and data extraction
  RUN pip install --no-cache-dir \
      beautifulsoup4==4.12.3 \
      requests==2.32.3 \
      lxml==5.2.1 \
      pdfminer.six==20241030 \
      pypdf==4.3.0

  # Install additional useful libraries
  RUN pip install --no-cache-dir \
      python-dateutil \
      openpyxl \
      xlrd \
      polars==0.20.21 \
      pyarrow==16.1.0 \
      pulp==2.8.0

  # Create verification script
  RUN cat > /tmp/verify_libraries.py << 'EOF'
  #!/usr/bin/env python3
  import sys
  import importlib

  def verify_library(name, import_name=None):
      """Verify a library is installed and get its version"""
      import_name = import_name or name
      try:
          module = importlib.import_module(import_name)
          version = getattr(module, '__version__', 'installed')
          print(f"✅ {name:25} {version}")
          return True
      except ImportError as e:
          print(f"❌ {name:25} MISSING - {e}")
          return False

  print("=" * 60)
  print("Ask Any Expert - E2B Template Library Verification")
  print("=" * 60)

  # Core libraries
  print("\n📊 Core Scientific Computing:")
  verify_library("numpy")
  verify_library("scipy")
  verify_library("pandas")
  verify_library("sympy")
  verify_library("mpmath")

  print("\n💰 Financial & Statistical:")
  verify_library("numpy_financial")
  verify_library("statsmodels")
  verify_library("yfinance")
  verify_library("QuantLib")

  print("\n📈 Visualization:")
  verify_library("matplotlib")
  verify_library("seaborn")
  verify_library("plotly")
  verify_library("PIL", "PIL")

  print("\n🤖 Machine Learning & AI:")
  verify_library("sklearn", "sklearn")
  verify_library("transformers")
  verify_library("sentence_transformers")
  verify_library("spacy")

  print("\n🧬 Domain-Specific:")
  verify_library("biopython", "Bio")
  verify_library("rdkit", "rdkit")
  verify_library("geopandas")
  verify_library("networkx")

  print("\n🌐 Web & Data:")
  verify_library("beautifulsoup4", "bs4")
  verify_library("requests")
  verify_library("lxml")
  verify_library("pdfminer", "pdfminer")

  print("\n" + "=" * 60)
  print("Template verification complete!")
  print("=" * 60)
  EOF

  # Run verification
  RUN python /tmp/verify_libraries.py

  # Create a marker file with template information
  RUN echo "Ask Any Expert E2B Template v2.0.0" > /tmp/template_info.txt && \
      echo "Build Date: $(date)" >> /tmp/template_info.txt && \
      echo "Python Version: $(python --version)" >> /tmp/template_info.txt && \
      echo "Total Packages: $(pip list | wc -l)" >> /tmp/template_info.txt

  # Set working directory
  WORKDIR /home/user

  # Final message
  RUN echo "✅ Ask Any Expert E2B Template Ready!"
