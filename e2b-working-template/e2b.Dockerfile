FROM e2bdev/code-interpreter:latest

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

# Update package list and install system dependencies
RUN apt-get update && apt-get install -y \
    libxrender1 \
    libxext6 \
    libgl1-mesa-glx \
    libglib2.0-0 \
    build-essential \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Ensure pip is up to date
RUN python -m pip install --upgrade pip

# Install Python packages
RUN pip install --no-cache-dir \
    "numpy<2" \
    scipy \
    sympy \
    mpmath \
    pandas \
    matplotlib \
    seaborn \
    scikit-learn \
    rdkit \
    rdchiral \
    requests \
    beautifulsoup4 \
    lxml \
    pillow \
    networkx \
    plotly \
    statsmodels \
    yfinance

# Verify installations
RUN python -c "import numpy; print('NumPy:', numpy.__version__)" && \
    python -c "import rdkit; print('RDKit:', rdkit.__version__)" && \
    python -c "import pandas; print('Pandas:', pandas.__version__)"

# Test Python execution
RUN python -c "print('Python execution test passed')"

WORKDIR /home/user