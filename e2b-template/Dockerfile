# Use the default E2B base image
FROM e2bdev/code-interpreter:latest

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Update system packages and install dependencies
RUN apt-get update -y && \
    apt-get install -y \
    build-essential \
    git \
    curl \
    ca-certificates \
    libgl1 \
    libglib2.0-0 \
    tzdata \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip and install Python scientific libraries
RUN python -m pip install --upgrade pip && \
    pip install --no-cache-dir \
    numpy \
    pandas \
    sympy \
    scipy \
    matplotlib \
    seaborn \
    plotly \
    networkx \
    scikit-learn \
    statsmodels \
    requests \
    beautifulsoup4 \
    lxml \
    pillow \
    mpmath \
    jupyter \
    ipython

# Verify installations
RUN python -c "import numpy, pandas, sympy, scipy, matplotlib, seaborn, plotly, networkx, sklearn, statsmodels, requests, bs4, lxml, PIL, mpmath; print('All libraries installed successfully')"

# Set working directory
WORKDIR /tmp