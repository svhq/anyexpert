# Use the official E2B base image which has Python already configured
FROM e2bdev/code-interpreter:latest

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Update pip first
RUN python -m pip install --upgrade pip

# Install core scientific libraries step by step with error checking
RUN echo "Installing numpy..." && \
    pip install --no-cache-dir numpy && \
    python -c "import numpy as np; print(f'✓ NumPy {np.__version__} installed')"

RUN echo "Installing pandas..." && \
    pip install --no-cache-dir pandas && \
    python -c "import pandas as pd; print(f'✓ Pandas {pd.__version__} installed')"

RUN echo "Installing sympy..." && \
    pip install --no-cache-dir sympy && \
    python -c "import sympy as sp; print(f'✓ SymPy {sp.__version__} installed')"

RUN echo "Installing scipy..." && \
    pip install --no-cache-dir scipy && \
    python -c "import scipy; print(f'✓ SciPy {scipy.__version__} installed')"

RUN echo "Installing matplotlib..." && \
    pip install --no-cache-dir matplotlib && \
    python -c "import matplotlib; print(f'✓ Matplotlib {matplotlib.__version__} installed')"

# Install additional useful libraries
RUN pip install --no-cache-dir \
    requests \
    beautifulsoup4 \
    lxml \
    pillow \
    mpmath \
    seaborn \
    plotly \
    networkx

# Install ML libraries (may take longer)
RUN pip install --no-cache-dir \
    scikit-learn \
    statsmodels

# Final verification - test core imports only
RUN python -c "import numpy as np, pandas as pd, sympy as sp, scipy; print('✅ Core scientific libraries installed!'); print(f'NumPy: {np.__version__}'); print(f'Pandas: {pd.__version__}'); print(f'SymPy: {sp.__version__}'); print(f'SciPy: {scipy.__version__}')"

# Set working directory
WORKDIR /tmp
