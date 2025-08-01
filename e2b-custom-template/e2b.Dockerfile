FROM e2bdev/code-interpreter:latest

# Install system dependencies for RDKit
RUN apt-get update && apt-get install -y \
    libxrender1 \
    libxext6 \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /home/user

# Copy requirements file
COPY requirements.txt /tmp/requirements.txt

# Install all packages
RUN pip install --no-cache-dir -r /tmp/requirements.txt && \
    rm /tmp/requirements.txt

# Verify installations
RUN python -c "import numpy; print(f'NumPy: {numpy.__version__}')" && \
    python -c "import scipy; print(f'SciPy: {scipy.__version__}')" && \
    python -c "import pandas; print(f'Pandas: {pandas.__version__}')" && \
    python -c "from rdkit import Chem; print('RDKit: OK')" && \
    python -c "import rdchiral; print('RDChiral: OK')"

# Create a marker file
RUN echo "E2B Custom Template for Ask Any Expert" > /home/user/template_info.txt && \
    echo "Built on: $(date)" >> /home/user/template_info.txt && \
    echo "Libraries included:" >> /home/user/template_info.txt && \
    echo "- NumPy < 2 (for RDKit compatibility)" >> /home/user/template_info.txt && \
    echo "- RDKit and RDChiral (chemistry)" >> /home/user/template_info.txt && \
    echo "- SciPy, SymPy, Pandas (scientific computing)" >> /home/user/template_info.txt && \
    echo "- Matplotlib, Seaborn, Plotly (visualization)" >> /home/user/template_info.txt && \
    echo "- Scikit-learn (machine learning)" >> /home/user/template_info.txt && \
    echo "- Requests, BeautifulSoup4 (web scraping)" >> /home/user/template_info.txt

# Set Python to run in unbuffered mode
ENV PYTHONUNBUFFERED=1