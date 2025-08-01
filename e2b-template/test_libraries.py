#!/usr/bin/env python3
"""
Test script to verify all scientific libraries are properly installed
"""

import sys
print(f"Python version: {sys.version}")
print("=" * 50)

# Test core scientific libraries
libraries = [
    ('numpy', 'np'),
    ('pandas', 'pd'), 
    ('sympy', 'sp'),
    ('scipy', None),
    ('matplotlib.pyplot', 'plt'),
    ('seaborn', 'sns'),
    ('plotly', None),
    ('networkx', 'nx'),
    ('sklearn', None),
    ('statsmodels', None),
    ('requests', None),
    ('bs4', None),
    ('lxml', None),
    ('PIL', None),
    ('mpmath', None)
]

print("Testing library imports:")
for lib, alias in libraries:
    try:
        if alias:
            exec(f'import {lib} as {alias}')
            if hasattr(eval(alias), '__version__'):
                version = eval(f'{alias}.__version__')
                print(f"✅ {lib} ({alias}): {version}")
            else:
                print(f"✅ {lib} ({alias}): imported successfully")
        else:
            exec(f'import {lib}')
            if hasattr(eval(lib), '__version__'):
                version = eval(f'{lib}.__version__')
                print(f"✅ {lib}: {version}")
            else:
                print(f"✅ {lib}: imported successfully")
    except ImportError as e:
        print(f"❌ {lib}: {e}")
    except Exception as e:
        print(f"⚠️  {lib}: {type(e).__name__}: {e}")

print("=" * 50)

# Test some basic functionality
print("Testing basic functionality:")

try:
    import numpy as np
    arr = np.array([1, 2, 3, 4, 5])
    print(f"✅ NumPy: mean([1,2,3,4,5]) = {np.mean(arr)}")
except Exception as e:
    print(f"❌ NumPy test failed: {e}")

try:
    import sympy as sp
    x = sp.Symbol('x')
    expr = x**2 - 4
    roots = sp.solve(expr, x)
    print(f"✅ SymPy: roots of x²-4 = {roots}")
except Exception as e:
    print(f"❌ SymPy test failed: {e}")

try:
    import pandas as pd
    df = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})
    print(f"✅ Pandas: DataFrame shape = {df.shape}")
except Exception as e:
    print(f"❌ Pandas test failed: {e}")

print("=" * 50)
print("🎉 Library test complete!")