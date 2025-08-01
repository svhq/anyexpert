import sys
print('Python version:', sys.version)
print('\nTesting library imports:')

# Test imports
libraries = [
    ('numpy', 'np'),
    ('pandas', 'pd'),
    ('sympy', 'sp'),
    ('scipy', None),
    ('matplotlib', None),
    ('sklearn', None),
    ('statsmodels', None),
    ('plotly', None),
    ('networkx', None),
    ('bs4', None),
    ('lxml', None),
    ('PIL', None),
    ('mpmath', None),
    ('requests', None),
    ('seaborn', None)
]

for lib, alias in libraries:
    try:
        if alias:
            exec(f'import {lib} as {alias}')
            exec(f'print("✓ {lib} ({alias}.__version__):", {alias}.__version__)')
        else:
            exec(f'import {lib}')
            if hasattr(eval(lib), '__version__'):
                exec(f'print("✓ {lib}:", {lib}.__version__)')
            else:
                print(f"✓ {lib}: imported successfully")
    except ImportError as e:
        print(f"✗ {lib}: {e}")
    except Exception as e:
        print(f"✗ {lib}: {type(e).__name__}: {e}")

print('\nTesting some functionality:')
try:
    import numpy as np
    arr = np.array([1, 2, 3])
    print(f"numpy array: {arr}")
    print(f"numpy mean: {np.mean(arr)}")
except Exception as e:
    print(f"numpy test failed: {e}")

try:
    import sympy as sp
    x = sp.Symbol('x')
    expr = x**2 + 2*x + 1
    print(f"sympy expression: {expr}")
    print(f"sympy factor: {sp.factor(expr)}")
except Exception as e:
    print(f"sympy test failed: {e}")