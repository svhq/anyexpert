const axios = require('axios');
const config = require('./config');

async function verifyE2BSetup() {
    console.log('🔍 Verifying E2B Setup and Configuration');
    console.log('=====================================\n');

    // Check configuration
    console.log('📋 Configuration Check:');
    console.log(`✓ E2B API Key: ${config.E2B_API_KEY ? 'Configured' : '❌ Missing'}`);
    console.log(`✓ E2B Service URL: ${config.codeExecution.serviceUrl}`);
    console.log(`✓ E2B Template ID: ${config.E2B_TEMPLATE_ID || 'Using free template'}`);
    console.log(`✓ E2B Timeout: ${config.E2B_TIMEOUT}ms`);
    
    // Test E2B service connection
    console.log('\n🔌 Testing E2B Service Connection...');
    try {
        const healthResponse = await axios.get(`${config.codeExecution.serviceUrl}/health`);
        console.log('✓ E2B service is running:', healthResponse.data);
    } catch (error) {
        console.error('❌ E2B service not reachable:', error.message);
        console.log('Please start the E2B service with: npm run start:e2b');
        return;
    }

    // Test Python with all libraries
    console.log('\n🐍 Testing Python Libraries...');
    const testCode = `
import sys
print(f"Python version: {sys.version}")

# Test all scientific libraries
libraries = {
    'numpy': 'np',
    'scipy': 'scipy',
    'sympy': 'sp',
    'pandas': 'pd',
    'matplotlib': 'plt',
    'decimal': 'decimal',
    'fractions': 'fractions',
    'mpmath': 'mpmath',
    'statistics': 'statistics'
}

for lib, alias in libraries.items():
    try:
        if lib == 'matplotlib':
            exec(f"import matplotlib.pyplot as {alias}")
        else:
            exec(f"import {lib} as {alias}" if alias != lib else f"import {lib}")
        print(f"✓ {lib:<15} - Imported successfully")
    except ImportError as e:
        print(f"❌ {lib:<15} - {str(e)}")

# Quick functionality test
import numpy as np
print(f"\\nNumPy test: np.pi = {np.pi}")
print(f"Array ops: np.array([1,2,3]) * 2 = {np.array([1,2,3]) * 2}")
`;

    try {
        const response = await axios.post(`${config.codeExecution.serviceUrl}/run`, {
            language: 'python',
            code: testCode,
            timeout: 30000
        });

        console.log('\n📦 Library Installation Status:');
        console.log(response.data.output);
        
        if (response.data.exitCode === 0) {
            console.log('\n✅ All libraries are working correctly!');
        } else {
            console.log('\n⚠️ Some libraries may need installation');
        }

    } catch (error) {
        console.error('❌ Error testing libraries:', error.message);
    }

    // Test a sample GPQA-style calculation
    console.log('\n🧮 Testing GPQA-style Calculation...');
    const gpqaTestCode = `
import numpy as np
from scipy import constants

# Example: Quantum mechanics calculation
hbar = constants.hbar
print(f"ℏ (reduced Planck constant) = {hbar} J⋅s")

# Example: Uncertainty principle calculation
delta_t = 1e-9  # 1 nanosecond
delta_E = hbar / (2 * delta_t)
delta_E_eV = delta_E / constants.e  # Convert to eV
print(f"\\nFor Δt = {delta_t} s:")
print(f"Minimum ΔE = {delta_E_eV:.2e} eV")
`;

    try {
        const response = await axios.post(`${config.codeExecution.serviceUrl}/run`, {
            language: 'python',
            code: gpqaTestCode,
            timeout: 10000
        });

        console.log('Result:');
        console.log(response.data.output);
        
    } catch (error) {
        console.error('❌ Error running calculation:', error.message);
    }

    console.log('\n✨ Verification Complete!');
}

// Run verification
verifyE2BSetup().catch(console.error);