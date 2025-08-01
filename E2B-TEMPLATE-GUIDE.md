# E2B Custom Template Guide

## Overview
We've created a custom E2B template (`u8w4jzut60mgigs4ofp9`) that includes all required libraries pre-installed with correct versions for chemistry problems.

## Template Details
- **Template ID**: `u8w4jzut60mgigs4ofp9`
- **Name**: `askanyexpert-numpy1`
- **Key Features**:
  - NumPy < 2 (for RDKit compatibility)
  - RDKit and RDChiral pre-installed
  - All scientific libraries (SciPy, SymPy, Pandas, Matplotlib)
  - No dynamic installation needed

## Configuration
The template ID is already configured in `.env`:
```
E2B_TEMPLATE_ID=u8w4jzut60mgigs4ofp9
```

## To Use the Custom Template

1. **Stop the current E2B service**:
   ```bash
   # Find the process using port 3001
   netstat -ano | findstr :3001
   # Kill the process (replace PID with actual process ID)
   taskkill /PID <PID> /F
   ```

2. **Start the custom E2B service**:
   ```bash
   cd microservices/run-code
   node server-e2b-custom.js
   ```

3. **Verify template is loaded**:
   ```bash
   curl http://localhost:3001/health
   ```
   Should show: `"templateId":"u8w4jzut60mgigs4ofp9"`

## Testing the Template

Run the test suite:
```bash
node test-simple.js
```

Expected output with custom template:
- Template info should show: "Custom E2B Template v3.0"
- NumPy version should be < 2
- No RDKit compatibility warnings

## Benefits
1. **Faster execution**: No library installation needed
2. **No compatibility issues**: NumPy version pinned for RDKit
3. **Consistent environment**: All libraries pre-installed
4. **Better for chemistry**: RDKit and RDChiral ready to use

## Troubleshooting
- If template not loading, check E2B_TEMPLATE_ID in .env
- Ensure using server-e2b-custom.js, not server-e2b-working.js
- Template requires E2B API key to be configured