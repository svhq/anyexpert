#!/bin/bash
# Build script for E2B custom template

echo "Building E2B template with advanced Python libraries..."

# Build and push the template
e2b template build -f e2b.Dockerfile -t askanyexpert-advanced

echo ""
echo "Template built successfully!"
echo "Copy the template ID from above and add it to your .env file as:"
echo "E2B_TEMPLATE_ID=tpl_xxxxxxxx"
echo ""
echo "Then restart your server to use the new template."