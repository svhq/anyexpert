@echo off
echo Starting GPQA Testing UI with Full Observability...
echo ================================================
echo.
echo Starting E2B service on port 8001...
start cmd /k "cd microservices\run-code && node server-e2b-custom.js"

echo.
echo Waiting for E2B service to start...
timeout /t 3 /nobreak > nul

echo.
echo Starting GPQA UI server on port 3456...
node server-gpqa-ui.js

pause