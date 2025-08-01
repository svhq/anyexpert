@echo off
echo ===================================
echo Starting Ask Any Expert System
echo ===================================
echo.

echo [1/3] Checking E2B service...
call node auto-start-e2b.js

if %errorlevel% neq 0 (
    echo Failed to start E2B service!
    echo Starting manually...
    start "E2B Service" cmd /k node start-e2b-service.js
    timeout /t 5 /nobreak > nul
)

echo.
echo [2/3] Starting main application...
echo.

echo ===================================
echo Ask Any Expert System Ready!
echo ===================================
echo.
echo Available commands:
echo - Type your questions directly
echo - Use Ctrl+C to exit
echo.

node index.js

pause