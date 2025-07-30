@echo off
echo Killing existing Node processes on port 3001...
for /f "tokens=5" %%a in ('netstat -aon ^| find ":3001" ^| find "LISTENING"') do taskkill /F /PID %%a 2>nul
timeout /t 2 /nobreak >nul
echo Starting E2B server...
cd microservices\run-code
node server-e2b-working.js