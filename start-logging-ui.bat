@echo off
echo Starting Logging UI Server...
echo.
echo UI will be available at http://localhost:3004
echo Logs will be saved to realtime-ui\session-logs\
echo.
cd realtime-ui
node server-with-logging.js