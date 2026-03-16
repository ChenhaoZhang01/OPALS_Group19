@echo off
setlocal
cd /d "%~dp0"

set "OUTPUT=metadata\master_metadata.csv"
set "PER_QUERY_LIMIT=5000"
set "SHOULD_PAUSE=1"

if not "%~1"=="" set "OUTPUT=%~1"
if not "%~2"=="" set "PER_QUERY_LIMIT=%~2"
if /i "%~3"=="--no-pause" set "SHOULD_PAUSE=0"

echo ===============================================
echo OPALS Phase 1: Metadata Quick Start
echo ===============================================
echo Output file: %OUTPUT%
echo Query limit: %PER_QUERY_LIMIT%
echo.

powershell -NoProfile -ExecutionPolicy Bypass -File ".\scripts\quickstart_phase1.ps1" -Output "%OUTPUT%" -PerQueryLimit %PER_QUERY_LIMIT%
if errorlevel 1 (
    echo.
    echo ERROR: Setup failed. Please read the message above.
    if "%SHOULD_PAUSE%"=="1" pause
    exit /b 1
)

echo.
echo Success. Metadata is ready at: %OUTPUT%
if "%SHOULD_PAUSE%"=="1" pause
exit /b 0
