@echo off
setlocal
set "TRIP_PYTHON=%~dp0.venv\Scripts\python.exe"
if not exist "%TRIP_PYTHON%" goto missing_environment
"%TRIP_PYTHON%" "%~dp0run_trip.py" %*
set "TRIP_RESULT=%ERRORLEVEL%"
goto finish

:missing_environment
echo The project Python environment is not ready.
echo Run setup_windows.cmd in this folder first, then start run_windows.cmd again.
set "TRIP_RESULT=2"

:finish
if "%~1"=="" pause
exit /b %TRIP_RESULT%
