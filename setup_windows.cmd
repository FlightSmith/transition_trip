@echo off
setlocal
where py >nul 2>&1
if errorlevel 1 goto try_python
py -3 -c "import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)" >nul 2>&1
if errorlevel 1 goto try_python
py -3 "%~dp0setup_environment.py" %*
set "SETUP_RESULT=%ERRORLEVEL%"
goto finish

:try_python
where python >nul 2>&1
if errorlevel 1 goto missing_python
python -c "import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)" >nul 2>&1
if errorlevel 1 goto missing_python
python "%~dp0setup_environment.py" %*
set "SETUP_RESULT=%ERRORLEVEL%"
goto finish

:missing_python
echo Python 3.10 or newer could not be found.
echo Install Python from https://www.python.org/downloads/ and enable its launcher or PATH option.
echo Then run setup_windows.cmd again.
set "SETUP_RESULT=2"

:finish
if not "%SETUP_RESULT%"=="0" echo Setup did not complete. Read the message above before retrying.
if "%~1"=="" pause
exit /b %SETUP_RESULT%
