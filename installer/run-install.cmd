@echo off
rem Called by FociAnalyzer-Setup.exe after it copied the files: runs installer\install.py with the uv that came
rem with the setup. The program folder is two levels up from this file ({app}\program\installer).
rem The setup runs it hidden with --progress FILE and shows the messages itself: then it must not pause.
setlocal
for %%I in ("%~dp0..\..") do set "TARGET=%%~fI"
set "UV=%TARGET%\bin\uv.exe"
set "UV_PYTHON_INSTALL_DIR=%TARGET%\python"
set "UV_CACHE_DIR=%TARGET%\cache"
set "UV_NO_PROGRESS=1"
set "PYTHONIOENCODING=utf-8"
title Foci Analyzer: installing
"%UV%" run --no-project --python 3.12 "%~dp0install.py" --target "%TARGET%" --uv "%UV%" %*
if errorlevel 1 goto failed
exit /b 0

:failed
echo %* | find "--progress" >nul && exit /b 1
echo.
echo  The installation did not finish. The messages above say what went wrong; the full log is
echo  %TARGET%\install.log
pause
exit /b 1
