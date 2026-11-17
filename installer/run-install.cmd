@echo off
rem Called by FociAnalyzer-Setup.exe after it copied the files: runs installer\install.py with the uv that came
rem with the setup. The program folder is two levels up from this file ({app}\program\installer).
rem The setup runs it hidden with --progress FILE and shows the messages itself: then everything is also
rem written to install-output.log (so a failure before install.py starts can be shown), and it must not pause.
setlocal
for %%I in ("%~dp0..\..") do set "TARGET=%%~fI"
set "UV=%TARGET%\bin\uv.exe"
set "UV_PYTHON_INSTALL_DIR=%TARGET%\python"
set "UV_CACHE_DIR=%TARGET%\cache"
set "UV_NO_PROGRESS=1"
set "PYTHONIOENCODING=utf-8"
title Foci Analyzer: installing
rem a Python that is only partly there (an uninstall while the program ran) is removed: uv would take it as installed
if exist "%TARGET%\python" for /d %%P in ("%TARGET%\python\cpython-*") do if not exist "%%P\Lib\os.py" rmdir /s /q "%%P"
set "HIDDEN="
rem the full path: with Git's tools on PATH, "find" is the Unix find
echo %* | "%SystemRoot%\System32\find.exe" "--progress" >nul && set "HIDDEN=1"
if defined HIDDEN goto hidden
"%UV%" run --no-project --python 3.12 "%~dp0install.py" --target "%TARGET%" --uv "%UV%" %*
if errorlevel 1 goto failed
exit /b 0

:hidden
"%UV%" run --no-project --python 3.12 "%~dp0install.py" --target "%TARGET%" --uv "%UV%" %* > "%TARGET%\install-output.log" 2>&1
set "CODE=%errorlevel%"
rem the setup waits for this file (it does not wait for this script: that would freeze its window)
> "%TARGET%\install-exit.txt" echo %CODE%
exit /b %CODE%

:failed
echo.
echo  The installation did not finish. The messages above say what went wrong; the full log is
echo  %TARGET%\install.log
pause
exit /b 1
