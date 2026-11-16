@echo off
rem Foci Analyzer: installation from this folder, for developers and for testing. Users install with
rem FociAnalyzer-Setup.exe, which runs the same steps (installer\install.py) with a normal setup window.
rem Run it again to update or to add networks. Options are passed on, for example:
rem     install.bat --networks stardist,spotiflow --own cellpose=C:\envs\cellpose4\Scripts\python.exe
setlocal
set "HERE=%~dp0"
if not defined FOCI_ANALYZER_INSTALL_DIR set "FOCI_ANALYZER_INSTALL_DIR=%LOCALAPPDATA%\FociAnalyzer"
set "TARGET=%FOCI_ANALYZER_INSTALL_DIR%"
echo.
echo  Foci Analyzer will be installed in %TARGET%
echo.
set "UV=%TARGET%\bin\uv.exe"
if not exist "%UV%" (
    echo  Fetching uv ^(a small tool that installs Python and the packages^)...
    powershell -NoProfile -ExecutionPolicy Bypass -Command "$env:UV_INSTALL_DIR='%TARGET%\bin'; $env:UV_NO_MODIFY_PATH='1'; irm https://astral.sh/uv/install.ps1 | iex"
    if errorlevel 1 goto :failed
)
set "UV_PYTHON_INSTALL_DIR=%TARGET%\python"
set "UV_CACHE_DIR=%TARGET%\cache"
"%UV%" run --no-project --python 3.12 "%HERE%installer\install.py" --target "%TARGET%" --uv "%UV%" --shortcut --desktop %*
if errorlevel 1 goto :failed
if "%~1"=="" pause
exit /b 0

:failed
echo.
echo  The installation did not finish. The messages above say what went wrong.
if "%~1"=="" pause
exit /b 1
