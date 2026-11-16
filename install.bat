@echo off
rem Foci Analyzer: installation. Double-click this file; run it again to update or to add the networks.
rem Everything goes into one folder (default %LOCALAPPDATA%\FociAnalyzer); removing that folder and the
rem desktop shortcut removes the program. No Python needs to be installed beforehand.
setlocal EnableDelayedExpansion
set "HERE=%~dp0"
if not defined FOCI_ANALYZER_INSTALL_DIR set "FOCI_ANALYZER_INSTALL_DIR=%LOCALAPPDATA%\FociAnalyzer"
set "TARGET=%FOCI_ANALYZER_INSTALL_DIR%"
echo.
echo  Foci Analyzer will be installed in %TARGET%
echo.

rem ---- 1. uv, the tool that fetches Python and the packages
set "UV=uv"
where uv >nul 2>nul
if errorlevel 1 (
    if exist "%USERPROFILE%\.local\bin\uv.exe" (
        set "UV=%USERPROFILE%\.local\bin\uv.exe"
    ) else (
        echo  Fetching uv ^(a small tool that installs Python and the packages^)...
        powershell -NoProfile -ExecutionPolicy Bypass -Command "irm https://astral.sh/uv/install.ps1 | iex"
        if errorlevel 1 goto :failed
        set "UV=%USERPROFILE%\.local\bin\uv.exe"
    )
)

rem ---- 2. the program itself
echo  Installing Foci Analyzer...
set "UV_PROJECT_ENVIRONMENT=%TARGET%\.venv"
"%UV%" sync --project "%HERE%." --no-dev --python 3.12
if errorlevel 1 goto :failed

rem ---- 3. the neural networks (optional)
if exist "%TARGET%\models-env\Scripts\python.exe" (
    echo  The neural networks are already installed.
    set "CHOICE=3"
) else if defined FOCI_ANALYZER_NETWORKS (
    set "CHOICE=%FOCI_ANALYZER_NETWORKS%"
) else (
    echo.
    echo  Neural networks: Cellpose and InstanSeg ^(nuclei and cells^), Spotiflow ^(foci^), Trackastra ^(tracking^).
    echo    1  Install them ^(about 6 GB; recommended^)
    echo    2  Not now. Threshold segmentation and watershed foci work without them, and an existing
    echo       environment can be given in the window under Nuclei / cells - Advanced.
    set /p "CHOICE=  Your choice [1]: "
    if not defined CHOICE set "CHOICE=1"
)
if "%CHOICE%"=="1" (
    echo  Installing the networks; this takes a while...
    "%UV%" venv "%TARGET%\models-env" --python 3.12
    if errorlevel 1 goto :failed
    set "MPY=%TARGET%\models-env\Scripts\python.exe"
    where nvidia-smi >nul 2>nul
    if errorlevel 1 (
        echo  No NVIDIA graphics card found: the networks will run on the processor ^(slower^).
        "%UV%" pip install --python "!MPY!" torch torchvision --index-url https://download.pytorch.org/whl/cpu
    ) else (
        "%UV%" pip install --python "!MPY!" torch torchvision --index-url https://download.pytorch.org/whl/cu124
    )
    if errorlevel 1 goto :failed
    "%UV%" pip install --python "!MPY!" "cellpose>=3.1,<4" spotiflow instanseg-torch trackastra
    if errorlevel 1 goto :failed
)

rem ---- 4. a shortcut on the desktop
if defined FOCI_ANALYZER_NO_SHORTCUT goto :ready
powershell -NoProfile -ExecutionPolicy Bypass -Command "$s = (New-Object -ComObject WScript.Shell).CreateShortcut([IO.Path]::Combine([Environment]::GetFolderPath('Desktop'), 'Foci Analyzer.lnk')); $s.TargetPath = '%TARGET%\.venv\Scripts\foci-analyzer-window.exe'; $s.WorkingDirectory = '%TARGET%'; $s.Description = 'Foci Analyzer'; $s.Save()"
echo.
echo.
:ready
echo  Ready. Start Foci Analyzer with the shortcut on the desktop; it opens in the browser.
echo  Stop it with the Quit button in the window.
echo.
if not defined FOCI_ANALYZER_NETWORKS pause
exit /b 0

:failed
echo.
echo  The installation did not finish. The message above says what went wrong; an internet connection is needed.
if not defined FOCI_ANALYZER_NETWORKS pause
exit /b 1
