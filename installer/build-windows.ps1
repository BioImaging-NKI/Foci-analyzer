# Builds installer\build\FociAnalyzer-Setup-<version>.exe. Needs Inno Setup 6 (ISCC.exe) and uv.exe (the one
# on PATH, or the path in $env:UV_EXE; the GitHub build downloads it). Run from anywhere:
#     powershell -ExecutionPolicy Bypass -File installer\build-windows.ps1
$ErrorActionPreference = "Stop"
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$root = Split-Path -Parent $here
$build = Join-Path $here "build"
$payload = Join-Path $build "payload"
if (Test-Path $payload) { Remove-Item -Recurse -Force $payload }
New-Item -ItemType Directory -Force (Join-Path $payload "installer") | Out-Null

# the program as the setup carries it: sources, lock file and the install script (no tests, no test data)
Copy-Item -Recurse (Join-Path $root "src") (Join-Path $payload "src")
Get-ChildItem -Recurse -Directory -Filter "__pycache__" (Join-Path $payload "src") | Remove-Item -Recurse -Force
foreach ($f in "pyproject.toml", "uv.lock", "README.md", "CHANGELOG.md") { Copy-Item (Join-Path $root $f) $payload }
foreach ($f in "install.py", "check.py", "run-install.cmd", "FociAnalyzer.ico", "FociAnalyzer.icns", "logo-512.png") { Copy-Item (Join-Path $here $f) (Join-Path $payload "installer") }

$uv = if ($env:UV_EXE) { $env:UV_EXE } else { (Get-Command uv).Source }
Copy-Item $uv (Join-Path $build "uv.exe")

$version = (Select-String -Path (Join-Path $root "pyproject.toml") -Pattern '^version\s*=\s*"([^"]+)"').Matches[0].Groups[1].Value
$iscc = @("$env:LOCALAPPDATA\Programs\Inno Setup 6\ISCC.exe", "${env:ProgramFiles(x86)}\Inno Setup 6\ISCC.exe",
          "$env:ProgramFiles\Inno Setup 6\ISCC.exe") | Where-Object { Test-Path $_ } | Select-Object -First 1
if (-not $iscc) { throw "Inno Setup 6 is not installed (ISCC.exe not found)" }
& $iscc "/DAppVersion=$version" (Join-Path $here "FociAnalyzer.iss")
if ($LASTEXITCODE -ne 0) { throw "ISCC failed" }
Get-ChildItem (Join-Path $build "FociAnalyzer-Setup-*.exe") | ForEach-Object { "Built $($_.FullName) ($([math]::Round($_.Length / 1MB, 1)) MB)" }
