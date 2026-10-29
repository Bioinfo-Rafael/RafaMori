# Run in PowerShell: powershell -ExecutionPolicy Bypass -File .\irc_eda\scripts\setup_windows.ps1
$ErrorActionPreference = 'Stop'
$adapterDirectory = Split-Path $PSScriptRoot -Parent
$repositoryDirectory = Split-Path $adapterDirectory -Parent
Set-Location $repositoryDirectory
if (-not (Get-Command git -ErrorAction SilentlyContinue)) {
    throw 'Git for Windows is required to fetch the pinned upstream implementations.'
}
& py -3 -m venv "$adapterDirectory\.venv"
if ($LASTEXITCODE -ne 0) { throw 'Python venv creation failed. Install Python 3.10 or later.' }
$adapterPython = "$adapterDirectory\.venv\Scripts\python.exe"
& $adapterPython -m pip install -r "$adapterDirectory\requirements.txt"
if ($LASTEXITCODE -ne 0) { throw 'Python dependency installation failed.' }
& $adapterPython "$PSScriptRoot\fetch_upstream.py"
if ($LASTEXITCODE -ne 0) { throw 'Pinned upstream checkout failed.' }
Write-Host 'Setup complete. Next: .\irc_eda\.venv\Scripts\python.exe -m irc_eda prepare'
