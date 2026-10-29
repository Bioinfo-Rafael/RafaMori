# Run from any directory; Python UTF-8 mode avoids Windows legacy-codepage failures.
$ErrorActionPreference = 'Stop'
$repositoryDirectory = Split-Path (Split-Path $PSScriptRoot -Parent) -Parent
$adapterPython = Join-Path $repositoryDirectory 'irc_eda\.venv\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $adapterPython)) { throw 'Run setup_windows.ps1 first.' }
Push-Location $repositoryDirectory
try {
    & $adapterPython -X utf8 -m irc_eda @args
    exit $LASTEXITCODE
} finally { Pop-Location }
