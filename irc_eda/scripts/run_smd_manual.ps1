param(
    [Parameter(Mandatory = $true)]
    [ValidateSet('RAFA14','MoriPF','RAFA15','RAFA16')]
    [string]$System
)
$ErrorActionPreference = 'Stop'
$repoRoot = Split-Path (Split-Path $PSScriptRoot -Parent) -Parent
Set-Location -LiteralPath $repoRoot
$pythonPath = Join-Path $repoRoot 'irc_eda\.venv\Scripts\python.exe'
$scriptPath = Join-Path $PSScriptRoot 'smd_distortion_pilot.py'
Write-Host "Running $System pilot only. New inputs: 4 cores / 4 GB EACH. Keep at most two such windows running."
& $pythonPath -X utf8 $scriptPath run $System
if ($LASTEXITCODE -ne 0) { throw "$System calculation failed; preserved output. No further jobs launched." }
& $pythonPath -X utf8 $scriptPath analyze $System
if ($LASTEXITCODE -ne 0) { throw "$System analysis failed." }
