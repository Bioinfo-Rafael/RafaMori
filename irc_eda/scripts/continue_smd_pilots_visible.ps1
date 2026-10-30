param(
    [Parameter(Mandatory = $true)]
    [int]$ExistingWorkerPid
)
$ErrorActionPreference = 'Stop'
$repoRoot = Split-Path (Split-Path $PSScriptRoot -Parent) -Parent
Set-Location -LiteralPath $repoRoot
$Host.UI.RawUI.WindowTitle = 'RafaMori - Gaussian SMD pilots (4 cores / 4 GB, max 2 systems)'
$pythonPath = Join-Path $repoRoot 'irc_eda\.venv\Scripts\python.exe'
$driverPath = Join-Path $PSScriptRoot 'continue_smd_pilots.py'
$logPath = Join-Path $repoRoot ('irc_eda\outputs\smd_distortion_interaction\visible_{0}.log' -f (Get-Date -Format 'yyyyMMdd_HHmmss'))
Write-Host ('PS {0}> & "{1}" -X utf8 "{2}" --wait-pid {3}' -f $repoRoot, $pythonPath, $driverPath, $ExistingWorkerPid)
Write-Host 'The existing RAFA14 worker occupies slot 1; MoriPF starts now in slot 2. Maximum two systems run at once.'
Write-Host 'Keep this window open while the continuation is running.'
& $pythonPath -X utf8 $driverPath --wait-pid $ExistingWorkerPid 2>&1 | Tee-Object -FilePath $logPath
Write-Host ('Controller exit code: {0}; log: {1}' -f $LASTEXITCODE, $logPath)
