param(
    [string]$Output = "metadata/master_metadata.csv",
    [int]$PerQueryLimit = 5000,
    [double]$SleepSeconds = 0.1
)

$ErrorActionPreference = "Stop"

$RepoRoot = Split-Path -Parent $PSScriptRoot
Push-Location $RepoRoot

try {
    if (-not (Test-Path "metadata")) {
        New-Item -ItemType Directory -Path "metadata" | Out-Null
    }

    if (-not (Test-Path ".venv")) {
        if (Get-Command py -ErrorAction SilentlyContinue) {
            py -3 -m venv .venv
        }
        elseif (Get-Command python -ErrorAction SilentlyContinue) {
            python -m venv .venv
        }
        else {
            Write-Error "Python was not found. Install Python 3.10+ and rerun this script."
        }
    }

    $VenvPython = Join-Path $RepoRoot ".venv\Scripts\python.exe"
    if (-not (Test-Path $VenvPython)) {
        Write-Error "Virtual environment python executable not found at $VenvPython"
    }

    $builderArgs = @(
        "scripts/build_master_metadata.py",
        "--output", $Output,
        "--per-query-limit", $PerQueryLimit,
        "--sleep", $SleepSeconds
    )

    & $VenvPython @builderArgs

    if ($LASTEXITCODE -ne 0) {
        Write-Error "Metadata build failed with exit code $LASTEXITCODE"
    }

    Write-Host "Phase 1 complete: $Output"
}
finally {
    Pop-Location
}
