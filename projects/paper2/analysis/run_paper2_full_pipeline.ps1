param(
    [string]$PythonExe = "python",
    [switch]$UseDemoBaselineData,
    [switch]$RunLowIdentity,
    [switch]$AllowSyntheticData,
    [string]$QueryFasta,
    [string]$DbFasta,
    [string]$QueryLabels,
    [string]$DbLabels,
    [bool]$FailOnSyntheticNonDemo = $true
)

$ErrorActionPreference = "Stop"

$repo = (Resolve-Path (Join-Path $PSScriptRoot "..\..\..")).Path
$paper2 = Join-Path $repo "projects\paper2"
$results = Join-Path $paper2 "results"
$analysis = Join-Path $paper2 "analysis"
$figures = Join-Path $analysis "figures"

$embeddings = Join-Path $results "embeddings\protein_embeddings.npy"
$labels = Join-Path $results "training_labels_template.csv"
$metrics = Join-Path $results "blast_vs_ml_metrics.csv"

Write-Host "[1/3] Running embedding model training"
& $PythonExe (Join-Path $analysis "train_arg_classifier.py") --embeddings $embeddings --labels $labels --metrics-out $metrics

Write-Host "[2/3] Running BLAST baseline (if required files exist)"
$blastRoot = Join-Path $repo "tools\blast\ncbi-blast-2.17.0+\bin"
$blastp = Join-Path $blastRoot "blastp.exe"
$makeblastdb = Join-Path $blastRoot "makeblastdb.exe"

$queryFasta = if ($QueryFasta) { $QueryFasta } else { Join-Path $paper2 "proteins\query_proteins.faa" }
$dbFasta = if ($DbFasta) { $DbFasta } else { Join-Path $paper2 "proteins\card_reference.faa" }
$queryLabels = if ($QueryLabels) { $QueryLabels } else { Join-Path $results "query_labels.csv" }
$dbLabels = if ($DbLabels) { $DbLabels } else { Join-Path $results "db_labels.csv" }
$blastMetricsOut = $metrics

function Get-FastaSequenceCount {
    param([Parameter(Mandatory = $true)][string]$Path)
    if (-not (Test-Path $Path)) { return 0 }
    return (Select-String -Path $Path -Pattern '^>' -AllMatches | Measure-Object).Count
}

if ($UseDemoBaselineData) {
    $demoDir = Join-Path $results "demo_baseline"
    New-Item -ItemType Directory -Force -Path $demoDir | Out-Null
    $blastMetricsOut = Join-Path $demoDir "blast_vs_ml_metrics_demo.csv"
}

$queryCount = Get-FastaSequenceCount -Path $queryFasta
$dbCount = Get-FastaSequenceCount -Path $dbFasta
$queryIsPlaceholder = ([System.IO.Path]::GetFileName($queryFasta) -eq "PLACE_PROTEIN_FASTA_HERE.faa")
$dbIsPlaceholder = ([System.IO.Path]::GetFileName($dbFasta) -eq "PLACE_PROTEIN_FASTA_HERE.faa")
$usesPlaceholderInput = $queryIsPlaceholder -or $dbIsPlaceholder
$likelySyntheticInputs = ($queryCount -lt 20) -or ($dbCount -lt 20) -or $usesPlaceholderInput
$blockedBySyntheticGuard = $likelySyntheticInputs -and -not $UseDemoBaselineData -and -not $AllowSyntheticData

if ($blockedBySyntheticGuard) {
    Write-Host "BLAST baseline and low-identity steps skipped in non-demo mode: synthetic/toy inputs detected."
    Write-Host "  Query sequences: $queryCount"
    Write-Host "  DB sequences:    $dbCount"
    Write-Host "  Placeholder selected as input: $usesPlaceholderInput"
    Write-Host "Use -UseDemoBaselineData for demo artifacts, or -AllowSyntheticData to force non-demo execution."
    Write-Host "You can also provide real inputs using -QueryFasta, -DbFasta, -QueryLabels, and -DbLabels."

    if ($FailOnSyntheticNonDemo) {
        throw "Non-demo run blocked: synthetic/toy BLAST inputs detected."
    }
}

if ((-not $blockedBySyntheticGuard) -and (Test-Path $blastp) -and (Test-Path $makeblastdb) -and (Test-Path $queryFasta) -and (Test-Path $dbFasta) -and (Test-Path $queryLabels) -and (Test-Path $dbLabels)) {
    & $PythonExe (Join-Path $analysis "run_blast_baseline.py") `
        --query-fasta $queryFasta `
        --db-fasta $dbFasta `
        --query-labels $queryLabels `
        --db-labels $dbLabels `
    --metrics-out $blastMetricsOut `
        --blast-bin $blastp `
        --makeblastdb-bin $makeblastdb
} else {
    if ($blockedBySyntheticGuard) {
        Write-Host "BLAST baseline skipped by synthetic-input guard in non-demo mode."
    } else {
        Write-Host "BLAST baseline skipped: add these files first:"
        Write-Host "  $queryFasta"
        Write-Host "  $dbFasta"
        Write-Host "  $queryLabels"
        Write-Host "  $dbLabels"
    }
}

Write-Host "[3/3] Generating Paper 2 figures"
& $PythonExe (Join-Path $analysis "generate_paper2_figures.py")

if ($RunLowIdentity) {
    Write-Host "[4/4] Running low-identity comparison"
    $lowIdentityRan = $false
    $lowOut = Join-Path $results "low_identity_comparison.csv"
    $identityBinOut = Join-Path $results "identity_bin_recall.csv"
    $perQueryOut = Join-Path $results "low_identity_per_query.csv"

    $rfRow = Import-Csv $metrics | Where-Object { $_.method -eq "RandomForest_embeddings" } | Select-Object -First 1
    if (-not $rfRow) {
        throw "Could not find RandomForest_embeddings row in $metrics"
    }

    if ($UseDemoBaselineData) {
        $demoDir = Join-Path $results "demo_baseline"
        New-Item -ItemType Directory -Force -Path $demoDir | Out-Null
        $lowOut = Join-Path $demoDir "low_identity_comparison.csv"
        $identityBinOut = Join-Path $demoDir "identity_bin_recall.csv"
        $perQueryOut = Join-Path $demoDir "low_identity_per_query.csv"
    }

    if ((-not $blockedBySyntheticGuard) -and (Test-Path $blastp) -and (Test-Path $makeblastdb) -and (Test-Path $queryFasta) -and (Test-Path $dbFasta) -and (Test-Path $queryLabels) -and (Test-Path $dbLabels)) {
        & $PythonExe (Join-Path $analysis "run_low_identity_experiment.py") `
            --query-fasta $queryFasta `
            --db-fasta $dbFasta `
            --query-labels $queryLabels `
            --db-labels $dbLabels `
            --out-csv $lowOut `
            --identity-bin-out $identityBinOut `
            --per-query-out $perQueryOut `
            --blast-bin $blastp `
            --makeblastdb-bin $makeblastdb `
            --low-identity-threshold 40 `
            --embedding-precision ([double]$rfRow.precision) `
            --embedding-recall ([double]$rfRow.recall) `
            --embedding-f1 ([double]$rfRow.f1)
        $lowIdentityRan = $true
    } else {
        if ($blockedBySyntheticGuard) {
            Write-Host "Low-identity run skipped by synthetic-input guard in non-demo mode."
        } else {
            Write-Host "Low-identity run skipped: required BLAST/query/db files not found"
        }
    }
}

Write-Host "Pipeline complete. Metrics: $metrics"
if ($UseDemoBaselineData) {
    Write-Host "Demo BLAST metrics: $blastMetricsOut"
}
if ($RunLowIdentity) {
    if ($lowIdentityRan) {
        if ($UseDemoBaselineData) {
            Write-Host "Demo low-identity metrics: $(Join-Path $results 'demo_baseline\low_identity_comparison.csv')"
        } else {
            Write-Host "Low-identity metrics: $(Join-Path $results 'low_identity_comparison.csv')"
        }
    }
}
Write-Host "Figures dir: $figures"
