param(
    [Parameter(Mandatory = $true)]
    [string]$ProjectName,
    [string]$MetadataPath = "metadata/master_metadata.csv",
    [int]$MinReads = 200000,
    [switch]$CopyReferenceMatrix
)

$ErrorActionPreference = "Stop"

$RepoRoot = Split-Path -Parent $PSScriptRoot
Push-Location $RepoRoot

try {
    $ProjectRoot = Join-Path $RepoRoot ("projects/{0}" -f $ProjectName)

    $Dirs = @(
        "metadata",
        "raw_reads",
        "assemblies",
        "proteins",
        "results",
        "results/arg_hits",
        "analysis",
        "analysis/figures"
    )

    foreach ($Dir in $Dirs) {
        $FullPath = Join-Path $ProjectRoot $Dir
        if (-not (Test-Path $FullPath)) {
            New-Item -ItemType Directory -Path $FullPath -Force | Out-Null
        }
    }

    $PrepareScript = Join-Path $PSScriptRoot "prepare_metadata_for_download.ps1"
    if (-not (Test-Path $PrepareScript)) {
        throw "Required script not found: $PrepareScript"
    }

    $CleanMetadataPath = Join-Path $ProjectRoot "metadata/metadata_final.csv"
    $DownloadListPath = Join-Path $ProjectRoot "metadata/download_list.txt"

    & $PrepareScript `
        -MetadataPath $MetadataPath `
        -Output $CleanMetadataPath `
        -DownloadList $DownloadListPath `
        -MinReads $MinReads

    if ($LASTEXITCODE -ne 0) {
        throw "Metadata preparation failed with exit code $LASTEXITCODE"
    }

    if ($CopyReferenceMatrix) {
        $ReferenceMatrix = Join-Path $RepoRoot "results/ARG_matrix.csv"
        if (Test-Path $ReferenceMatrix) {
            Copy-Item -Path $ReferenceMatrix -Destination (Join-Path $ProjectRoot "results/ARG_matrix_reference.csv") -Force
        }
    }

    $ProjectReadmePath = Join-Path $ProjectRoot "README.md"
    $ProjectReadme = @"
# $ProjectName Workflow

This folder isolates one research project.

## What is already created

- metadata/metadata_final.csv (read_count numeric, cutoff >= $MinReads)
- metadata/download_list.txt (accession IDs only)
- raw_reads/
- assemblies/
- proteins/
- results/arg_hits/
- analysis/figures/

## Start here (intern checklist)

1. Download reads from metadata/download_list.txt into raw_reads/ using prefetch and fastq-dump.
2. Run FastQC on raw_reads/*.fastq and review quality metrics.
3. Assemble each sample with MEGAHIT into assemblies/<sample_id>/.
4. Run Prodigal on each final.contigs.fa and write proteins/<sample_id>.faa.
5. Run DIAMOND against CARD and write results/arg_hits/<sample_id>.tsv.
6. Build ARG matrix with scripts/build_arg_matrix.py from results/arg_hits/*.tsv.
7. Normalize and merge with metadata using scripts/build_arg_dataset.py.
8. Generate summary figures with scripts/plot_first_summary.py.
9. Validate quality with scripts/check_dataset_quality.py.
"@

    Set-Content -Path $ProjectReadmePath -Value $ProjectReadme -Encoding UTF8

    Write-Host "Project initialized:" $ProjectRoot
    Write-Host "Metadata file:" $CleanMetadataPath
    Write-Host "Download list:" $DownloadListPath
}
finally {
    Pop-Location
}
