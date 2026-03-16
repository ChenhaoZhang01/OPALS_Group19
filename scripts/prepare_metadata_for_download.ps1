param(
    [string]$MetadataPath = "metadata/master_metadata.csv",
    [string]$Output = "metadata_final.csv",
    [string]$DownloadList = "download_list.txt",
    [int]$MinReads = 200000
)

$ErrorActionPreference = "Stop"

function Convert-ReadCountToInt {
    param([string]$Value)

    if ([string]::IsNullOrWhiteSpace($Value)) {
        return $null
    }

    $normalized = $Value.Trim().ToUpperInvariant()

    if ($normalized -match '^\d+$') {
        return [int64]$normalized
    }

    if ($normalized -match '^(\d+(?:\.\d+)?)K$') {
        return [int64][Math]::Round([double]$Matches[1] * 1000)
    }

    if ($normalized -match '^(\d+(?:\.\d+)?)M$') {
        return [int64][Math]::Round([double]$Matches[1] * 1000000)
    }

    if ($normalized -match '^(\d+(?:\.\d+)?)B$') {
        return [int64][Math]::Round([double]$Matches[1] * 1000000000)
    }

    return $null
}

if (-not (Test-Path $MetadataPath)) {
    throw "Input metadata file not found: $MetadataPath"
}

$rows = Import-Csv $MetadataPath
if (-not $rows) {
    throw "Input metadata file is empty: $MetadataPath"
}

$kept = foreach ($row in $rows) {
    $readCountInt = Convert-ReadCountToInt -Value $row.read_count
    if ($null -eq $readCountInt) {
        continue
    }

    if ($readCountInt -ge $MinReads) {
        $row.read_count = [string]$readCountInt
        $row
    }
}

if (-not $kept) {
    throw "No rows passed read_count cutoff >= $MinReads"
}

$kept | Export-Csv -Path $Output -NoTypeInformation -Encoding UTF8
$kept | Select-Object -ExpandProperty sample_id | Set-Content -Path $DownloadList -Encoding UTF8

Write-Host "Input rows:" $rows.Count
Write-Host "Rows kept (read_count >= $MinReads):" $kept.Count
Write-Host "Wrote cleaned metadata:" $Output
Write-Host "Wrote download list:" $DownloadList
