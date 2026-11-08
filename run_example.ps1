# Run the example with a local Python installation or the Codex bundled runtime.
$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $PSScriptRoot
$pythonCandidates = @(
    (Join-Path $PSScriptRoot '.venv\Scripts\python.exe'),
    (Join-Path $env:USERPROFILE '.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe')
)
$pathPython = Get-Command python -ErrorAction SilentlyContinue
if ($pathPython -and $pathPython.Source -notlike '*WindowsApps*') {
    $pythonCandidates = @($pathPython.Source) + $pythonCandidates
}
$tripPython = $pythonCandidates | Where-Object { Test-Path -LiteralPath $_ } | Select-Object -First 1
if (-not $tripPython) {
    throw 'Install Python 3.10 or later and run: python -m pip install -r requirements.txt'
}
$runtimePackages = Join-Path $PSScriptRoot '.runtime'
if (Test-Path -LiteralPath $runtimePackages) {
    if ($env:PYTHONPATH) { $env:PYTHONPATH = $runtimePackages + [IO.Path]::PathSeparator + $env:PYTHONPATH }
    else { $env:PYTHONPATH = $runtimePackages }
}
& $tripPython trip_design.py --config examples/m6_synthetic/config.json --output results/m6_cp
if ($LASTEXITCODE -ne 0) { throw 'Cp example failed; inspect its report.' }
& $tripPython trip_design.py --config examples/m6_synthetic/config_profiles.json --output results/m6_profiles
if ($LASTEXITCODE -ne 0) { throw 'Profile example failed; inspect its report.' }
Write-Host 'Done. Open results\m6_cp\report.html or results\m6_profiles\report.html.'
