# restart-run-code.ps1 — idempotent, safe, verbose
param(
  [int]$Port = 3001,
  [string]$Script = "microservices\run-code\server-e2b-working.js"
)

Write-Host "== Preflight: env & proxies =="
if (-not $env:E2B_API_KEY)     { throw "E2B_API_KEY is not set" }
if (-not $env:E2B_TEMPLATE_ID) { throw "E2B_TEMPLATE_ID is not set" }

# Ensure local requests aren't sent through a proxy
$env:NO_PROXY        = "127.0.0.1,localhost"
Remove-Item Env:HTTP_PROXY   -ErrorAction SilentlyContinue
Remove-Item Env:HTTPS_PROXY  -ErrorAction SilentlyContinue
Remove-Item Env:http_proxy   -ErrorAction SilentlyContinue
Remove-Item Env:https_proxy  -ErrorAction SilentlyContinue

Write-Host "== Kill anything on :$Port =="
$lines = (netstat -ano | Select-String ":$Port\s+LISTENING")
if ($lines) {
  $pids = $lines -replace '.*\s+(\d+)$','$1' | Select-Object -Unique
  foreach ($pid in $pids) {
    try { taskkill /PID $pid /F | Out-Null } catch {}
  }
  Start-Sleep -Milliseconds 300
}

Write-Host "== Start server on 127.0.0.1:$Port =="
$env:HOST = "127.0.0.1"
$env:PORT = "$Port"
# Start Node detached so this script can continue
Start-Process -FilePath "node" -ArgumentList "$Script" -WindowStyle Hidden

Write-Host "== Wait for /health (30s) =="
$ok = $false
for ($i=0; $i -lt 30; $i++) {
  try {
    $resp = Invoke-WebRequest -UseBasicParsing "http://127.0.0.1:$Port/health" -TimeoutSec 2
    if ($resp.StatusCode -eq 200 -and $resp.Content -match "OK|healthy") { $ok = $true; break }
  } catch {}
  Start-Sleep -Seconds 1
}
if (-not $ok) { throw "Service failed health check on http://127.0.0.1:$Port/health" }

Write-Host "✅ run-code is up at http://127.0.0.1:$Port"