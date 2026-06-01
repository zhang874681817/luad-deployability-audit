param(
    [Parameter(Mandatory=$false)]
    [string]$RepoUrl,

    [Parameter(Mandatory=$false)]
    [string]$CommitMessage = "Release LUAD deployability audit package",

    [Parameter(Mandatory=$false)]
    [string]$Branch = "main"
)

$ErrorActionPreference = "Stop"

function Require-Command($Name) {
    if (-not (Get-Command $Name -ErrorAction SilentlyContinue)) {
        throw "Required command '$Name' was not found. Install it or add it to PATH."
    }
}

Require-Command git

$repoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $repoRoot

if (-not $RepoUrl) {
    $RepoUrl = Read-Host "Paste the HTTPS GitHub repository URL, for example https://github.com/USERNAME/REPO.git"
}

if (-not $RepoUrl) {
    throw "Repository URL is required."
}

$name = git config user.name 2>$null
$email = git config user.email 2>$null

if (-not $name) {
    $name = Read-Host "Git author name for this repository"
    git config user.name $name
}

if (-not $email) {
    $email = Read-Host "Git author email for this repository"
    git config user.email $email
}

if (-not (Test-Path ".git")) {
    git init
}

git branch -M $Branch

$remoteExists = (git remote 2>$null) -contains "origin"

if ($remoteExists) {
    git remote set-url origin $RepoUrl
} else {
    git remote add origin $RepoUrl
}

git add .

$status = git status --porcelain
if ($status) {
    git commit -m $CommitMessage
} else {
    Write-Host "No file changes to commit."
}

$commit = git rev-parse HEAD
if ($LASTEXITCODE -ne 0 -or -not $commit) {
    throw "No commit was created. Check git status and author identity."
}
Write-Host "Frozen commit hash: $commit"

git push -u origin $Branch
if ($LASTEXITCODE -ne 0) {
    throw "git push failed. Check GitHub credentials and repository URL."
}

Write-Host ""
Write-Host "Upload complete."
Write-Host "Repository: $RepoUrl"
Write-Host "Commit: $commit"
