# One-Click GitHub Upload

## 1. Create an empty GitHub repository

Create a new empty GitHub repository in your account. Do not add a README, license or `.gitignore` on GitHub because this folder already contains them.

Copy the HTTPS repository URL, for example:

```text
https://github.com/USERNAME/luad-deployability-audit.git
```

## 2. Run the upload script

Open PowerShell in this folder and run:

```powershell
.\upload_to_github.ps1 -RepoUrl "https://github.com/USERNAME/luad-deployability-audit.git"
```

If Git asks for credentials, sign in with your GitHub username and a GitHub personal access token, or use Git Credential Manager when prompted.

## 3. Copy the frozen commit hash

The script prints:

```text
Frozen commit hash: <hash>
```

Use that hash in the manuscript Code Availability statement.

## 4. Archive for DOI

After the GitHub upload is correct, archive the GitHub release with Zenodo or OSF and add the DOI to the manuscript.
