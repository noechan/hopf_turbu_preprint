# GitHub publication workflow

## Before making the repository public

1. Select an institutionally approved licence for the authors' original code
   and add it as `LICENSE`. The bundled GPL-3.0 `permutationTest.m` files and
   other attributed MATLAB helpers retain their own terms; review
   `THIRD_PARTY_NOTICES.md` before choosing a repository-wide licence.
2. Confirm that redistribution of the derived Neurosynth maps and Schaefer
   atlas is permitted and add the final resource citations to the manuscript
   or README.
3. Run the privacy and staging checks below. Keep the remote private until the
   staged file list has been inspected.

`CITATION.cff`, `.github/`, and `config/` are deliberately not required for
this publication release.

## Initialize the local history

From this repository directory:

```bash
git init -b main
git add .gitignore .gitattributes
git commit -m "chore: add repository safeguards"

git add README.md DATA.md REPRODUCIBILITY.md GITHUB_UPLOAD.md
git add THIRD_PARTY_NOTICES.md turbulence hopf machine_learning neuromaps_analysis
git status --short
git diff --cached --stat
git commit -m "release: add manuscript analysis code"
```

Before the second commit, run:

```bash
git grep -n -E '/Users/|/Volumes/|[0-9]{3}_S_[0-9]{4}' --cached || true
git diff --cached --name-only
git diff --cached --stat
```

## Create the private GitHub remote

Replace `<github-user>` with the account or organisation that will own the
repository:

```bash
gh auth status
gh repo create <github-user>/hopf_turbu_ms \
  --private --source=. --remote=origin --push
```

Review the complete staged file list before each commit. Do not use `git add .`
until the identifier, private-path, and file-size checks have passed.

After verifying the private remote, change its visibility to public in the
GitHub repository settings or run:

```bash
gh repo edit <github-user>/hopf_turbu_ms --visibility public
```

## Submission release

```bash
git tag -a v1.0.0-submission.1 \
  -m "Code corresponding to the manuscript submission"
git push origin v1.0.0-submission.1
gh release create v1.0.0-submission.1 \
  --title "Manuscript submission release" \
  --notes "Code corresponding to the submitted manuscript."
```

Reference the tag URL, not `tree/main`, in the manuscript.

## Peer-review updates

For each reviewer request:

```bash
git switch main
git pull --ff-only
git switch -c review/r1-short-description

# edit and validate
git add path/to/changed/files
git commit -m "analysis: address reviewer request"
git push -u origin review/r1-short-description
gh pr create --fill
```

Merge reviewed changes into `main`, update `REPRODUCIBILITY.md` when an
analysis changes, and create a new immutable release such as
`v1.0.0-review.1`. Never move or overwrite an existing manuscript tag.
