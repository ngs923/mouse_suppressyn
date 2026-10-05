#!/usr/bin/env bash
# Creates the first commit and configures the remote.
# It stops short of pushing, since that needs your credentials.
#
#   bash setup_repo.sh
#
# Then run the git push command it prints.

set -euo pipefail

REMOTE="https://github.com/ngs923/mouse_suppressyn.git"

if [ -d .git ]; then
  echo "Already a git repository; skipping init."
else
  git init
fi

# Use main as the default branch, matching GitHub's default
git symbolic-ref HEAD refs/heads/main 2>/dev/null || git branch -M main

git add .
echo
echo "--- files to be committed ---"
git status --short
echo

git commit -m "Add ERV env RNA-seq pipeline and suppressyn dot plot analysis"

if git remote get-url origin >/dev/null 2>&1; then
  git remote set-url origin "$REMOTE"
else
  git remote add origin "$REMOTE"
fi

cat <<'MSG'

The local repository is ready. Push it yourself with:

    git push -u origin main

When prompted for a password, use a Personal Access Token rather than your
account password. The GitHub CLI makes this easier:

    brew install gh
    gh auth login
    git push -u origin main
MSG
