#!/bin/bash
# save as push_epiclock.sh
# Push /home/n.atallah/epiclock/ec01/final_code to GitHub Epiclock repo

PROJECT_DIR="/home/n.atallah/epiclock/ec01/final_code"
REPO_SSH="git@github.com:Nabilatallah/Epiclock.git"

echo "📦 Navigating to project directory..."
cd "$PROJECT_DIR" || { echo "❌ Directory not found: $PROJECT_DIR"; exit 1; }

# Initialize repo if needed
if [ ! -d ".git" ]; then
    echo "🆕 Initializing Git repository..."
    git init
    git remote add origin "$REPO_SSH"
else
    echo "✅ Git repository already initialized."
fi

# Ensure remote is correct
git remote set-url origin "$REPO_SSH"

# === Auto-stage & commit before pull ===
if ! git diff --quiet || ! git diff --cached --quiet; then
    echo "💾 Auto-committing local changes before pull..."
    git add .
    git commit -m "Auto-commit before rebase on $(date '+%Y-%m-%d %H:%M:%S')"
else
    echo "✅ No local changes to commit before pull."
fi

# Fetch latest from GitHub & rebase
echo "🔄 Pulling latest changes from GitHub..."
git fetch origin main
git pull --rebase origin main || { echo "⚠️ Pull failed — check for merge conflicts."; exit 1; }

# Stage changes again after rebase
git add .

# Commit changes after pull (if any new ones)
git commit -m "Update from final_code on $(date '+%Y-%m-%d %H:%M:%S')" || echo "⚠️ No new changes to commit after pull."

# Push to GitHub
echo "🚀 Pushing to GitHub..."
git branch -M main
git push -u origin main

echo "✅ Done! Changes pushed to: $REPO_SSH"
