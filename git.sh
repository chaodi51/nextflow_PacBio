#!/usr/bin/bash
echo "pushing changes to Github..."
change=$1

git add .
git commit -m $change
git push origin main