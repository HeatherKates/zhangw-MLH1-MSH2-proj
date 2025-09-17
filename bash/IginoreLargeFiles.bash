#!/bin/bash
# Find and ignore files over 100MB
find . -type f -size +100M -exec echo "{}" \; >> .gitignore
# Remove duplicates
sort -u -o .gitignore .gitignore
# Optional: Add and commit changes
#git add .gitignore
#git commit -m "Automatically ignore large files"

