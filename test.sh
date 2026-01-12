#!/bin/bash
# Local test script to verify publications workflow
# Run this before pushing changes to ensure CI will pass

set -e

echo "🔍 Testing publications workflow..."
echo

echo "1️⃣  Checking Python environment..."
if ! command -v python &> /dev/null; then
    echo "❌ Python not found. Please install Python 3."
    exit 1
fi
echo "✅ Python found: $(python --version)"
echo

echo "2️⃣  Installing dependencies..."
pip install -q -r requirements.txt
echo "✅ Dependencies installed"
echo

echo "3️⃣  Running bib_to_md.py..."
python bib_to_md.py
echo "✅ publications.md generated"
echo

echo "4️⃣  Updating README..."
awk '/<!-- PUBLICATIONS START -->/,/<!-- PUBLICATIONS END -->/ {next} 1' README.md > tmp.md
echo '<!-- PUBLICATIONS START -->' >> tmp.md
echo '## 📚 Publications' >> tmp.md
tail -n +2 publications.md >> tmp.md
echo '<!-- PUBLICATIONS END -->' >> tmp.md
mv tmp.md README.md
echo "✅ README updated"
echo

echo "5️⃣  Verifying all bib keys are in README..."
python scripts/verify_readme_bibs.py
echo

echo "✅ All checks passed! Safe to commit and push."
