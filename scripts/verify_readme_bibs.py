#!/usr/bin/env python3
"""Verify that every bib key under bib/ appears in README.md's publications block.

Exit codes:
 0 - OK
 1 - missing README or markers
 2 - missing bib keys found
"""
from pathlib import Path
import re
import sys

repo_root = Path(__file__).resolve().parent.parent
bib_dir = repo_root / 'bib'
readme = repo_root / 'README.md'

if not readme.exists():
    print('ERROR: README.md not found')
    sys.exit(1)

text = readme.read_text(encoding='utf-8')
start = text.find('<!-- PUBLICATIONS START -->')
end = text.find('<!-- PUBLICATIONS END -->')
if start == -1 or end == -1:
    print('ERROR: publications markers not found in README.md')
    sys.exit(1)

block = text[start:end]

# Collect bib keys from bib/ files
bib_keys = []
for p in sorted(bib_dir.glob('*.bib')) + sorted(bib_dir.glob('*.bibtex')):
    data = p.read_text(encoding='utf-8')
    # find keys like @article{key,
    m = re.findall(r"@\w+\{([^,\n]+)", data)
    for k in m:
        bib_keys.append(k.strip())

missing = []
for k in sorted(set(bib_keys)):
    # check if key appears literally in the publications block
    if k not in block:
        missing.append(k)

if missing:
    print('Missing bib keys in README publications block:')
    for k in missing:
        print('-', k)
    sys.exit(2)

print('OK: all bib keys referenced in README publications block')
sys.exit(0)
