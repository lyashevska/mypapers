DEVNOTES
=======

This file contains short developer notes for setting up a local environment and previewing the generated publications block.

Local setup (recommended):

1. Create and activate a virtual environment, then install pinned dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

2. Run the generator and verifier to preview the publications block locally:

```bash
python bib_to_md.py
python verify_publications.py
```

What the scripts do:

- `bib_to_md.py` reads all `.bib` / `.bibtex` files from the `bib/` directory, normalizes author names, gathers available metadata (year, volume, number, pages) and writes `publications.md` with a title line `## 📚 Publications` (the GitHub Actions workflow expects this line and uses `tail -n +2` when injecting into the README).

- `scripts/verify_readme_bibs.py` checks that every BibTeX key in `bib/` appears in the publications block in `README.md` (between `<!-- PUBLICATIONS START -->` and `<!-- PUBLICATIONS END -->`).

Notes and conventions:

- Keep BibTeX keys consistent between the `bib/` and `pdf/` directories. Example: for a publication with key `lyashevska2016a`, provide `bib/lyashevska2016a.bib` and `pdf/lyashevska2016a.pdf` (the generator will only emit a [PDF] link when the corresponding `pdf/{key}.pdf` file exists).

- `publications.md` is an intermediate generated file and does not need to be tracked in Git; CI generates it during the workflow.

- If you rename BibTeX keys or PDF filenames, re-run the generator and `scripts/verify_readme_bibs.py` locally to ensure the README block remains consistent.

- The CI workflow `.github/workflows/update-readme.yml` installs dependencies, runs the generator, verifies keys, replaces the publications block in `README.md`, and commits the updated README back to the repository.

Contact / notes:
- If you encounter mismatched or ambiguous PDF ↔ BibTeX mappings, prefer manual resolution (rename either the bib key or the PDF) rather than broad automatic renames.
