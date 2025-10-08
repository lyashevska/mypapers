## Quick orientation for AI coding agents

This repo is a small, single-author publications workspace: a Python script (`bib_to_md.py`) converts BibTeX files into `publications.md` and a GitHub Actions workflow (`.github/workflows/update-readme.yml`) injects that generated content into `README.md` between the HTML markers `<!-- PUBLICATIONS START -->` and `<!-- PUBLICATIONS END -->`.

Keep edits minimal and predictable: change `bib_to_md.py` only when you need to adjust parsing/formatting or file-path handling. The action runs `pip install bibtexparser` and `python bib_to_md.py` on pushes to `**.bib` and `bib_to_md.py`.

Key files to read and update together
- `bib_to_md.py` — core generator. Sorts entries by `year` (descending), formats authors/titles and emits `publications.md` with a title line `## 📚 Publications`.
- `.github/workflows/update-readme.yml` — CI automation that runs the script and replaces the README block. It expects the generated `publications.md` to have a title line which the workflow skips with `tail -n +2`.
- `README.md` — contains the publications block delimited by `<!-- PUBLICATIONS START -->` / `<!-- PUBLICATIONS END -->`.
- `publications.md` — intermediate artifact created by `bib_to_md.py`.

Repository layout requirement
- All BibTeX files must live in the `bib/` directory (accepted extensions: `.bib`, `.bibtex`).
- All PDFs must live in the `pdf/` directory. The generator now reads only `bib/` and builds PDF links as `pdf/{key}.pdf`.

Migration tip
- To move existing root-level files into the required layout (run from repo root):

```bash
mkdir -p bib pdf
mv *.bib *.bibtex bib/ 2>/dev/null || true
mv *.pdf pdf/ 2>/dev/null || true
```

After moving files, run the generator and verify:

```bash
pip install bibtexparser
python bib_to_md.py
python verify_publications.py
```

Repository-specific conventions and gotchas
- The workflow replaces only the content between the two HTML comments; preserve those markers exactly.
- `bib_to_md.py` currently builds links as `bib/{key}.bib` and `pdf/{key}.pdf` and checks `os.path.exists` for those paths. In this repository many `.bib` and `.pdf` files are at the repo root (not under `bib/` or `pdf/`). Before changing file-path logic, verify whether files are moved into `bib/`/`pdf/` on purpose or update the script to glob for `*.bib` at the repo root.
- The script uses `bibtexparser` — CI and the workflow install that. Locally install the same package to reproduce results.
- Formatting behavior to preserve when changing the generator: authors converted from `and` to `, `; title braces `{}` are stripped; missing `year` sorts as empty string (older entries may shift position).

How to run & test locally
1) Install dependency and run generator:

```bash
pip install bibtexparser
python bib_to_md.py
```

2) Mirror the workflow's README injection (the workflow runs this automatically; run locally to preview):

```bash
awk '/<!-- PUBLICATIONS START -->/,/<!-- PUBLICATIONS END -->/ {next} 1' README.md > tmp.md
echo '<!-- PUBLICATIONS START -->' >> tmp.md
echo '## 📚 Publications' >> tmp.md
tail -n +2 publications.md >> tmp.md
echo '<!-- PUBLICATIONS END -->' >> tmp.md
mv tmp.md README.md
git add README.md
git commit -m "Update publications" || echo "No changes to commit"
git push
```

Patterns worth preserving (examples)
- Sorting: `all_entries.sort(key=lambda e: e.get("year", ""), reverse=True)` — keep descending-year behavior.
- Link building: the code checks `os.path.exists` and only emits `[PDF]` / `[BibTeX]` links when the corresponding files exist.
- Publications block lifecycle: generator -> `publications.md` -> workflow inserts into `README.md`.

When you modify code
- Run the generator locally and confirm `publications.md` looks correct and that the workflow injection steps (above) produce the intended `README.md` output.
- If you change dependencies, update the workflow to `pip install` the new packages.
- If you change the `publications.md` title line, also update the workflow's `tail -n +2` logic.

If anything in this file is unclear or you want more examples (e.g., patch showing how to make the generator accept top-level `.bib` files), tell me which area to expand and I will update this file.
