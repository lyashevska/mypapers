import os
import glob
import bibtexparser


def format_entry(entry):
    """Format one bib entry nicely in Markdown."""
    authors = entry.get('author', '').replace('\n', ' ').replace(' and ', ', ')
    title = entry.get('title', '')
    # strip surrounding braces commonly used in bibtex titles
    title = title.strip().strip('{}')
    year = entry.get('year', '')
    journal = entry.get('journal') or entry.get('booktitle', '')
    key = entry.get('ID', '')

    # Remove line breaks and extra spaces
    authors = ' '.join(authors.split())

    # Build file paths: repository enforces bib/ and pdf/ directories
    bib_candidates = [f"bib/{key}.bib", f"bib/{key}.bibtex"]
    pdf_candidates = [f"pdf/{key}.pdf"]

    bib_link = ''
    for p in bib_candidates:
        if os.path.exists(p):
            bib_link = f"[BibTeX]({p})"
            break

    pdf_link = ''
    for p in pdf_candidates:
        if os.path.exists(p):
            pdf_link = f"[PDF]({p})"
            break

    links = " | ".join([link for link in [pdf_link, bib_link] if link])

    # Format Markdown line
    parts = [f"- **{authors}**"]
    if year:
        parts.append(f"({year}).")
    parts.append(f"*{title}.*")
    if journal:
        parts.append(f"_{journal}._")
    if links:
        parts.append(links)

    return ' '.join(parts).strip()


def main():
    # Collect .bib and .bibtex files from bib/ directory only
    bib_files = []
    if os.path.isdir('bib'):
        bib_files.extend(glob.glob('bib/*.bib'))
        bib_files.extend(glob.glob('bib/*.bibtex'))

    # remove duplicates and keep stable order
    bib_files = list(dict.fromkeys(bib_files))

    all_entries = []

    for bibfile in bib_files:
        try:
            with open(bibfile, encoding='utf-8') as bibtex_file:
                bib_database = bibtexparser.load(bibtex_file)
                all_entries.extend(bib_database.entries)
        except Exception as e:
            print(f"Warning: failed to read {bibfile}: {e}")

    # Sort newest first (missing years come last)
    all_entries.sort(key=lambda e: e.get('year', ''), reverse=True)

    lines = ['## 📚 Publications']
    for entry in all_entries:
        lines.append(format_entry(entry))

    with open('publications.md', 'w', encoding='utf-8') as out:
        out.write('\n'.join(lines))


if __name__ == '__main__':
    main()
