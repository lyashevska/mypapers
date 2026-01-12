import os
import glob
import bibtexparser


def format_entry(entry):
    """Format one bib entry nicely in Markdown."""
    raw_authors = entry.get('author', '').replace('\n', ' ')
    # Parse authors into a list and normalize to 'Lastname, Firstname' where possible
    def normalize_name(name: str) -> str:
        name = name.strip()
        if not name:
            return ''
        # If name already contains a comma, assume 'Last, First' format
        if ',' in name:
            parts = [p.strip() for p in name.split(',', 1)]
            last = parts[0]
            first = parts[1]
            return f"{last}, {first}" if first else last
        # Otherwise, assume 'First Middle Last' -> 'Last, First Middle'
        parts = name.split()
        if len(parts) == 1:
            return parts[0]
        last = parts[-1]
        first = ' '.join(parts[:-1])
        return f"{last}, {first}"
    authors_list = [normalize_name(a) for a in raw_authors.split(' and ') if a.strip()]
    # separate multiple authors with a semicolon for clarity
    authors = '; '.join(authors_list)
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

    # Journal + optional volume/number/pages
    if journal:
        vol = entry.get('volume', '').strip()
        num = entry.get('number', '').strip()
        pages = entry.get('pages', '').strip()

        journal_part = f"_{journal}._"
        volparts = ''
        if vol:
            volparts += vol
            if num:
                volparts += f"({num})"
        if pages:
            if volparts:
                volparts += f": {pages}"
            else:
                volparts = f"{pages}"

        if volparts:
            journal_part = f"{journal_part} {volparts}."

        parts.append(journal_part)

    if links:
        parts.append(links)

    return ' '.join(parts).strip()


def main():
    # Collect .bib and .bibtex files from bib/ directory only
    bib_files = []
    if not os.path.isdir('bib'):
        print("ERROR: bib/ directory not found")
        print("Please create a bib/ directory and add .bib or .bibtex files to it.")
        return 1
    
    bib_files.extend(glob.glob('bib/*.bib'))
    bib_files.extend(glob.glob('bib/*.bibtex'))

    if not bib_files:
        print("ERROR: No .bib or .bibtex files found in bib/ directory")
        return 1

    # remove duplicates and keep stable order
    bib_files = list(dict.fromkeys(bib_files))
    print(f"Found {len(bib_files)} bib file(s)")

    all_entries = []

    for bibfile in bib_files:
        try:
            with open(bibfile, encoding='utf-8') as bibtex_file:
                bib_database = bibtexparser.load(bibtex_file)
                if not bib_database.entries:
                    print(f"Warning: {bibfile} contains no entries")
                else:
                    print(f"  ✓ {bibfile}: {len(bib_database.entries)} entry(ies)")
                all_entries.extend(bib_database.entries)
        except FileNotFoundError:
            print(f"ERROR: File not found: {bibfile}")
            return 1
        except Exception as e:
            print(f"ERROR: Failed to read {bibfile}: {e}")
            return 1

    if not all_entries:
        print("ERROR: No valid entries found in any bib file")
        return 1

    # Sort newest first (missing years come last)
    def year_key(e):
        try:
            return int(e.get('year', ''))
        except (ValueError, TypeError):
            return -1  # Missing/invalid years sort last (reverse=True)
    all_entries.sort(key=year_key, reverse=True)

    lines = ['## 📚 Publications']
    for entry in all_entries:
        lines.append(format_entry(entry))

    try:
        with open('publications.md', 'w', encoding='utf-8') as out:
            out.write('\n'.join(lines))
        print(f"\n✅ Generated publications.md with {len(all_entries)} publication(s)")
        return 0
    except Exception as e:
        print(f"ERROR: Failed to write publications.md: {e}")
        return 1


if __name__ == '__main__':
    exit(main())
