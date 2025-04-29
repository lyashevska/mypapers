import bibtexparser
import os

def bibtex_to_markdown(entry):
    authors = entry.get("author", "Unknown").replace("\n", " ")
    title = entry.get("title", "No title").strip("{}")
    year = entry.get("year", "n.d.")
    journal = entry.get("journal") or entry.get("booktitle", "")
    doi = entry.get("doi", "")
    url = entry.get("url", "")
    abstract = entry.get("abstract", "").strip()

    citation = f"- **{title}**  \n  {authors} ({year})."
    if journal:
        citation += f" _{journal}_."

    if doi:
        citation += f" [DOI](https://doi.org/{doi})"
    elif url:
        citation += f" [Link]({url})"

    if abstract:
        citation += f"\n\n  > {abstract}"

    return citation

def convert_bib_files_to_md():
    # Get all .bib files in the current directory
    bib_files = [f for f in os.listdir('.') if f.endswith('.bib')]
    all_entries = []

    for bib_file in bib_files:
        with open(bib_file, encoding="utf-8") as bibtex_file:
            bib_database = bibtexparser.load(bibtex_file)
            all_entries.extend(bib_database.entries)

    # Sort by year (descending), fallback to 0 if year missing or invalid
    def safe_year(entry):
        try:
            return int(entry.get("year", 0))
        except ValueError:
            return 0

    all_entries.sort(key=safe_year, reverse=True)

    # Write to Markdown file
    with open("publications.md", "w", encoding="utf-8") as out_md:
        out_md.write("# Publications\n\n")
        for entry in all_entries:
            out_md.write(bibtex_to_markdown(entry) + "\n\n")

if __name__ == "__main__":
    convert_bib_files_to_md()

