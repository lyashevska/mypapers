from pathlib import Path


def main():
    p = Path('publications.md')
    if not p.exists():
        print('ERROR: publications.md not found')
        return 2

    text = p.read_text(encoding='utf-8')
    if '## 📚 Publications' not in text:
        print('ERROR: title line "## 📚 Publications" not found')
        return 3

    # very small heuristic: at least one bullet line starting with '- '
    if '\n- ' not in text and not text.strip().startswith('- '):
        print('ERROR: no entries found in publications.md')
        return 4

    print('OK: publications.md looks valid')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
