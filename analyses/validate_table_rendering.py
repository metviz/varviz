#!/usr/bin/env python3
"""Verify that every table cell in a markdown source survives into its .docx.

Supplementary Table S3 shipped corrupted for an unknown number of builds. It
used pandoc's fixed-width multiline table format, where the dash rule defines
the column boundaries in characters. Adding reference markers lengthened the
tool names - `InterVar` became `InterVar[8]`, `MARRVEL` became `MARRVEL[10]` -
pushing cell text past those boundaries, so pandoc sliced words at the column
edges and rendered "Open-source ACMG" as "Ope" in one cell and "n-source ACMG"
in the next. The source looked fine to read; only the built document showed it.

Neither the citation validator nor a visual skim of the markdown catches that,
so this checks the two failure modes directly:

  1. fixed-width tables whose rows put text in the gaps between columns, which
     is what silently splits a cell;
  2. pipe tables with ragged rows, where a wrong delimiter count merges or
     drops cells;

and then confirms every source cell appears in the rendered document. The last
check is the one that matters, because it holds whatever the table format is.

Comparison is on normalised text: emphasis markers, superscript carets and
escaped backslashes removed, XML entities decoded, whitespace collapsed.

Usage:
    python3 analyses/validate_table_rendering.py \
        analyses/humu/supp_humu.md docs/Metpally_VarViz_HumMutat_Supplementary.docx
"""
import argparse
import html
import re
import sys
import zipfile

PIPE_ROW = re.compile(r"\s*\|.*\|\s*$")


def is_rule(line):
    """A fixed-width table's dash rule, which defines its column boundaries."""
    return bool(line and line.strip()) and set(line.strip()) <= set("- ")


def normalise(cell):
    cell = cell.replace("**", "").replace("*", "").replace("^", "")
    cell = cell.replace("\\\\", "\\")
    return re.sub(r"\s+", " ", html.unescape(cell)).strip()


def docx_cells(path):
    xml = zipfile.ZipFile(path).read("word/document.xml").decode("utf8", "ignore")
    xml = re.sub(r"</w:p>", "\n", xml)
    xml = re.sub(r"<[^>]+>", "", xml)
    return {normalise(s) for s in html.unescape(xml).split("\n") if s.strip()}


def fixed_width_violations(lines):
    """Rows of a fixed-width table with text in the gaps between columns."""
    out = []
    for r, line in enumerate(lines):
        if not is_rule(line):
            continue
        spans = [(m.start(), m.end()) for m in re.finditer(r"-+", line)]
        if len(spans) < 2:
            continue
        gaps = [(spans[i][1], spans[i + 1][0]) for i in range(len(spans) - 1)]
        for k in range(r + 1, len(lines)):
            if is_rule(lines[k]):
                break
            row = lines[k]
            if not row.strip():
                continue
            for a, b in gaps:
                if row[a:b].strip():
                    out.append((k + 1, a, b, row[a:b], row.strip()[:70]))
    return out


def pipe_tables(lines):
    rows = [k for k, l in enumerate(lines) if PIPE_ROW.match(l)]
    groups = []
    for k in rows:
        if groups and k == groups[-1][-1] + 1:
            groups[-1].append(k)
        else:
            groups.append([k])
    return groups


def cells_of(line):
    s = line.strip().strip("|")
    return [c for c in s.split("|")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("markdown")
    ap.add_argument("docx")
    args = ap.parse_args()

    lines = open(args.markdown).read().splitlines()
    rendered = docx_cells(args.docx)
    problems = 0

    viol = fixed_width_violations(lines)
    print(f"fixed-width tables: {len(viol)} boundary violation(s)")
    for ln, a, b, text, row in viol:
        problems += 1
        print(f"  line {ln}: gap {a}-{b} holds {text!r}  in  {row!r}")

    groups = pipe_tables(lines)
    ragged = 0
    for g in groups:
        widths = {len(cells_of(lines[k])) for k in g
                  if not set(lines[k].strip()) <= set("-| :")}
        if len(widths) > 1:
            ragged += 1
            problems += 1
            print(f"  RAGGED pipe table at line {g[0] + 1}: cell counts {sorted(widths)}")
    print(f"pipe tables: {len(groups)} found, {ragged} ragged")

    # The check that holds regardless of table format.
    total = found = 0
    for g in groups:
        for k in g:
            if set(lines[k].strip()) <= set("-| :"):
                continue
            for c in cells_of(lines[k]):
                c = normalise(c)
                if not c:
                    continue
                total += 1
                if c in rendered:
                    found += 1
                else:
                    problems += 1
                    print(f"  line {k + 1}: cell {c!r} does not appear in the .docx")
    print(f"pipe-table cells present in the rendered document: {found}/{total}")

    if problems:
        print(f"\nFAIL - {problems} problem(s)")
        sys.exit(1)
    print("\nPASS - every table cell survives into the rendered document")


if __name__ == "__main__":
    main()
