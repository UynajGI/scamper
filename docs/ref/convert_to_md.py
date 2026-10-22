#!/usr/bin/env python3
"""Convert every PDF in docs/ref/ to Markdown under docs/ref/md/ (incremental).

Uses PyMuPDF (any local pdf->md text extraction is acceptable per project
decision; this script replaces the MinerU-only flow until a token exists).
Page boundaries become `<!-- page N -->` comments so equations and figure
captions stay traceable to their source page.
"""
import sys
from pathlib import Path

import pymupdf

REF = Path(__file__).resolve().parent
OUT = REF / "md"

PAPERS = {
    "NewmanZiff_2000_PRL_percolation": (
        "Efficient Monte Carlo Algorithm and High-Precision Results for Percolation",
        "Phys. Rev. Lett. 85, 4104 (2000)",
        "https://arxiv.org/abs/cond-mat/0005264",
    ),
    "Callaway_2000_PRL_network_robustness": (
        "Network Robustness and Fragility: Percolation on Random Graphs",
        "Phys. Rev. Lett. 85, 5468 (2000)",
        "https://arxiv.org/abs/cond-mat/0007300",
    ),
    "Cohen_2000_PRL_internet_resilience": (
        "Resilience of the Internet to Random Breakdowns",
        "Phys. Rev. Lett. 85, 4626 (2000)",
        "https://arxiv.org/abs/cond-mat/0007048",
    ),
    "Cohen_2001_PRL_intentional_attack": (
        "Breakdown of the Internet under Intentional Attack",
        "Phys. Rev. Lett. 86, 3682 (2001)",
        "https://arxiv.org/abs/cond-mat/0010251",
    ),
}

def main() -> int:
    OUT.mkdir(exist_ok=True)
    for stem, (title, journal, url) in PAPERS.items():
        pdf = REF / f"{stem}.pdf"
        target = OUT / f"{stem}.md"
        if not pdf.exists():
            print(f"skip {stem}: no pdf")
            continue
        if target.exists() and target.stat().st_mtime >= pdf.stat().st_mtime:
            print(f"skip {stem}: up to date")
            continue
        doc = pymupdf.open(pdf)
        parts = [
            f"# {title}",
            "",
            f"> {journal} · [arXiv]({url})",
            f"> Local copy: `docs/ref/{stem}.pdf` ({doc.page_count} pages)",
            "",
        ]
        for index, page in enumerate(doc, start=1):
            text = page.get_text("text").strip()
            parts.append(f"<!-- page {index} -->")
            parts.append("")
            parts.append(text)
            parts.append("")
        target.write_text("\n".join(parts), encoding="utf-8")
        print(f"wrote {target.name} ({doc.page_count} pages)")
    return 0

if __name__ == "__main__":
    sys.exit(main())
