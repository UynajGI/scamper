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
    "Hinrichsen_2000_AdvPhys_absorbing_states": (
        "Nonequilibrium Critical Phenomena and Phase Transitions into Absorbing States",
        "Adv. Phys. 49, 815 (2000)",
        "https://arxiv.org/abs/cond-mat/0001070",
    ),
    "Cohen_2001_PRL_intentional_attack": (
        "Breakdown of the Internet under Intentional Attack",
        "Phys. Rev. Lett. 86, 3682 (2001)",
        "https://arxiv.org/abs/cond-mat/0010251",
    ),
    "Holroyd_2003_PTRF_bootstrap_sharp_threshold": (
        "Sharp Metastability Threshold for Two-Dimensional Bootstrap Percolation",
        "Probab. Theory Relat. Fields 125, 195-224 (2003)",
        "https://arxiv.org/abs/math/0206132",
    ),
    "Janson_2007_kcore_simple_solution": (
        "A Simple Solution to the k-Core Problem",
        "Random Structures Algorithms 30, 50-62 (2007)",
        "https://arxiv.org/abs/math/0508453",
    ),
    "Riordan_2012_Science_achlioptas_continuous": (
        "Achlioptas Process Phase Transitions Are Continuous",
        "Science 335, 1189-1191 (2012)",
        "https://arxiv.org/abs/1102.5306",
    ),
    "Mertens_2012_PRE_continuum_percolation_thresholds_2d": (
        "Continuum Percolation Thresholds in Two Dimensions",
        "Phys. Rev. E 86, 061109 (2012)",
        "https://arxiv.org/abs/1209.4936",
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
