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
    "Auffinger_2017_50yrs_first_passage_percolation": (
        "50 years of first passage percolation",
        "AMS Univ. Lecture Series 68 (2017)",
        "https://arxiv.org/abs/1511.03262",
    ),
    "Biskup_2003_long_range_chemical_distance": (
        "On the scaling of the chemical distance in long-range percolation",
        "Ann. Probab. 32, 2938 (2004)",
        "https://arxiv.org/abs/math/0304418",
    ),
    "Geyer_2015_sweeny_random_cluster": (
        "Efficient simulation of the random-cluster model",
        "J. Stat. Mech. (2015) P08010",
        "https://arxiv.org/abs/1307.6647",
    ),
    "Schrenk_2013_correlated_disorder_percolation": (
        "Percolation with long-range correlated disorder",
        "Phys. Rev. E 88, 052102 (2013)",
        "https://arxiv.org/abs/1309.2994",
    ),
    "Pruessner_2018_correlated_thresholds": (
        "Percolation thresholds and fractal dimensions for square/cubic lattices with correlated disorder",
        "Phys. Rev. E 97, 062141 (2018)",
        "https://arxiv.org/abs/1708.02296",
    ),
    "Saberi_2015_percolation_review": (
        "Recent advances and open challenges in percolation",
        "Eur. Phys. J. B 89, 253 (2016)",
        "https://arxiv.org/abs/1404.5325",
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
