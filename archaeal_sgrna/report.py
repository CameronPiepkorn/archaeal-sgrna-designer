
from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

from .designer import GuideRNA


TSV_HEADER = "\t".join([
    "seq_id", "strand", "start", "spacer", "pam",
    "gc_content", "score", "warnings",
])


def to_tsv(guides: Iterable[GuideRNA], path: str | Path | None = None) -> str:
    rows = [TSV_HEADER]
    for g in guides:
        rows.append("\t".join([
            g.chromosome,
            g.strand,
            str(g.start),
            g.sequence,
            g.pam,
            f"{g.gc_content:.3f}",
            f"{g.score:.2f}",
            "; ".join(g.warnings) if g.warnings else "",
        ]))
    content = "\n".join(rows) + "\n"

    if path:
        Path(path).write_text(content)
        return str(path)
    return content


def to_json(guides: Iterable[GuideRNA], path: str | Path | None = None) -> str:
    records = [
        {
            "seq_id": g.chromosome,
            "strand": g.strand,
            "start": g.start,
            "spacer": g.sequence,
            "pam": g.pam,
            "gc_content": round(g.gc_content, 4),
            "score": round(g.score, 2),
            "warnings": g.warnings,
        }
        for g in guides
    ]
    content = json.dumps(records, indent=2)

    if path:
        Path(path).write_text(content)
        return str(path)
    return content


def summary(guides: list[GuideRNA], top_n: int = 10) -> str:
    lines = [
        "=" * 72,
        f"  Archaeal sgRNA Designer — Top {min(top_n, len(guides))} Guides",
        "=" * 72,
        f"{'#':<4} {'Spacer (5→3)':<22} {'PAM':<6} "
        f"{'GC%':<7} {'Score':<7} Warnings",
        "-" * 72,
    ]
    for i, g in enumerate(guides[:top_n], 1):
        warn = ", ".join(g.warnings) if g.warnings else "—"
        lines.append(
            f"{i:<4} {g.sequence:<22} {g.pam:<6} "
            f"{g.gc_content:>5.0%}   {g.score:>6.1f}   {warn}"
        )
    lines.append("=" * 72)
    return "\n".join(lines)
