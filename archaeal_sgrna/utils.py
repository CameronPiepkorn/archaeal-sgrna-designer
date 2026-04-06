
from __future__ import annotations

from pathlib import Path
from typing import Iterator


_COMPLEMENT: dict[str, str] = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def gc_content(seq: str) -> float:
    seq = seq.upper()
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def parse_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    path = Path(path)
    current_id: str | None = None
    buffer: list[str] = []

    with path.open() as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line or line.startswith(";"):
                continue
            if line.startswith(">"):
                if current_id is not None:
                    yield current_id, "".join(buffer)
                current_id = line[1:].split()[0]
                buffer = []
            else:
                buffer.append(line.upper())

    if current_id is not None:
        yield current_id, "".join(buffer)


def validate_dna(seq: str) -> bool:
    import re
    return bool(re.fullmatch(r"[ACGTacgtNnRrYySwWkKmMbBdDhHvV]+", seq))
