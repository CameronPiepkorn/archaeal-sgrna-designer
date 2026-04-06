
from __future__ import annotations


_BOX_A_CONSENSUS = "TTTATA"

_BOX_A_PENALTY = 15.0
_HOMOPOLYMER_PENALTY = 10.0
_GC_PENALTY_PER_UNIT = 25.0

_GC_OPTIMAL_LOW = 0.30
_GC_OPTIMAL_HIGH = 0.70


def score_guide(
    spacer: str,
    gc: float,
    archaeal_promoter_bias: bool = True,
) -> float:
    spacer = spacer.upper()
    score = 60.0

    if gc < _GC_OPTIMAL_LOW:
        score -= (_GC_OPTIMAL_LOW - gc) * _GC_PENALTY_PER_UNIT
    elif gc > _GC_OPTIMAL_HIGH:
        score -= (gc - _GC_OPTIMAL_HIGH) * _GC_PENALTY_PER_UNIT

    for base in ("AAAA", "CCCC", "GGGG", "TTTT"):
        if base in spacer:
            score -= _HOMOPOLYMER_PENALTY

    if archaeal_promoter_bias and _BOX_A_CONSENSUS in spacer:
        score -= _BOX_A_PENALTY

    return float(max(0.0, min(100.0, score)))
