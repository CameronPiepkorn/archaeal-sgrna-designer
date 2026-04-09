
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Iterator, Optional

from .scoring import score_guide
from .utils import reverse_complement, parse_fasta


PAM_PATTERNS: dict[str, str] = {
    "SpCas9":  r"(?=[ACGT]GG)",
    "SaCas9":  r"(?=[ACGT][ACGT]G[AG][AG]T)",
    "AsCas12a": r"(?=TTT[ACG])",
}

DEFAULT_GUIDE_LENGTH = 20


@dataclass
class GuideRNA:

    sequence: str
    pam: str
    chromosome: str
    start: int
    strand: str
    gc_content: float
    score: float
    warnings: list[str] = field(default_factory=list)

    @property
    def full_spacer_with_pam(self) -> str:
        # For 5' PAM nucleases (Cas12a), PAM precedes the spacer
        return self.pam + self.sequence if len(self.pam) == 4 else self.sequence + self.pam

    def __repr__(self) -> str:
        warn = f"  ⚠ {'; '.join(self.warnings)}" if self.warnings else ""
        return (
            f"<GuideRNA {self.sequence}|{self.pam}  "
            f"chr={self.chromosome} pos={self.start}{self.strand} "
            f"GC={self.gc_content:.0%} score={self.score:.1f}{warn}>"
        )


class ArchaealGuideDesigner:

    def __init__(
        self,
        pam_type: str = "SpCas9",
        guide_length: int = DEFAULT_GUIDE_LENGTH,
        gc_min: float = 0.30,
        gc_max: float = 0.80,
        avoid_tttt: bool = True,
        archaeal_promoter_bias: bool = True,
    ) -> None:
        if pam_type not in PAM_PATTERNS:
            raise ValueError(
                f"Unknown PAM type '{pam_type}'. "
                f"Choose from: {list(PAM_PATTERNS)}"
            )
        self.pam_type = pam_type
        self.guide_length = guide_length
        self.gc_min = gc_min
        self.gc_max = gc_max
        self.avoid_tttt = avoid_tttt
        self.archaeal_promoter_bias = archaeal_promoter_bias
        self._pam_re = re.compile(PAM_PATTERNS[pam_type], re.IGNORECASE)


    def design_from_sequence(
        self,
        sequence: str,
        seq_id: str = "input",
        top_n: Optional[int] = None,
    ) -> list[GuideRNA]:
        sequence = sequence.upper().replace(" ", "").replace("\n", "")
        guides: list[GuideRNA] = []

        for guide in self._extract_guides(sequence, seq_id, strand="+"):
            guides.append(guide)

        rc_seq = reverse_complement(sequence)
        for guide in self._extract_guides(rc_seq, seq_id, strand="-"):
            guides.append(guide)

        guides.sort(key=lambda g: g.score, reverse=True)
        return guides[:top_n] if top_n else guides

    def design_from_fasta(
        self,
        fasta_path: str,
        top_n: Optional[int] = None,
    ) -> dict[str, list[GuideRNA]]:
        results: dict[str, list[GuideRNA]] = {}
        for seq_id, sequence in parse_fasta(fasta_path):
            results[seq_id] = self.design_from_sequence(
                sequence, seq_id=seq_id, top_n=top_n
            )
        return results

    def design_from_gene_region(
        self,
        sequence: str,
        target_start: int,
        target_end: int,
        seq_id: str = "input",
        flank: int = 200,
        top_n: Optional[int] = None,
    ) -> list[GuideRNA]:
        start = max(0, target_start - flank)
        end = min(len(sequence), target_end + flank)
        region = sequence[start:end]
        guides = self.design_from_sequence(region, seq_id=seq_id, top_n=top_n)
        for g in guides:
            g.start += start
        return guides


    # PAMs that come AFTER the spacer (3' PAM)
    _PAM_IS_3PRIME = {"SpCas9", "SaCas9"}
    # PAMs that come BEFORE the spacer (5' PAM)
    _PAM_IS_5PRIME = {"AsCas12a"}

    def _extract_guides(
        self, sequence: str, seq_id: str, strand: str
    ) -> Iterator[GuideRNA]:
        pam_len = self._pam_length()
        gl = self.guide_length
        five_prime_pam = self.pam_type in self._PAM_IS_5PRIME

        for match in self._pam_re.finditer(sequence):
            pam_start = match.start()

            if five_prime_pam:
                # PAM is 5' of spacer: [PAM][spacer]
                spacer_start = pam_start + pam_len
                spacer_end = spacer_start + gl
                if spacer_end > len(sequence):
                    continue
                spacer = sequence[spacer_start:spacer_end]
                pam_seq = sequence[pam_start:pam_start + pam_len]
                guide_start = pam_start
            else:
                # PAM is 3' of spacer: [spacer][PAM]
                guide_start = pam_start - gl
                if guide_start < 0:
                    continue
                if pam_start + pam_len > len(sequence):
                    continue
                spacer = sequence[guide_start:pam_start]
                pam_seq = sequence[pam_start:pam_start + pam_len]

            if len(spacer) != gl:
                continue
            if not re.fullmatch(r"[ACGT]+", spacer):
                continue

            gc = (spacer.count("G") + spacer.count("C")) / gl

            warnings: list[str] = []
            if gc < self.gc_min:
                warnings.append(f"Low GC ({gc:.0%})")
            if gc > self.gc_max:
                warnings.append(f"High GC ({gc:.0%})")
            if self.avoid_tttt and "TTTT" in spacer:
                warnings.append("TTTT run (Pol-III terminator)")

            composite = score_guide(
                spacer,
                gc=gc,
                archaeal_promoter_bias=self.archaeal_promoter_bias,
            )

            yield GuideRNA(
                sequence=spacer,
                pam=pam_seq,
                chromosome=seq_id,
                start=guide_start,
                strand=strand,
                gc_content=gc,
                score=composite,
                warnings=warnings,
            )

    def _pam_length(self) -> int:
        lengths = {"SpCas9": 3, "SaCas9": 6, "AsCas12a": 4}
        return lengths[self.pam_type]
