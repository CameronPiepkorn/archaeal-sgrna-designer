[![CI](https://github.com/...)]
[![Python 3.9+](...)]
[![License: MIT](...)]
[![PyPI](https://img.shields.io/pypi/v/archaeal-sgrna-designer.svg)](https://pypi.org/project/archaeal-sgrna-designer/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19506160.svg)](https://doi.org/10.5281/zenodo.19506160)
# archaeal-sgrna-designer

[![CI](https://github.com/CameronPiepkorn/archaeal-sgrna-designer/actions/workflows/ci.yml/badge.svg)](https://github.com/CameronPiepkorn/archaeal-sgrna-designer/actions)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A Python library for designing and ranking sgRNA candidates for archaeal
CRISPR-Cas9 genome editing, with scoring heuristics tuned for *Methanosarcina acetivorans*
and related methanogens.

Built to support research following:

> Nayak & Metcalf (2017). *Cas9-mediated genome editing in the methanogenic
> archaeon Methanosarcina acetivorans.* PNAS 114(11):2976–2981.

> Chadwick, Joiner, Ramesh, Mitchell & Nayak (2023). *McrD binds asymmetrically
> to methyl-coenzyme M reductase improving active-site accessibility during
> assembly.* PNAS 120(25):e2302815120.

---

## Features

- **PAM-aware scanning** — supports SpCas9 (NGG), SaCas9 (NNGRRT), and AsCas12a (TTTV)
- **Archaeal-specific scoring** — penalises internal Box A (TATA-box) promoter motifs and
  TTTT runs that terminate archaeal transcription
- **Strand-aware** — finds guides on both strands and returns genomic coordinates
- **Zero mandatory dependencies** — pure Python stdlib; no BioPython required
- **Flexible input** — raw sequence string, FASTA file, or windowed gene region
- **Export** — TSV, JSON, or formatted text summary

---

## Installation

```bash
# From PyPI (once published)
pip install archaeal-sgrna-designer

# From source
git clone https://github.com/CameronPiepkorn/archaeal-sgrna-designer
cd archaeal-sgrna-designer
pip install -e ".[dev]"
```

---

## Quick Start

```python
from archaeal_sgrna import ArchaealGuideDesigner, report

# Initialise with SpCas9 (validated in M. acetivorans, Nayak & Metcalf 2017)
designer = ArchaealGuideDesigner(
    pam_type="SpCas9",
    guide_length=20,
    gc_min=0.35,
    gc_max=0.75,
)

# Design from a raw sequence
guides = designer.design_from_sequence(my_dna_string, seq_id="mcrA", top_n=10)

# Pretty-print top results
print(report.summary(guides))

# Export
report.to_tsv(guides, path="mcrA_guides.tsv")
report.to_json(guides, path="mcrA_guides.json")
```

### Design from a FASTA file

```python
results = designer.design_from_fasta("my_genes.fasta", top_n=5)
for seq_id, guides in results.items():
    print(f"\n=== {seq_id} ===")
    print(report.summary(guides))
```

### Target a specific gene region

```python
guides = designer.design_from_gene_region(
    full_genome_seq,
    target_start=125_000,   # 0-based CDS start
    target_end=126_500,     # 0-based CDS end
    seq_id="MA0528_mcrA",
    flank=200,              # also search ±200 bp of coding region
    top_n=10,
)
```

---

## Scoring

Each candidate is assigned a heuristic filter score in **[0, 100]** (higher = fewer
red-flag properties). This is **not** a validated on-target efficiency predictor —
no quantitative model has been published for Cas9 activity in *Methanosarcina* or
any other archaeon. Use the score to deprioritise guides with known problems;
experimental validation is always required.

| Rule | Biological basis |
|---|---|
| GC 30–70% window | Guides outside this range show reduced stability or secondary structure issues (general CRISPR literature; no Methanosarcina-specific data) |
| Homopolymer penalty (≥4 identical bases) | TTTT terminates archaeal Pol-III transcription (Santangelo et al. 2009, J Bacteriol 191:7102); other homopolymers risk secondary structure |
| Archaeal Box A (TTTATA) penalty | Exact consensus TTTATA is the archaeal TATA-box recognised by TBP; a spacer containing it may produce competing transcripts (Qureshi & Jackson 1998, Mol Cell 1:389) |

> **What was removed:** An earlier version of this tool included a position-weight
> matrix falsely attributed to Doench et al. (2016), a fabricated mismatch-tolerance
> claim about *Methanosarcina* Cas9, and a "5' G bonus" incorrectly referencing the
> eukaryotic U6 promoter. Those have been removed. If archaeal-specific training data
> become available, a proper model can replace these heuristics.

---

## Archaeal-Specific Considerations

### Box A (TATA-box) promoter elements
The archaeal TATA-box has the consensus sequence **TTTATA** and is bound by TBP
to initiate RNA Pol transcription (Qureshi & Jackson 1998, Mol Cell 1:389).
When the sgRNA spacer contains this hexamer, it may drive aberrant transcription
that competes with the intended guide RNA. This tool flags and penalises spacers
containing the exact TTTATA motif. Note: Box B is a separate element found in
tRNA genes, not relevant to most coding-gene targets.

### TTTT runs
Runs of four or more T residues act as Pol-III transcription terminators and
reduce sgRNA expression (Santangelo et al. 2009, J Bacteriol 191:7102). Guides
containing TTTT are flagged and penalised.

### PAM validation
SpCas9 (NGG PAM) is the only system experimentally validated for editing in
*Methanosarcina* (Nayak & Metcalf 2017). SaCas9 and AsCas12a PAM options are
included as convenience features but have not been tested in any archaeon.

---

## Example Files

| File | Description |
|---|---|
| `examples/mcrA_mcrB_example.fasta` | Partial *mcrA* and *mcrB* sequences from *M. acetivorans* C2A |
| `examples/basic_usage.py` | End-to-end usage walkthrough |

Run the example:
```bash
python examples/basic_usage.py
```

---

## Running Tests

```bash
pytest                        # all tests
pytest --cov=archaeal_sgrna  # with coverage
```

---

## Repository Structure

```
archaeal-sgrna-designer/
├── archaeal_sgrna/
│   ├── __init__.py       # public API
│   ├── designer.py       # ArchaealGuideDesigner, GuideRNA
│   ├── scoring.py        # composite on-target scoring
│   ├── report.py         # TSV / JSON / text export
│   └── utils.py          # sequence utilities, FASTA parser
├── tests/
│   └── test_designer.py
├── examples/
│   ├── mcrA_mcrB_example.fasta
│   └── basic_usage.py
├── .github/workflows/ci.yml
├── pyproject.toml
└── README.md
```

---

## Contributing

Pull requests welcome. Please:
1. Add tests for new features
2. Keep external dependencies optional (stdlib-first)
3. Cite any new scoring heuristics you introduce

---

## Citation

If you use this tool in published research, please cite:

```bibtex
@article{nayak2017cas9,
  title   = {Cas9-mediated genome editing in the methanogenic archaeon
             Methanosarcina acetivorans},
  author  = {Nayak, Dipti D and Metcalf, William W},
  journal = {Proceedings of the National Academy of Sciences},
  volume  = {114},
  number  = {11},
  pages   = {2976--2981},
  year    = {2017}
}
```

---

## License

MIT © 2024. See [LICENSE](LICENSE).
