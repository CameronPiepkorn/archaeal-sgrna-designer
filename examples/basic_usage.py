
from archaeal_sgrna import ArchaealGuideDesigner, report

print("=" * 60)
print("Example 1: Design from a raw DNA sequence")
print("=" * 60)

mcrA_fragment = (
    "ATGAGCAACAAACCGGAAACCGAGGGCAAAGACAAACTCCCGGAGCAGTT"
    "CGAAATCGCGACCGATTTCGCCAAGGCGCTCAAGGTCGACGGCGCCGAG"
    "GCGGTCATCGTCAACGGCACCTCCGACGCCGGCATCGTCAAGGCCGTGG"
    "TGGAAGGCGAGGTCATGAACCTCGACACCATCCGGCGCATGGACGCCGT"
)

designer = ArchaealGuideDesigner(
    pam_type="SpCas9",
    guide_length=20,
    gc_min=0.35,
    gc_max=0.75,
)

guides = designer.design_from_sequence(mcrA_fragment, seq_id="mcrA_frag", top_n=5)
print(report.summary(guides))

print("\n" + "=" * 60)
print("Example 2: Design from a FASTA file")
print("=" * 60)

fasta_path = "examples/mcrA_mcrB_example.fasta"
results = designer.design_from_fasta(fasta_path, top_n=3)

for seq_id, guides in results.items():
    print(f"\n--- {seq_id} ---")
    print(report.summary(guides, top_n=3))

print("\n" + "=" * 60)
print("Example 3: Export results")
print("=" * 60)

all_guides = []
for guides in results.values():
    all_guides.extend(guides)

tsv_str = report.to_tsv(all_guides)
print("TSV output (first 3 lines):")
for line in tsv_str.split("\n")[:3]:
    print(" ", line)

import json
json_str = report.to_json(all_guides)
records = json.loads(json_str)
print(f"\nJSON: {len(records)} guides serialised")
if records:
    print(f"  First guide: {records[0]}")

print("\n" + "=" * 60)
print("Example 4: Target a specific gene region")
print("=" * 60)

full_seq = (
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    "ATGAGCAACAAACCGGAAACCGAGGGCAAAGACAAACTCCCGGAGCAGTT"
    "CGAAATCGCGACCGATTTCGCCAAGGCGCTCAAGGTCGACGGCGCCGAG"
    "GCGGTCATCGTCAACGGCACCTCCGACGCCGGCATCGTCAAGGCCGTGG"
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
)

guides = designer.design_from_gene_region(
    full_seq,
    target_start=57,
    target_end=57 + 147,
    seq_id="mcrA_full",
    flank=50,
    top_n=3,
)
print(report.summary(guides, top_n=3))
print("\nDone! See report.to_tsv() / report.to_json() to save results.")
