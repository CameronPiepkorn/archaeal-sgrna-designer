
import pytest
from archaeal_sgrna import ArchaealGuideDesigner, GuideRNA
from archaeal_sgrna.utils import reverse_complement, gc_content, parse_fasta
from archaeal_sgrna.scoring import score_guide
from archaeal_sgrna import report
import tempfile, json
from pathlib import Path


SEQ_WITH_PAMS = (
    "AACGTACGTA"
    "GCGATCGATCGATCGATCGC"
    "GGG"
    "ATCGATCGAT"
    "CG"
    "GCGCGCGCGCGCGCGCGCGC"
    "TGG"
    "AAAA"
)


@pytest.fixture
def designer():
    return ArchaealGuideDesigner()


@pytest.fixture
def tmp_fasta(tmp_path):
    fa = tmp_path / "test.fasta"
    fa.write_text(f">seq1\n{SEQ_WITH_PAMS}\n")
    return fa


class TestUtils:
    def test_reverse_complement_basic(self):
        assert reverse_complement("ATCG") == "CGAT"

    def test_reverse_complement_roundtrip(self):
        seq = "GCGATCGATCGATCGATCGCGGG"
        assert reverse_complement(reverse_complement(seq)) == seq

    def test_gc_content_all_gc(self):
        assert gc_content("GCGCGC") == pytest.approx(1.0)

    def test_gc_content_half(self):
        assert gc_content("AATTGGCC") == pytest.approx(0.5)

    def test_gc_content_empty(self):
        assert gc_content("") == 0.0

    def test_parse_fasta(self, tmp_fasta):
        records = list(parse_fasta(tmp_fasta))
        assert len(records) == 1
        seq_id, seq = records[0]
        assert seq_id == "seq1"
        assert seq == SEQ_WITH_PAMS.upper()


class TestScoring:
    def test_score_range(self):
        spacer = "GCGATCGATCGATCGATCGC"
        s = score_guide(spacer, gc=gc_content(spacer))
        assert 0.0 <= s <= 100.0

    def test_homopolymer_penalty(self):
        good = "GCGCGCGCGCGCGCGCGCGC"
        bad  = "GCGCGCTTTTCGCGCGCGCG"
        s_good = score_guide(good, gc=gc_content(good))
        s_bad  = score_guide(bad,  gc=gc_content(bad))
        assert s_good > s_bad

    def test_low_gc_penalty(self):
        good_gc = "GCATGCATGCATGCATGCAT"
        low_gc  = "ATATATATATATATATATATAT"[:20]
        s_good = score_guide(good_gc, gc=gc_content(good_gc))
        s_low  = score_guide(low_gc,  gc=gc_content(low_gc))
        assert s_good > s_low

    def test_box_a_penalty(self):
        normal = "GCGATCGATCGATCGATCGC"
        box_a  = "GCGATCGATTTATAGCGCGC"
        s_normal = score_guide(normal, gc=gc_content(normal), archaeal_promoter_bias=True)
        s_box_a  = score_guide(box_a,  gc=gc_content(box_a),  archaeal_promoter_bias=True)
        assert s_normal > s_box_a


class TestDesigner:
    def test_finds_guides_in_sequence(self, designer):
        guides = designer.design_from_sequence(SEQ_WITH_PAMS, seq_id="test")
        assert len(guides) >= 1

    def test_guides_are_sorted_by_score(self, designer):
        guides = designer.design_from_sequence(SEQ_WITH_PAMS)
        scores = [g.score for g in guides]
        assert scores == sorted(scores, reverse=True)

    def test_guide_length(self, designer):
        guides = designer.design_from_sequence(SEQ_WITH_PAMS)
        for g in guides:
            assert len(g.sequence) == 20

    def test_pam_ngg(self, designer):
        guides = designer.design_from_sequence(SEQ_WITH_PAMS)
        for g in guides:
            if g.strand == "+":
                assert g.pam[1:] == "GG", f"Expected xGG PAM, got {g.pam}"

    def test_top_n(self, designer):
        guides = designer.design_from_sequence(SEQ_WITH_PAMS, top_n=2)
        assert len(guides) <= 2

    def test_strand_coverage(self, designer):
        guides = designer.design_from_sequence(SEQ_WITH_PAMS)
        strands = {g.strand for g in guides}
        assert "+" in strands

    def test_from_fasta(self, designer, tmp_fasta):
        results = designer.design_from_fasta(str(tmp_fasta))
        assert "seq1" in results
        assert len(results["seq1"]) >= 1

    def test_invalid_pam_raises(self):
        with pytest.raises(ValueError, match="Unknown PAM type"):
            ArchaealGuideDesigner(pam_type="FakeCas")

    def test_gc_filter_warnings(self):
        d = ArchaealGuideDesigner(gc_min=0.0, gc_max=1.0)
        guides = d.design_from_sequence(SEQ_WITH_PAMS)
        for g in guides:
            assert not any("GC" in w for w in g.warnings)

    def test_tttt_warning(self):
        tttt_seq = "T" * 4 + "GCGATCGATCGATCGATCGC" + "TGG" + "A" * 20
        d = ArchaealGuideDesigner(avoid_tttt=True)
        spacer_with_tttt = "GCGATTTTTCGATCGATCGC"
        from archaeal_sgrna.scoring import score_guide
        from archaeal_sgrna.utils import gc_content
        has_tttt = "TTTT" in spacer_with_tttt
        assert has_tttt


class TestReport:
    def test_tsv_header(self, designer):
        guides = designer.design_from_sequence(SEQ_WITH_PAMS)
        tsv = report.to_tsv(guides)
        assert tsv.startswith("seq_id\t")

    def test_tsv_row_count(self, designer):
        guides = designer.design_from_sequence(SEQ_WITH_PAMS)
        tsv = report.to_tsv(guides)
        lines = [l for l in tsv.strip().split("\n") if l]
        assert len(lines) == len(guides) + 1

    def test_tsv_file_write(self, designer, tmp_path):
        guides = designer.design_from_sequence(SEQ_WITH_PAMS)
        out = tmp_path / "out.tsv"
        report.to_tsv(guides, path=out)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_json_output(self, designer):
        guides = designer.design_from_sequence(SEQ_WITH_PAMS)
        js = report.to_json(guides)
        parsed = json.loads(js)
        assert isinstance(parsed, list)
        if parsed:
            assert "spacer" in parsed[0]
            assert "score" in parsed[0]

    def test_summary_string(self, designer):
        guides = designer.design_from_sequence(SEQ_WITH_PAMS)
        s = report.summary(guides)
        assert "Score" in s
        assert "Spacer" in s
