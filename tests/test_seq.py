import pytest

from crispr_toolkit.seq import gc_content, has_homopolymer, normalize_seq, reverse_complement


def test_normalize_seq_basic():
    assert normalize_seq("acgt") == "ACGT"
    assert normalize_seq(" ac gt \n") == "ACGT"


def test_normalize_seq_rna_to_dna():
    assert normalize_seq("augu") == "ATGT"


def test_reverse_complement_simple():
    assert reverse_complement("ACGTN") == "NACGT"


def test_gc_content_edges():
    assert gc_content("") == 0.0
    assert gc_content("AAAA") == 0.0
    assert gc_content("GGGG") == 1.0
    assert pytest.approx(gc_content("ACGT"), 0.001) == 0.5


def test_gc_content_with_odd_chars():
    # Non-ACGT bases reduce the apparent GC fraction.
    val = gc_content("ACGTNNNN")
    assert 0.1 < val < 0.5


def test_has_homopolymer_default_length():
    assert has_homopolymer("AAAAC") is True
    assert has_homopolymer("CAAAG") is False
    assert has_homopolymer("TTTT") is True


def test_has_homopolymer_custom_length():
    assert has_homopolymer("AAAG", length=3) is True
    assert has_homopolymer("AAGA", length=3) is False


def test_has_homopolymer_invalid_length():
    with pytest.raises(ValueError):
        has_homopolymer("ACGT", length=1)

