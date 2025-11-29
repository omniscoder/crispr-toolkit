from crispr_toolkit.grna import GuideCheckResult, validate_guide


def test_good_guide_passes():
    seq = "ACGTACGTACGTACGTACGT"  # 20 nt, 50% GC, no long homopolymer
    result = validate_guide(seq)
    assert isinstance(result, GuideCheckResult)
    assert result.is_valid
    assert 0.3 <= result.gc_content <= 0.8
    assert result.length == len(seq)


def test_too_short_and_too_long():
    short = validate_guide("ACGT")
    assert not short.is_valid
    assert any("length" in r for r in short.reasons)

    long = validate_guide("ACGT" * 10)
    assert not long.is_valid
    assert any("length" in r for r in long.reasons)


def test_gc_out_of_bounds():
    low_gc = validate_guide("ATATATATATATATATATAT")
    assert not low_gc.is_valid
    assert any("GC content" in r for r in low_gc.reasons)

    high_gc = validate_guide("GCGCGCGCGCGCGCGCGCGC")
    assert not high_gc.is_valid
    assert any("GC content" in r for r in high_gc.reasons)


def test_homopolymer_fail():
    seq = "AAAACGTACGTACGTACGT"  # leading A4 homopolymer
    result = validate_guide(seq)
    assert not result.is_valid
    assert any("homopolymer" in r for r in result.reasons)


def test_invalid_bases():
    result = validate_guide("ACGTNNNNACGT")
    assert not result.is_valid
    assert any("invalid bases" in r for r in result.reasons)

