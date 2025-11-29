from crispr_toolkit.scoring_basic import basic_off_target_penalty, basic_on_target_score


def test_basic_on_target_score_prefers_reasonable_guides():
    ideal = "ACGTACGTACGTACGTACGT"  # 20 nt, ~50% GC, no homopolymers
    messy = "GGGGGGGGGGGGGGGGGGGG"  # extreme GC and homopolymer

    ideal_score = basic_on_target_score(ideal)
    messy_score = basic_on_target_score(messy)

    assert 0.0 <= ideal_score <= 1.0
    assert 0.0 <= messy_score <= 1.0
    assert ideal_score > messy_score


def test_basic_off_target_penalty_increases_with_mismatches_and_bulges():
    base = basic_off_target_penalty(0, 0)
    more_mismatch = basic_off_target_penalty(3, 0)
    more_bulge = basic_off_target_penalty(0, 2)
    both = basic_off_target_penalty(3, 2)

    assert 0.0 <= base <= 1.0
    assert base < more_mismatch <= both
    assert base < more_bulge <= both

