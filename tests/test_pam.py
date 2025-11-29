from crispr_toolkit.pam import find_pam_sites, match_pam


def test_match_pam_simple():
    assert match_pam("NGG", "AGG")
    assert not match_pam("NGG", "AAA")


def test_match_pam_degenerate():
    # R = A/G, Y = C/T, N = any
    assert match_pam("RYN", "ACT")
    assert match_pam("RYN", "GCC")
    assert not match_pam("RYN", "TTT")


def test_find_pam_spcas9_plus_strand():
    seq = "AAGTCCTAGGCTTGGG"
    hits = find_pam_sites(seq, "NGG", strand="plus")
    # Expect NGG at positions corresponding to "AGG", "TGG", and "GGG".
    positions = {pos for (pos, s) in hits}
    assert positions == {7, 12, 13}
    assert all(s == "+" for (_, s) in hits)


def test_find_pam_both_strands_ngg():
    # On the minus strand, NGG corresponds to CCN on the plus sequence.
    seq = "CCCAAGGTT"
    hits = find_pam_sites(seq, "NGG", strand="both")
    plus_hits = {(p, s) for (p, s) in hits if s == "+"}
    minus_hits = {(p, s) for (p, s) in hits if s == "-"}

    # Plus strand NGG: "AGG" starting at position 4.
    assert (4, "+") in plus_hits

    # Minus strand NGG: look for CCN (pattern RC) starting at position 0.
    assert (0, "-") in minus_hits
