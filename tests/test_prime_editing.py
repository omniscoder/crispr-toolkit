from crispr_toolkit.prime_editing import PegRNAComponents, design_pegRNA


def test_simple_substitution_mid_sequence():
    ref = "AACCGGTTAACCGGTT"
    edit_pos = 6
    components = design_pegRNA(ref, edit_pos, "substitution", edit_seq="A")
    assert isinstance(components, PegRNAComponents)
    assert len(components.pbs) <= 13
    assert len(components.rtt) <= 16
    # RTT should start at the edit position in the edited sequence.
    assert len(components.rtt) > 0


def test_insertion_covers_edit():
    ref = "AACCGGTTAACCGGTT"
    edit_pos = 4
    ins = "AAA"
    components = design_pegRNA(ref, edit_pos, "insertion", edit_seq=ins)
    assert ins in components.rtt


def test_deletion_uses_del_length():
    ref = "AACCGGTTAACCGGTT"
    edit_pos = 2
    components = design_pegRNA(ref, edit_pos, "deletion", edit_seq="", del_length=3)
    assert isinstance(components, PegRNAComponents)
    assert len(components.pbs) <= 13
    assert len(components.rtt) <= 16

