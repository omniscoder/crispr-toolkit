"""
Prime editing scaffold design (minimal, rule-based).

This module intentionally provides only a very simple pegRNA component
builder that is easy to reason about and does not attempt to mimic any
proprietary design pipeline.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List

from .seq import normalize_seq


@dataclass
class PegRNAComponents:
    spacer: str
    pbs: str
    rtt: str
    notes: List[str]


def design_pegRNA(
    ref_seq: str,
    edit_pos: int,
    edit_type: str,
    edit_seq: str = "",
    del_length: int = 1,
    *,
    pbs_length: int = 13,
    rtt_length: int = 16,
) -> PegRNAComponents:
    """
    Design a very simple pegRNA template around an intended edit.

    Parameters
    ----------
    ref_seq:
        Reference DNA sequence (plus strand). Normalized internally.
    edit_pos:
        0-based index in the normalized reference sequence indicating the
        position of the edit.
    edit_type:
        One of ``\"substitution\"``, ``\"insertion\"``, or ``\"deletion\"``.
    edit_seq:
        For substitutions: the replacement bases.
        For insertions: the inserted bases.
        For deletions: typically empty; see ``del_length``.
    del_length:
        For substitutions: number of reference bases replaced (default 1).
        For deletions: number of bases removed.
        For insertions: ignored.
    pbs_length, rtt_length:
        Target lengths for the primer binding site (PBS) and reverse
        transcriptase template (RTT). If the sequence ends are reached, the
        actual lengths may be shorter; a note is added in that case.

    Notes
    -----
    This is a deliberately simple design helper. The spacer is chosen as a
    20-nt window around ``edit_pos`` where possible. The PBS is the region
    immediately upstream of the edit on the reference. The RTT is taken from
    the edited sequence starting at ``edit_pos``.
    """
    norm_ref = normalize_seq(ref_seq)
    n = len(norm_ref)
    notes: List[str] = []

    if not (0 <= edit_pos <= n):
        raise ValueError("edit_pos must be within the reference sequence bounds")

    etype = edit_type.lower()
    norm_edit = normalize_seq(edit_seq)

    if etype == "insertion":
        effective_del = 0
    elif etype == "substitution":
        effective_del = max(del_length, 1)
    elif etype == "deletion":
        if del_length <= 0:
            raise ValueError("del_length must be > 0 for deletions")
        effective_del = del_length
    else:
        raise ValueError("edit_type must be 'substitution', 'insertion', or 'deletion'")

    if edit_pos + effective_del > n:
        raise ValueError("edit spans beyond the end of the reference sequence")

    left = norm_ref[:edit_pos]
    right = norm_ref[edit_pos + effective_del :]

    edited = left + norm_edit + right

    # Spacer: try to center around the edit with length 20 where possible.
    spacer_len = 20
    spacer_start = max(0, min(edit_pos - spacer_len // 2, n - spacer_len))
    spacer_end = max(spacer_start, min(spacer_start + spacer_len, n))
    spacer = norm_ref[spacer_start:spacer_end]
    if len(spacer) < spacer_len:
        notes.append("spacer truncated at sequence boundary")

    # PBS: upstream of the edit site on the reference.
    pbs_start = max(0, edit_pos - pbs_length)
    pbs = norm_ref[pbs_start:edit_pos]
    if len(pbs) < pbs_length:
        notes.append("PBS truncated at sequence boundary")

    # RTT: downstream of the edit in the edited sequence.
    rtt = edited[edit_pos : edit_pos + rtt_length]
    if len(rtt) < rtt_length:
        notes.append("RTT truncated at sequence boundary")

    return PegRNAComponents(spacer=spacer, pbs=pbs, rtt=rtt, notes=notes)

