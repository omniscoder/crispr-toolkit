"""
Low-level sequence utilities for CRISPR workflows.

These helpers intentionally stay "boring" and predictable so they can be
reused in different pipelines and tools.
"""

from __future__ import annotations

from typing import Dict


_COMPLEMENT_TABLE: Dict[str, str] = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    # IUPAC ambiguity codes – keep symmetric where possible.
    "R": "Y",  # A/G <-> C/T
    "Y": "R",
    "S": "S",  # G/C
    "W": "W",  # A/T
    "K": "M",  # G/T <-> A/C
    "M": "K",
    "B": "V",  # not A <-> not T
    "D": "H",  # not C <-> not G
    "H": "D",
    "V": "B",
    "N": "N",
}


def normalize_seq(seq: str) -> str:
    """
    Normalize a nucleotide sequence string.

    - Strip whitespace.
    - Convert to uppercase.
    - Replace RNA U with DNA T.
    """
    # Remove all ASCII whitespace characters.
    cleaned = "".join(ch for ch in seq if not ch.isspace())
    cleaned = cleaned.upper().replace("U", "T")
    return cleaned


def reverse_complement(seq: str) -> str:
    """
    Return the reverse complement of a DNA sequence.

    Unknown characters are preserved but reversed.
    IUPAC ambiguity codes are complemented where possible.
    """
    norm = normalize_seq(seq)
    # For characters without a defined complement, fall back to themselves.
    return "".join(_COMPLEMENT_TABLE.get(base, base) for base in reversed(norm))


def gc_content(seq: str) -> float:
    """
    Compute GC content as a fraction in [0.0, 1.0].

    All characters contribute to the length; only G/C contribute to the GC
    count. Non-ACGT characters reduce the apparent GC fraction.
    """
    norm = normalize_seq(seq)
    if not norm:
        return 0.0
    gc = sum(1 for b in norm if b in {"G", "C"})
    return gc / len(norm)


def has_homopolymer(seq: str, length: int = 4) -> bool:
    """
    Return True if the sequence contains a homopolymer stretch of at least
    ``length`` identical bases.
    """
    if length <= 1:
        raise ValueError("length must be >= 2 for homopolymer detection")

    norm = normalize_seq(seq)
    if not norm:
        return False

    current_base = norm[0]
    run_length = 1
    for base in norm[1:]:
        if base == current_base:
            run_length += 1
            if run_length >= length:
                return True
        else:
            current_base = base
            run_length = 1
    return False

