"""
PAM pattern helpers.

Provides simple scanning of sequences for PAM-like motifs with support for
IUPAC degenerate bases (N, R, Y, etc.).
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Tuple

from .seq import normalize_seq, reverse_complement


_IUPAC_CODES: Dict[str, Iterable[str]] = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}


def _match_base(pattern_base: str, sequence_base: str) -> bool:
    allowed = _IUPAC_CODES.get(pattern_base, pattern_base)
    return sequence_base in allowed


def match_pam(pam_pattern: str, window: str) -> bool:
    """
    Return True if a PAM pattern matches the provided window.

    Both the pattern and window are normalized (uppercase, U->T). The pattern
    may contain IUPAC degenerate bases such as N, R, Y, etc.
    """
    pat = normalize_seq(pam_pattern)
    win = normalize_seq(window)
    if len(pat) != len(win):
        return False
    return all(_match_base(pat_b, win_b) for pat_b, win_b in zip(pat, win))


def find_pam_sites(seq: str, pam_pattern: str, strand: str = "both") -> List[Tuple[int, str]]:
    """
    Scan a sequence for PAM occurrences.

    Parameters
    ----------
    seq:
        Genomic sequence to scan (any case, optional whitespace).
    pam_pattern:
        PAM motif using IUPAC codes (e.g. ``\"NGG\"`` for SpCas9).
    strand:
        ``\"plus\"``, ``\"minus\"``, or ``\"both\"``. For the minus strand,
        the reverse-complement of ``pam_pattern`` is used, but positions are
        always reported as 0-based indices on the plus strand corresponding
        to the PAM start.

    Returns
    -------
    list[tuple[int, str]]
        List of (position, strand) tuples where strand is ``\"+\"`` or ``\"-\"``.
    """
    sequence = normalize_seq(seq)
    pattern = normalize_seq(pam_pattern)
    k = len(pattern)

    if k == 0:
        return []

    strand = strand.lower()
    results: List[Tuple[int, str]] = []

    scan_plus = strand in {"plus", "+", "both", "both_strands"}
    scan_minus = strand in {"minus", "-", "both", "both_strands"}

    if scan_plus:
        for i in range(0, len(sequence) - k + 1):
            window = sequence[i : i + k]
            if match_pam(pattern, window):
                results.append((i, "+"))

    if scan_minus:
        pattern_rc = reverse_complement(pattern)
        for i in range(0, len(sequence) - k + 1):
            window = sequence[i : i + k]
            if match_pam(pattern_rc, window):
                results.append((i, "-"))

    return results

