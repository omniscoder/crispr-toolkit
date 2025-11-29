"""
Guide RNA validation helpers.

Simple heuristics for turning a raw sequence into a \"sane guide?\" verdict.
No ML models, no proprietary scoring – just transparent rule-based checks.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List

from .seq import gc_content, has_homopolymer, normalize_seq


@dataclass
class GuideCheckResult:
    is_valid: bool
    reasons: List[str]
    length: int
    gc_content: float


def validate_guide(
    seq: str,
    *,
    min_len: int = 19,
    max_len: int = 23,
    min_gc: float = 0.30,
    max_gc: float = 0.80,
    max_homopolymer: int = 4,
) -> GuideCheckResult:
    """
    Validate a candidate CRISPR guide sequence using simple heuristics.

    Checks performed:
    - Sequence normalization to A/C/G/T only (U->T).
    - Only A/C/G/T bases are allowed.
    - Length must lie within [min_len, max_len].
    - GC content must lie within [min_gc, max_gc].
    - No homopolymer stretches of length >= max_homopolymer.
    """
    norm = normalize_seq(seq)
    reasons: List[str] = []

    if not norm:
        reasons.append("sequence is empty after normalization")

    valid_bases = {"A", "C", "G", "T"}
    invalid_bases = sorted({b for b in norm if b not in valid_bases})
    if invalid_bases:
        reasons.append(f"contains invalid bases: {''.join(invalid_bases)}")

    length = len(norm)
    if length < min_len:
        reasons.append(f"length {length} < minimum {min_len}")
    if length > max_len:
        reasons.append(f"length {length} > maximum {max_len}")

    gc = gc_content(norm)
    if gc < min_gc:
        reasons.append(f"GC content {gc:.3f} < minimum {min_gc:.3f}")
    if gc > max_gc:
        reasons.append(f"GC content {gc:.3f} > maximum {max_gc:.3f}")

    if has_homopolymer(norm, max_homopolymer):
        reasons.append(f"contains homopolymer of length >= {max_homopolymer}")

    is_valid = not reasons
    return GuideCheckResult(
        is_valid=is_valid,
        reasons=reasons,
        length=length,
        gc_content=gc,
    )

