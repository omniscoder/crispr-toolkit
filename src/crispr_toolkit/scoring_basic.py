"""
Simple, transparent scoring helpers.

These functions are intentionally basic and heuristic. They are *not*
intended to represent production-grade on/off-target models, only to give
users a rough, explainable signal.
"""

from __future__ import annotations

from .seq import gc_content, has_homopolymer, normalize_seq


def basic_on_target_score(seq: str) -> float:
    """
    Compute a crude on-target score in [0.0, 1.0].

    Penalizes:
    - Length deviations from 20 nt.
    - GC content far from 50%.
    - Presence of homopolymers (>=4, >=5).
    """
    norm = normalize_seq(seq)
    if not norm:
        return 0.0

    score = 1.0

    # Length penalty: prefer 20 nt.
    ideal_len = 20
    length_diff = abs(len(norm) - ideal_len)
    score -= min(length_diff * 0.02, 0.3)

    # GC penalty: prefer ~50% GC.
    gc = gc_content(norm)
    gc_diff = abs(gc - 0.5)
    score -= min(gc_diff * 0.6, 0.3)

    # Homopolymer penalties.
    if has_homopolymer(norm, 4):
        score -= 0.2
    if has_homopolymer(norm, 5):
        score -= 0.1

    return max(0.0, min(1.0, score))


def basic_off_target_penalty(mismatches: int, bulges: int) -> float:
    """
    Compute a crude off-target penalty in [0.0, 1.0].

    Higher values indicate more severe divergence from the intended target.
    Intended usage is to subtract this from an on-target score, or to use it
    qualitatively when ranking candidates.
    """
    if mismatches < 0 or bulges < 0:
        raise ValueError("mismatches and bulges must be non-negative")

    penalty = mismatches * 0.08 + bulges * 0.12
    return max(0.0, min(1.0, penalty))

