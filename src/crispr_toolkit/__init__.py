"""
crispr_toolkit
==============

A small, pure-Python utility library for CRISPR and prime-editing sequence
design, validation, and simple heuristic scoring.

This package is independent from Helix Studio internals and is intended as
general-purpose FOSS infrastructure for browser tools and scripting.
"""

from .seq import gc_content, normalize_seq, reverse_complement
from .pam import find_pam_sites
from .grna import GuideCheckResult, validate_guide
from .prime_editing import PegRNAComponents, design_pegRNA
from .scoring_basic import basic_off_target_penalty, basic_on_target_score

__all__ = [
    "normalize_seq",
    "reverse_complement",
    "gc_content",
    "find_pam_sites",
    "GuideCheckResult",
    "validate_guide",
    "PegRNAComponents",
    "design_pegRNA",
    "basic_on_target_score",
    "basic_off_target_penalty",
]

__version__ = "0.1.0"
