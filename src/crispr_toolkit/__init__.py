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
from .scoring_basic import basic_on_target_score

__all__ = [
    "normalize_seq",
    "reverse_complement",
    "gc_content",
    "find_pam_sites",
    "GuideCheckResult",
    "validate_guide",
    "basic_on_target_score",
]

__version__ = "0.1.0"

