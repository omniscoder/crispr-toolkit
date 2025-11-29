# crispr-toolkit

A small, pure-Python utility library for CRISPR and prime-editing sequence design, validation, and basic scoring. This package is intended as a clean, well-tested FOSS backend for browser tools and scripting workflows.

> This project is fully independent from Helix Studio internals and is provided as a standalone, open-source utility library.

## Installation

```bash
python -m pip install crispr-toolkit
```

## Quickstart

```python
from crispr_toolkit import normalize_seq, find_pam_sites, validate_guide

seq = normalize_seq("acgtacgtGGGcttacggt")
sites = find_pam_sites(seq, "NGG")
result = validate_guide(seq[2:24])

print(sites)
print(result.is_valid, result.gc_content)
```

## Roadmap

Planned future features include:

- More PAM presets and convenience helpers.
- Simple off-target enumeration utilities.
- Additional rule-based prime editing recipes and presets.

## License

This project is released under the MIT License. See `LICENSE` for details.

