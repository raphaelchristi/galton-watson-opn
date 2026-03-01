# The Cyclotomic Factor Cascade as a Supercritical Branching Process

**Computational companion to the paper by Raphael Christi (2026)**

## Overview

We model the factor chain induced by the cyclotomic polynomial Φ₃(p) = p² + p + 1 as a Galton–Watson branching process. Over 300 seed configurations generating 12,927 parent–child observations, we find a global mean offspring rate **μ̄ = 1.1747** (supercritical), with a phase transition from subcritical to supercritical near p ≈ 200.

This provides a quantitative explanation for why no finite set of odd primes is closed under the σ-induced factor map — probabilistic evidence against the existence of odd perfect numbers.

## Files

| File | Description |
|------|-------------|
| `galton_watson.py` | Main computation script (~0.5s on Mac Mini M4) |
| `paper.tex` | LaTeX source of the paper |
| `paper.pdf` | Compiled paper |

## Quick Start

```bash
python3 galton_watson.py
```

**Requirements:** Python 3.8+ (standard library only, no dependencies).

**Output:**
- Terminal: full statistics (μ̄, offspring distribution, μ(p) by prime size, residue classes, Stewart ratio, multi-Φ comparison)
- `galton_watson_plot.html`: interactive Chart.js plots (open in browser)

## Key Results

| Measurement | Value |
|---|---|
| Global μ̄ | **1.1747** (supercritical) |
| Phase transition | p ≈ 200 |
| μ for p ≡ 1 (mod 3) | 1.182 |
| μ for Φ₅ | 2.987 |
| μ for Φ₇ | 1.707 |
| Digit ratio (median) | 1.345 (Stewart predicts → 2.0) |
| Cascades closed | **0 / 300** |

## Citation

```bibtex
@article{christi2026galtonwatson,
  title={The Cyclotomic Factor Cascade as a Supercritical Branching Process},
  author={Christi, Raphael},
  year={2026},
  note={Preprint}
}
```

## License

MIT
