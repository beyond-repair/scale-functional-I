# Claim status — scale-functional-I

**Classification:** RESEARCH  
**Claim level:** 1 (declared functional, computed on lattice dumbbells)  
**experimental_validation:** false  
**thrust_validated:** false  
**energy_extraction_validated:** false

I = Perimeter / sqrt(Area) is declared. beta_I has no zero on the tested dumbbells. Comparison with 0.08 is not reached.

Does not modify the 2026-10-01 lock. Parent: [-ware-constant-derivation](https://github.com/beyond-repair/-ware-constant-derivation).

## Tooling note — 2026-10-08

Added packaging (`pyproject.toml`, NumPy declared), a `scale-functional` command (`table`, `proxy`), and `tests/test_cli.py`. `python3 scale_functional.py` with no arguments prints the same output as before. `proxy` restores the rejected least-squares fit described in FUNCTIONAL.md; it is a reconstruction and reproduces a0 < 0 on both boundary conditions. Claim level unchanged.
