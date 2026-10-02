<div align="center">

# Scale functional I

### Defined before any comparison with 0.08

[![RESEARCH](https://img.shields.io/badge/claim_%E2%89%A41-7c3aed?style=for-the-badge)](https://github.com/beyond-repair/ADL-Governance)

</div>

---

This repo holds the next object after the 2026-10-01 pinch falsification. It does not reopen that lock.

Parent lock: [STATUS_LOCK_2026-10-01.md](https://github.com/beyond-repair/-ware-constant-derivation/blob/main/STATUS_LOCK_2026-10-01.md)

## Declared functional

For a planar pinch domain, the first two Seeley coefficients are bulk volume and boundary length, up to universal factors that do not depend on the pinch parameter. The shape content is

```text
I(ell) = Perimeter(ell) / sqrt(Area(ell))
beta_I = dI / d ln ell
```

0.08 is not an input. 1/6 and 1/12 are not candidates.

A lattice least-squares proxy `K ~ a0/t + a_half/sqrt(t) + a1` was tried and rejected: the fitted leading coefficient was negative, so it is not a Seeley extraction.

## Result on the tested families

No zero of beta_I. I is monotone in neck width for Neumann and Dirichlet dumbbells at two chamber sizes. Stationary-point comparison with 0.08 does not arise.

```text
experimental_validation     = false
thrust_validated            = false
energy_extraction_validated = false
```

Run: `python3 scale_functional.py`
