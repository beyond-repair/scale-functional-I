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

## Run it

Needs Python 3.9+ and NumPy.

```bash
git clone https://github.com/beyond-repair/scale-functional-I.git
cd scale-functional-I
python3 -m venv .venv && . .venv/bin/activate
pip install -e ".[test]"

python3 scale_functional.py          # reproduces the RESULT.md tables, exit 0
scale-functional table --chamber 12 --neck 7          # any dumbbell family
scale-functional table --chamber 10 --neck 6 --json   # machine-readable rows
scale-functional proxy               # the rejected Seeley proxy fit, Neumann and Dirichlet
python3 -m pytest                    # 13 tests
```

Without installing, `pip install numpy` and run `python3 scale_functional.py` (no arguments) or `python3 scale_functional.py table ...` from the checkout.

The default run ends with `no stationary point: I strictly decreases with w on these families`. `table` reports whether beta_I has a zero or sign change on the family you ask for; it does not assert the answer. `proxy` refits `K(t) ~ a0/t + a_half/sqrt(t) + a1` on 17 evenly spaced t in [0.04, 0.20] with the 5-point grid Laplacian and prints `REJECTED` when a0 < 0. The original fit code was not committed; this reconstruction reproduces the recorded range (a0 about -5 to -10, relative residual about 0.006). Bad input exits 2.
