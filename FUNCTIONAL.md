# Functional declaration

**Date:** 2026-10-02  
**Assumption class:** A4 model definition, fixed before the numerical comparison.  
**Forbidden input:** 0.08, beta = -0.005888, 1/6, 1/12 as targets.

## Definition

Let X_ell be a planar dumbbell with neck width ell, fixed chamber size, fixed neck length. Let L_ell be the grid Laplacian (Dirichlet or Neumann).

Continuum d=2 heat trace:

```text
K(t) ~ a0 t^{-1} + a_{1/2} t^{-1/2} + a1 + ...
a0 proportional to Area
a_{1/2} proportional to Perimeter
```

The universal factors cancel in the shape ratio, so the declared functional is

```text
I(ell) = Perimeter(ell) / sqrt(Area(ell))
beta_I(ell) = dI / d ln ell
```

Stationary test: beta_I(ell_*) = 0. Only after a zero survives neck length, chamber size, boundary condition, and mesh refinement is I(ell_*) compared with anything.

## Rejected extraction

Fitting K ~ a0/t + a_half/sqrt(t) + a1 on t in [0.04, 0.20] gave a0 < 0 on every tested dumbbell. A negative leading coefficient is not a0. That proxy is not used as I.

## What this does not claim

This I is the shape content of the first two coefficients. It is not the unique scale functional. A later F(a0, a_{1/2}, a1, ...) may be declared the same way: write F first, then test. Do not choose F because its stationary value is near 0.08.
