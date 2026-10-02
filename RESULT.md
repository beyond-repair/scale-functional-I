# Result — 2026-10-02

**Classification:** RESEARCH  
**experimental_validation:** false

I = Perimeter / sqrt(Area). beta_I = dI / d ln w. Target 0.08 not used.

## Neumann, chamber 8, neck length 5

| w | Area | Perimeter | I | beta_I |
|---|------|-----------|---|--------|
| 1 | 133 | 72 | 6.2432 | n/a (no previous w) |
| 2 | 138 | 70 | 5.9588 | -0.4103 |
| 3 | 143 | 68 | 5.6864 | -0.6717 |
| 4 | 148 | 66 | 5.4252 | -0.9082 |
| 5 | 153 | 64 | 5.1741 | -1.1252 |
| 6 | 158 | 62 | 4.9325 | -1.3253 |
| 7 | 163 | 60 | 4.6996 | -1.5108 |
| 8 | 168 | 58 | 4.4748 | -1.6832 |

Table beta_I values are the finite difference `(I(w)-I(w-1))/ln(w/(w-1))` printed by `scale_functional.py` on head `23c00dd701f7c2178ccc4bc5d1a7ee58075d3656`. The previous column (about -0.62 to -0.86) did not match that definition. Sign and absence of a zero are unchanged.

beta_I zeros: none.

## Dirichlet, same geometry

I and beta_I match the Neumann row. The declared functional is geometric. Boundary condition does not enter I. It does enter the spectrum, which is why the earlier spectral ratios moved and this ratio does not.

## Neumann, chamber 10, neck length 6

I falls from 6.271 at w=1 to 4.465 at w=10. beta_I stays negative. No zero.

## Comparison step

Not reached. There is no ell_* with beta_I(ell_*)=0 on these families, so no dimensionless stationary value is compared with 0.08.

## Lattice proxy

Rejected. Fitted a0 was negative (about -5 to -10) with relative residual about 0.006. Small residual does not repair the sign.

## Conclusion

This independently defined I has no universal fixed point on the tested pinch families. That is outcome 3 for this functional. It does not revive 0.08, and it does not close every future F.
