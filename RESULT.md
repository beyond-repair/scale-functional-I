# Result — 2026-10-02

**Classification:** RESEARCH  
**experimental_validation:** false

I = Perimeter / sqrt(Area). beta_I = dI / d ln w. Target 0.08 not used.

## Neumann, chamber 8, neck length 5

| w | Area | Perimeter | I | beta_I |
|---|------|-----------|---|--------|
| 1 | 133 | 72 | 6.243 | -0.62 |
| 2 | 138 | 70 | 5.959 | -0.65 |
| 3 | 143 | 68 | 5.686 | -0.69 |
| 4 | 148 | 66 | 5.425 | -0.72 |
| 5 | 153 | 64 | 5.174 | -0.76 |
| 6 | 158 | 62 | 4.932 | -0.80 |
| 7 | 163 | 60 | 4.700 | -0.84 |
| 8 | 168 | 58 | 4.475 | -0.86 |

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
