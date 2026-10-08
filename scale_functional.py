"""Declared scale functional. 0.08 is not an input.

I = Perimeter / sqrt(Area)
beta_I = dI / d ln w

A least-squares Seeley proxy is computed only to record that a0 < 0,
which rejects it as an extraction.
"""
from __future__ import annotations

import numpy as np
from numpy.linalg import eigvalsh


def dumbbell_mask(w: int, chamber: int, neck_len: int) -> np.ndarray:
    H = chamber
    W = chamber + neck_len + chamber
    mask = np.zeros((H, W), dtype=bool)
    mask[:, :chamber] = True
    mask[:, chamber + neck_len :] = True
    y0 = (H - w) // 2
    mask[y0 : y0 + w, chamber : chamber + neck_len] = True
    return mask


def perimeter(mask: np.ndarray) -> int:
    H, W = mask.shape
    p = 0
    for r, c in np.argwhere(mask):
        for dr, dc in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            rr, cc = r + dr, c + dc
            if not (0 <= rr < H and 0 <= cc < W and mask[rr, cc]):
                p += 1
    return p


def scale_I(w: int, chamber: int, neck_len: int) -> float:
    """Declared functional I = Perimeter / sqrt(Area) for one dumbbell."""
    mask = dumbbell_mask(w, chamber, neck_len)
    return perimeter(mask) / float(np.sqrt(mask.sum()))


def family_rows(chamber: int, neck_len: int) -> list[dict]:
    """I and the finite-difference beta_I for w = 1..chamber."""
    rows: list[dict] = []
    prev = None
    for w in range(1, chamber + 1):
        mask = dumbbell_mask(w, chamber, neck_len)
        area = int(mask.sum())
        peri = perimeter(mask)
        I = peri / float(np.sqrt(area))
        beta = None if prev is None else (I - prev[0]) / float(np.log(w / prev[1]))
        rows.append({"w": w, "area": area, "perimeter": peri, "I": I, "beta_I": beta})
        prev = (I, w)
    return rows


def beta_sign_changes(rows: list[dict]) -> list[int]:
    """Widths w where beta_I is zero or changes sign from the previous step."""
    hits = []
    betas = [(r["w"], r["beta_I"]) for r in rows if r["beta_I"] is not None]
    for i, (w, b) in enumerate(betas):
        if b == 0 or (i > 0 and np.sign(b) != np.sign(betas[i - 1][1])):
            hits.append(w)
    return hits


def grid_laplacian(mask: np.ndarray, bc: str = "neumann") -> np.ndarray:
    """5-point grid Laplacian on the mask. Dirichlet adds 1 per missing neighbour."""
    if bc not in ("neumann", "dirichlet"):
        raise ValueError("bc must be 'neumann' or 'dirichlet'")
    H, W = mask.shape
    pts = np.argwhere(mask)
    idx = -np.ones(mask.shape, dtype=int)
    idx[tuple(pts.T)] = np.arange(len(pts))
    L = np.zeros((len(pts), len(pts)))
    for k, (r, c) in enumerate(pts):
        for dr, dc in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            rr, cc = r + dr, c + dc
            if 0 <= rr < H and 0 <= cc < W and mask[rr, cc]:
                L[k, k] += 1.0
                L[k, idx[rr, cc]] -= 1.0
            elif bc == "dirichlet":
                L[k, k] += 1.0
    return L


def seeley_proxy(
    w: int, chamber: int, neck_len: int, bc: str = "neumann",
    t_min: float = 0.04, t_max: float = 0.20, n_t: int = 17,
) -> dict:
    """Least-squares fit K(t) ~ a0/t + a_half/sqrt(t) + a1 on [t_min, t_max].

    Recorded only to show the rejection: the fitted a0 is negative, so this
    is not a Seeley extraction and is never used as I.
    """
    ev = eigvalsh(grid_laplacian(dumbbell_mask(w, chamber, neck_len), bc))
    t = np.linspace(t_min, t_max, n_t)
    K = np.exp(-np.outer(t, ev)).sum(axis=1)
    A = np.column_stack([1.0 / t, 1.0 / np.sqrt(t), np.ones_like(t)])
    coef, *_ = np.linalg.lstsq(A, K, rcond=None)
    rel = float(np.linalg.norm(A @ coef - K) / np.linalg.norm(K))
    return {
        "w": w, "chamber": chamber, "neck": neck_len, "bc": bc,
        "a0": float(coef[0]), "a_half": float(coef[1]), "a1": float(coef[2]),
        "relative_residual": rel, "rejected": bool(coef[0] < 0),
    }


DEFAULT_FAMILIES = ((8, 5), (10, 6))


def main() -> None:
    print("I = Perimeter/sqrt(Area); 0.08 is not an input")
    for chamber, neck in DEFAULT_FAMILIES:
        print(f"\nchamber={chamber} neck={neck}")
        for r in family_rows(chamber, neck):
            beta = float("nan") if r["beta_I"] is None else r["beta_I"]
            print(f"w={r['w']:2d} area={r['area']:4d} peri={r['perimeter']:3d} I={r['I']:.4f} dI/dlnw={beta:.4f}")
        assert all(
            scale_I(w, chamber, neck) > scale_I(w + 1, chamber, neck)
            for w in range(1, chamber)
        )
    print("\nno stationary point: I strictly decreases with w on these families")


def _check_family(chamber: int, neck: int) -> None:
    if chamber < 2 or neck < 1:
        raise ValueError("need chamber >= 2 and neck >= 1")


def cli(argv: list[str] | None = None) -> int:
    """Command line: no arguments reproduces the RESULT tables; see --help."""
    import argparse
    import json
    import sys

    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        main()
        return 0
    parser = argparse.ArgumentParser(
        prog="scale-functional",
        description="Declared scale functional I = Perimeter/sqrt(Area) on lattice dumbbells. 0.08 is not an input.",
    )
    sub = parser.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("table", help="I and beta_I for w = 1..chamber on one family")
    t.add_argument("--chamber", type=int, default=8)
    t.add_argument("--neck", type=int, default=5)
    t.add_argument("--json", action="store_true")
    p = sub.add_parser("proxy", help="rejected least-squares Seeley proxy (records a0 < 0)")
    p.add_argument("--chamber", type=int, default=8)
    p.add_argument("--neck", type=int, default=5)
    p.add_argument("--w", type=int, action="append", help="neck width (repeatable; default 1, mid, chamber)")
    p.add_argument("--bc", choices=("neumann", "dirichlet", "both"), default="both")
    p.add_argument("--json", action="store_true")
    args = parser.parse_args(argv)
    try:
        _check_family(args.chamber, args.neck)
        if args.cmd == "table":
            rows = family_rows(args.chamber, args.neck)
            hits = beta_sign_changes(rows)
            out = {"chamber": args.chamber, "neck": args.neck, "rows": rows,
                   "beta_zero_or_sign_change_at_w": hits, "stationary_point_found": bool(hits)}
            if args.json:
                print(json.dumps(out, indent=2))
            else:
                print(f"chamber={args.chamber} neck={args.neck}; I = Perimeter/sqrt(Area); 0.08 is not an input")
                for r in rows:
                    beta = "n/a" if r["beta_I"] is None else f"{r['beta_I']:.4f}"
                    print(f"w={r['w']:2d} area={r['area']:4d} peri={r['perimeter']:3d} I={r['I']:.4f} beta_I={beta}")
                print("beta_I zero or sign change at w=" + ",".join(map(str, hits)) if hits
                      else "no zero or sign change of beta_I on this family")
            return 0
        widths = args.w or sorted({1, max(1, args.chamber // 2), args.chamber})
        for w in widths:
            if not 1 <= w <= args.chamber:
                raise ValueError(f"--w must be between 1 and chamber ({args.chamber})")
        bcs = ("neumann", "dirichlet") if args.bc == "both" else (args.bc,)
        fits = [seeley_proxy(w, args.chamber, args.neck, bc) for bc in bcs for w in widths]
        if args.json:
            print(json.dumps(fits, indent=2))
        else:
            print("K(t) ~ a0/t + a_half/sqrt(t) + a1 fitted on t in [0.04, 0.20]; a0 < 0 rejects it")
            for f in fits:
                print(f"{f['bc']:9s} w={f['w']:2d} a0={f['a0']:8.3f} a_half={f['a_half']:8.3f} "
                      f"a1={f['a1']:8.3f} rel_res={f['relative_residual']:.2e} "
                      f"{'REJECTED' if f['rejected'] else 'not rejected'}")
        return 0
    except ValueError as exc:
        print(f"scale-functional: error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(cli())
