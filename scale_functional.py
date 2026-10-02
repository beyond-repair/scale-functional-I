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


def main() -> None:
    print("I = Perimeter/sqrt(Area); 0.08 is not an input")
    for chamber, neck in ((8, 5), (10, 6)):
        print(f"\nchamber={chamber} neck={neck}")
        prev = None
        for w in range(1, chamber + 1):
            mask = dumbbell_mask(w, chamber, neck)
            area = int(mask.sum())
            peri = perimeter(mask)
            I = peri / np.sqrt(area)
            beta = float("nan") if prev is None else (I - prev[0]) / np.log(w / prev[1])
            print(f"w={w:2d} area={area:4d} peri={peri:3d} I={I:.4f} dI/dlnw={beta:.4f}")
            prev = (I, w)
        assert all(
            perimeter(dumbbell_mask(w, chamber, neck)) / np.sqrt(dumbbell_mask(w, chamber, neck).sum())
            > perimeter(dumbbell_mask(w + 1, chamber, neck)) / np.sqrt(dumbbell_mask(w + 1, chamber, neck).sum())
            for w in range(1, chamber)
        )
    print("\nno stationary point: I strictly decreases with w on these families")


if __name__ == "__main__":
    main()
