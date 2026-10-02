"""Witness that declared I has no stationary point on the two dumbbell families."""
from __future__ import annotations

import unittest

import numpy as np

from scale_functional import dumbbell_mask, perimeter


def scale_I(w: int, chamber: int, neck: int) -> float:
    mask = dumbbell_mask(w, chamber, neck)
    return perimeter(mask) / float(np.sqrt(mask.sum()))


class ScaleFunctionalTests(unittest.TestCase):
    def test_strict_decrease_and_negative_beta(self) -> None:
        for chamber, neck in ((8, 5), (10, 6)):
            values = [scale_I(w, chamber, neck) for w in range(1, chamber + 1)]
            self.assertTrue(all(values[i] > values[i + 1] for i in range(len(values) - 1)))
            betas = [
                (values[i] - values[i - 1]) / np.log((i + 1) / i)
                for i in range(1, len(values))
            ]
            self.assertTrue(all(beta < 0 for beta in betas))


if __name__ == "__main__":
    unittest.main()
