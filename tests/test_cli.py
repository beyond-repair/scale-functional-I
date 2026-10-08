"""CLI and library checks added with the runnable packaging (2026-10-08)."""
from __future__ import annotations

import io
import json
import unittest
from contextlib import redirect_stderr, redirect_stdout

from scale_functional import (
    beta_sign_changes,
    cli,
    dumbbell_mask,
    family_rows,
    grid_laplacian,
    perimeter,
    scale_I,
    seeley_proxy,
)


def run(argv):
    out, err = io.StringIO(), io.StringIO()
    with redirect_stdout(out), redirect_stderr(err):
        code = cli(argv)
    return code, out.getvalue(), err.getvalue()


class GeometryTests(unittest.TestCase):
    def test_result_table_values(self) -> None:
        rows = family_rows(8, 5)
        self.assertEqual([r["area"] for r in rows], [133, 138, 143, 148, 153, 158, 163, 168])
        self.assertEqual([r["perimeter"] for r in rows], [72, 70, 68, 66, 64, 62, 60, 58])
        self.assertAlmostEqual(rows[0]["I"], 6.2432, places=4)
        self.assertAlmostEqual(rows[-1]["beta_I"], -1.6832, places=4)

    def test_chamber_10_endpoints(self) -> None:
        self.assertAlmostEqual(scale_I(1, 10, 6), 6.271, places=3)
        self.assertAlmostEqual(scale_I(10, 10, 6), 4.465, places=3)

    def test_no_sign_change_on_declared_families(self) -> None:
        for chamber, neck in ((8, 5), (10, 6)):
            self.assertEqual(beta_sign_changes(family_rows(chamber, neck)), [])

    def test_sign_change_detector(self) -> None:
        rows = [{"w": 1, "beta_I": None}, {"w": 2, "beta_I": -1.0}, {"w": 3, "beta_I": 0.5}]
        self.assertEqual(beta_sign_changes(rows), [3])

    def test_perimeter_of_full_rectangle(self) -> None:
        mask = dumbbell_mask(2, 2, 1)
        self.assertEqual(int(mask.sum()), 10)  # full 2x5 rectangle
        self.assertEqual(perimeter(mask), 14)


class ProxyTests(unittest.TestCase):
    def test_laplacian_neumann_rows_sum_to_zero(self) -> None:
        L = grid_laplacian(dumbbell_mask(3, 8, 5), "neumann")
        self.assertAlmostEqual(float(abs(L.sum(axis=1)).max()), 0.0)

    def test_proxy_a0_negative_both_bcs(self) -> None:
        for bc in ("neumann", "dirichlet"):
            for chamber, neck in ((8, 5), (10, 6)):
                for w in (1, chamber):
                    fit = seeley_proxy(w, chamber, neck, bc)
                    self.assertLess(fit["a0"], 0)
                    self.assertTrue(fit["rejected"])
                    self.assertLess(fit["relative_residual"], 0.01)
                    self.assertTrue(-11 < fit["a0"] < -5)

    def test_bad_bc(self) -> None:
        with self.assertRaises(ValueError):
            grid_laplacian(dumbbell_mask(1, 4, 2), "robin")


class CliTests(unittest.TestCase):
    def test_no_args_reproduces_result(self) -> None:
        code, out, _ = run([])
        self.assertEqual(code, 0)
        self.assertIn("no stationary point", out)
        self.assertIn("w= 8 area= 168 peri= 58 I=4.4748 dI/dlnw=-1.6832", out)

    def test_table_json(self) -> None:
        code, out, _ = run(["table", "--chamber", "10", "--neck", "6", "--json"])
        self.assertEqual(code, 0)
        data = json.loads(out)
        self.assertEqual(len(data["rows"]), 10)
        self.assertFalse(data["stationary_point_found"])

    def test_proxy_json(self) -> None:
        code, out, _ = run(["proxy", "--w", "2", "--bc", "dirichlet", "--json"])
        self.assertEqual(code, 0)
        fits = json.loads(out)
        self.assertEqual(len(fits), 1)
        self.assertTrue(fits[0]["rejected"])

    def test_bad_input_exit_2(self) -> None:
        self.assertEqual(run(["proxy", "--w", "99"])[0], 2)
        self.assertEqual(run(["table", "--chamber", "1"])[0], 2)


if __name__ == "__main__":
    unittest.main()
