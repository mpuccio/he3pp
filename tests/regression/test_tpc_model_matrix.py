import unittest

from he3pp.tasks_common import resolve_tpc_model_matrix


class TestTPCModelMatrix(unittest.TestCase):
    def test_resolve_matrix_default_order(self) -> None:
        matrix = resolve_tpc_model_matrix(
            ["GausGaus", "ExpGaus", "ExpTailGaus", "LognormalLognormal"],
            ["GausGaus", "ExpGaus", "ExpTailGaus", "LognormalLognormal"],
        )
        self.assertEqual(matrix, [("GausGaus", 0), ("ExpGaus", 1), ("ExpTailGaus", 2), ("LognormalLognormal", 3)])

    def test_resolve_matrix_reordered_subset(self) -> None:
        matrix = resolve_tpc_model_matrix(
            ["ExpGaus", "GausGaus"],
            ["GausGaus", "ExpGaus", "ExpTailGaus", "LognormalLognormal"],
        )
        self.assertEqual(matrix, [("ExpGaus", 1), ("GausGaus", 0)])

    def test_resolve_matrix_rejects_duplicates(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "Duplicated TPC function name"):
            resolve_tpc_model_matrix(
                ["ExpGaus", "ExpGaus"],
                ["GausGaus", "ExpGaus", "ExpTailGaus", "LognormalLognormal"],
            )

    def test_resolve_matrix_rejects_unknown_name(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "Unsupported TPC function name"):
            resolve_tpc_model_matrix(
                ["ExpGaus", "NotAModel"],
                ["GausGaus", "ExpGaus", "ExpTailGaus", "LognormalLognormal"],
            )


if __name__ == "__main__":
    unittest.main()
