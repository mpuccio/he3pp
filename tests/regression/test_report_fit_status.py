import unittest

from he3pp.report_logic import status_from_metrics


class TestReportFitStatus(unittest.TestCase):
    def test_missing_when_plot_not_available(self) -> None:
        self.assertEqual(status_from_metrics(False, {}, 0.05), ("MISSING", "missing"))

    def test_unknown_when_metrics_missing(self) -> None:
        self.assertEqual(status_from_metrics(True, {}, 0.05), ("UNK", "unknown"))

    def test_ok_when_pvalue_meets_threshold(self) -> None:
        self.assertEqual(status_from_metrics(True, {"p_value": 0.05}, 0.05), ("OK", "ok"))
        self.assertEqual(status_from_metrics(True, {"p_value": 0.5}, 0.05), ("OK", "ok"))

    def test_ko_when_pvalue_below_threshold(self) -> None:
        self.assertEqual(status_from_metrics(True, {"p_value": 0.049}, 0.05), ("KO", "ko"))


if __name__ == "__main__":
    unittest.main()
