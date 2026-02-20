import os
from pathlib import Path
import unittest


def _import_root_or_none():
    try:
        import ROOT  # type: ignore

        return ROOT
    except Exception:
        return None


class TestMCTrialVariation(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.ROOT = _import_root_or_none()

    def _max_abs_diff(self, h_default, h_trial) -> float:
        nbins = int(h_default.GetNbinsX()) + 2
        max_diff = 0.0
        for ibin in range(nbins):
            max_diff = max(max_diff, abs(float(h_default.GetBinContent(ibin)) - float(h_trial.GetBinContent(ibin))))
        return max_diff

    def test_trial_variations_change_mc_histograms(self) -> None:
        if self.ROOT is None:
            self.skipTest("ROOT not available")

        nuclei_input = os.environ.get("NUCLEI_INPUT", "")
        if not nuclei_input:
            self.skipTest("NUCLEI_INPUT is not set")

        ref_path = Path(nuclei_input) / "smoke_references" / "python_single" / "LHC24_apass1__LHC25b9" / "MChistos_he3_ref.root"
        if not ref_path.is_file():
            self.skipTest(f"Missing MC reference file: {ref_path}")

        root_file = self.ROOT.TFile(str(ref_path))
        self.assertTrue(root_file and not root_file.IsZombie(), f"Cannot open ROOT file: {ref_path}")

        default_hist_names = ["TPCAHe3", "TOFAHe3", "TPCMHe3", "TOFMHe3"]
        trial_dirs = []
        for key in root_file.GetListOfKeys():
            name = str(key.GetName())
            if name.startswith("nuclei") and name[6:].isdigit():
                trial_dirs.append(name)
        self.assertTrue(trial_dirs, "No MC trial directories found (nuclei0, nuclei1, ...).")

        changed = False
        for hist_name in default_hist_names:
            h_default = root_file.Get(f"nuclei/{hist_name}")
            if not h_default:
                continue
            for trial_dir in trial_dirs:
                h_trial = root_file.Get(f"{trial_dir}/{hist_name}")
                if not h_trial:
                    continue
                if self._max_abs_diff(h_default, h_trial) > 0.0:
                    changed = True
                    break
            if changed:
                break

        root_file.Close()
        self.assertTrue(changed, "Expected at least one MC trial histogram to differ from default output.")


if __name__ == "__main__":
    unittest.main()
