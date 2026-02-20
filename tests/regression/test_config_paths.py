import unittest

from he3pp import settings as s


class TestConfigAndPathDerivation(unittest.TestCase):
    def test_runtime_paths_derive_from_common(self) -> None:
        cfg = {
            "common": {
                "period": "LHC99",
                "reco_pass": "apassX",
                "mc_production": "LHC99mc",
                "variant": "vreg",
                "base_input_dir": "/tmp/in/",
                "base_output_root": "/tmp/out/",
            },
        }
        runtime = s.current_runtime_config(cfg)

        self.assertEqual(runtime.paths.base_output_dir, "/tmp/out/LHC99/apassX/")
        self.assertEqual(runtime.paths.base_variant_output_dir, "/tmp/out/LHC99/apassX/vreg/")
        self.assertEqual(runtime.paths.data_tree_filename, "/tmp/in/data/LHC99/apassX/AO2D_skimmed.root")
        self.assertEqual(runtime.paths.mc_tree_filename, "/tmp/in/MC/LHC99mc/AO2D_coalescence.root")

        species_paths = runtime.paths.species_stage_paths("he3")
        self.assertEqual(species_paths["data_output"], "/tmp/out/LHC99/apassX/vreg/he3/DataHistos.root")
        self.assertEqual(species_paths["report_dir"], "/tmp/out/LHC99/apassX/vreg/he3/report")

    def test_merge_config_keeps_nested_defaults(self) -> None:
        merged = s.merge_config({"selections": {"he3": {"nsigma_tof": 2.4}}})

        self.assertIn("he3", merged["selections"])
        self.assertNotIn("he4", merged["selections"])
        self.assertIn("mc_reco_append", merged["selections"]["he3"])
        self.assertEqual(merged["selections"]["he3"]["nsigma_tof"], 2.4)

    def test_merge_config_selects_he4_defaults_when_requested(self) -> None:
        merged = s.merge_config({"run": {"species": "he4"}})

        self.assertEqual(merged["run"]["species"], "he4")
        self.assertIn("he4", merged["selections"])
        self.assertNotIn("he3", merged["selections"])
        self.assertIn("he4", merged["particle"])
        self.assertNotIn("he3", merged["particle"])

    def test_particle_overlay_applies_to_runtime_profile(self) -> None:
        runtime = s.current_runtime_config({"selections": {"he3": {"nsigma_tof": 2.1}}})
        he3 = runtime.get_particle_config("he3")
        self.assertAlmostEqual(float(he3["tof_nsigma_cut"]), 2.1)


if __name__ == "__main__":
    unittest.main()
