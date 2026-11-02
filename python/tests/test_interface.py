import itertools
from pathlib import Path
import tempfile
import unittest
import h5py
import numpy as np
import yaml

from sentinel_interface.config import ROOT, load
from sentinel_interface.field import prepare
from sentinel_interface.grid import cells, initial_batches
from sentinel_interface.runner import run

class InterfaceTests(unittest.TestCase):
    def config(self):
        return yaml.safe_load((ROOT / "configs/smoke_test.yaml").read_text())

    def test_volume_sampling_and_partial_cells(self):
        config = self.config()
        config["bounds"] = {"r": [0, 1.5], "z": [0, 1.5]}
        config["grid"]["electrons_per_cell"] = 10000
        table = cells(config)
        self.assertEqual(len(table), 4)
        batches = list(initial_batches(table, 10000, np.random.default_rng(42), batch_size=701))
        data = np.concatenate(batches)
        self.assertEqual(len(data), 40000)
        np.testing.assert_array_equal(data["event"], np.arange(40000))
        np.testing.assert_array_equal(np.bincount(data["cell_id"].astype(int)), [10000]*4)
        self.assertTrue((data["x0"] >= 0).all())
        self.assertTrue((data["y0"] == 0).all())
        for cell in table:
            rows = data[data["cell_id"] == cell["cell_id"]]
            self.assertTrue((rows["x0"] >= cell["r_min_mm"]).all())
            self.assertTrue((rows["x0"] <= cell["r_max_mm"]).all())
            self.assertTrue((rows["z0"] >= cell["z_min_mm"]).all())
            self.assertTrue((rows["z0"] <= cell["z_max_mm"]).all())
            expected = (cell["r_min_mm"]**2+cell["r_max_mm"]**2)/2
            self.assertAlmostEqual(np.mean(rows["x0"]**2), expected, delta=0.02)

    def test_units_and_validation(self):
        with tempfile.TemporaryDirectory() as work:
            path = Path(work) / "config.yaml"
            config = self.config()
            config["bounds"] = {"unit": "cm", "r": [0, .1], "z": [.2, .3]}
            path.write_text(yaml.safe_dump(config))
            parsed = load(path)
            self.assertEqual(parsed["bounds"]["r"], [0, 1])
            self.assertEqual(parsed["bounds"]["z"], [2, 3])
            for wrong in (False, 0, -1, float("nan")):
                config["grid"]["size_mm"] = wrong
                path.write_text(yaml.safe_dump(config))
                with self.assertRaises(ValueError):
                    load(path)

    def test_comsol_units_columns_and_missing_ez(self):
        with tempfile.TemporaryDirectory() as work:
            source, destination = Path(work)/"field.txt", Path(work)/"normalized.dat"
            header = "% Length unit: m\n% r  z  es.Er (V/m)  es.Ephi (V/m)  es.Ez (V/m)\n"
            source.write_text(header + "0.001 0.002 300 0 -400\n")
            mode, _ = prepare(source, destination)
            self.assertEqual(mode, "vector")
            np.testing.assert_allclose(np.loadtxt(destination), [1, 2, 3, -4])
            source.write_text(header.replace("es.Ez", "missing") + "0.001 0.002 300 0 -400\n")
            with self.assertRaisesRegex(ValueError, "Ez"):
                prepare(source, destination)

    def check_run(self, output, total, endpoint, trajectory, summary):
        with h5py.File(output/"initial_positions.h5") as f:
            positions = f["initial_positions"][:]
            mapping = f["event_cells"][:]
            self.assertEqual(len(positions), total)
            self.assertTrue((positions["y0"] == 0).all())
        for name, enabled in (("endpoints", endpoint), ("trajectories", trajectory), ("summary", summary)):
            self.assertEqual((output/f"{name}.h5").exists(), enabled)
            if enabled:
                with h5py.File(output/f"{name}.h5") as f:
                    np.testing.assert_array_equal(f["event_cells"][:], mapping)
                    if name == "summary":
                        self.assertEqual(f["summary"]["n_electrons_processed"].sum(), total)
        if endpoint:
            with h5py.File(output/"endpoints.h5") as f:
                endpoints = f["endpoints"][:]
                self.assertEqual(len(endpoints), total)
                np.testing.assert_array_equal(endpoints["event"], positions["event"])
                for key in ("x0", "y0", "z0"):
                    np.testing.assert_allclose(endpoints[key], positions[key], rtol=1e-6)
                self.assertTrue((endpoints["t"] > 0).all())
        if trajectory:
            with h5py.File(output/"trajectories.h5") as f:
                points = f["trajectories"][:]
            self.assertTrue(np.any(points["y"] != 0))  # diffusion is three dimensional
            self.assertEqual(len(np.unique(points["event"])), total)
            if endpoint:
                for row in endpoints:
                    track = points[points["event"] == row["event"]]
                    for end_key, point_key in (("x1","x"), ("y1","y"), ("z1","z"), ("t","t")):
                        self.assertEqual(row[end_key], track[-1][point_key])

    def test_all_output_switches_and_repeat_runs(self):
        with tempfile.TemporaryDirectory() as work:
            path = Path(work)/"config.yaml"
            for endpoint, trajectory, summary in itertools.product((False, True), repeat=3):
                config = self.config()
                config["grid"]["size_mm"] = 0.25
                config["output"] = dict(directory=str(Path(work)/"outputs"), save_endpoints=endpoint,
                                        save_trajectories=trajectory, save_summary=summary)
                path.write_text(yaml.safe_dump(config))
                output = run(path)
                self.check_run(output, 6, endpoint, trajectory, summary)
                self.assertTrue((output/"configs/initial_conditions.yaml").is_file())
            self.assertEqual(len(list((Path(work)/"outputs").iterdir())), 8)

    def test_real_comsol(self):
        source = ROOT/"data/fields/Full_output.txt"
        if not source.is_file():
            self.skipTest("Real COMSOL field unavailable")
        with tempfile.TemporaryDirectory() as work:
            config = self.config()
            config["grid"]["electrons_per_cell"] = 4
            config["bounds"] = {"unit":"mm", "r":[.5,.9], "z":[50,50.5]}
            config["electric_field"]["file"] = str(source)
            config["output"]["directory"] = str(Path(work)/"outputs")
            path = Path(work)/"config.yaml"
            path.write_text(yaml.safe_dump(config))
            self.check_run(run(path), 4, True, True, True)

if __name__ == "__main__":
    unittest.main()
