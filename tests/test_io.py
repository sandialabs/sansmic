# coding: utf-8

import os
import unittest
from os.path import abspath, dirname, join

import numpy as np
import pandas as pd
import sansmic, sansmic.io, sansmic.model
import tempfile

testdir = dirname(abspath(str(__file__)))


class TestReadAndWriteScenarios(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Set up temporary directory
        cls.tempdir = tempfile.TemporaryDirectory()
        cls.tempdirname = cls.tempdir.name
        # cls.tempdirname = abspath('tests')

        # Set temporary filenames for the static data
        cls.withdrawal_dat = os.path.join(cls.tempdirname, "withdrawal.dat")
        cls.ordinary_dir_dat = os.path.join(cls.tempdirname, "ordinary-dir.dat")
        cls.ordinary_rev_dat = os.path.join(cls.tempdirname, "ordinary-rev.dat")
        cls.leach_fill_dat = os.path.join(cls.tempdirname, "leach-fill.txt")

        # Write the .DAT formatted data
        with open(cls.withdrawal_dat, "w") as fin:
            fin.write(cls._WITHDRAWAL_DAT_DATA)
        with open(cls.ordinary_dir_dat, "w") as fin:
            fin.write(cls._ORDINARY_DIRECT_DATA)
        with open(cls.ordinary_rev_dat, "w") as fin:
            fin.write(cls._ORDINARY_REVERSE_DATA)
        with open(cls.leach_fill_dat, "w") as fin:
            fin.write(cls._LEACH_FILL_DATA)

        # Set temporary filenames
        cls.withdrawal_toml = os.path.join(cls.tempdirname, "withdrawal.toml")
        cls.ordinary_dir_json = os.path.join(cls.tempdirname, "ordinary-dir.json")
        cls.ordinary_rev_yaml = os.path.join(cls.tempdirname, "ordinary-rev.yaml")
        cls.leach_fill_format = os.path.join(cls.tempdirname, "leach-fill-new.txt")

    def test_read_scenario_leachfill(self):
        scenario = sansmic.io.read_scenario(self.leach_fill_dat, format="dat")
        self.assertEqual(len(scenario.stages), 1)
        self.assertEqual(scenario.stages[0].title, "Leach fill")
        self.assertEqual(scenario.num_cells, 100)
        self.assertEqual(scenario.stages[0].simulation_mode, 2)
        self.assertIs(
            scenario.stages[0].simulation_mode, sansmic.model.SimulationMode.LEACH_FILL
        )
        self.assertEqual(scenario.stages[0].save_frequency, 240)
        self.assertTrue(scenario.stages[0].set_initial_conditions)
        self.assertEqual(scenario.stages[0].rest_duration, 1440)
        self.assertIn(scenario.advanced.coallescing_wells, [1, None])
        self.assertEqual(scenario.geometry_format, 0)
        self.assertIs(
            scenario.geometry_format, sansmic.model.GeometryFormat.RADIUS_LIST
        )
        self.assertEqual(scenario.stages[0].stop_value, -200.0)
        self.assertIs(
            scenario.stages[0].stop_condition, sansmic.model.StopCondition.DEPTH
        )
        self.assertEqual(scenario.cavern_height, 1000.0)
        self.assertEqual(scenario.stages[0].brine_injection_depth, 32.0)
        self.assertEqual(scenario.stages[0].brine_production_depth, 350.0)
        self.assertEqual(scenario.stages[0].brine_interface_depth, 900.0)
        self.assertEqual(scenario.ullage_standoff, 20.0)
        self.assertEqual(scenario.stages[0].brine_injection_rate, 240000.0)
        self.assertEqual(scenario.stages[0].inner_tbg_inside_diam / 2, 4.925)
        self.assertEqual(scenario.stages[0].inner_tbg_outside_diam / 2, 5.375)
        self.assertEqual(scenario.stages[0].outer_csg_inside_diam / 2, 4.925)
        self.assertEqual(scenario.stages[0].outer_csg_outside_diam / 2, 5.375)
        self.assertEqual(scenario.stages[0].brine_injection_sg, 1.03)
        self.assertEqual(scenario.stages[0].set_cavern_sg, 1.2019)
        self.assertEqual(scenario.stages[0].solver_timestep, 0.1)
        self.assertEqual(scenario.stages[0].injection_duration, 72)
        self.assertEqual(scenario.stages[0].product_injection_rate, 240000.0)
        self.assertIsInstance(scenario.geometry_data, dict)
        self.assertIn("radii", scenario.geometry_data)
        self.assertEqual(len(scenario.geometry_data["radii"]), 101)
        self.assertEqual(scenario.floor_depth, 4000.0)
        self.assertEqual(scenario.advanced.dissolution_factor, 1.01)
        self.assertEqual(scenario.insolubles_ratio, 0.05)

    def test_read_scenario_withdrawal(self):
        scenario = sansmic.io.read_scenario(self.withdrawal_dat)
        self.assertEqual(len(scenario.stages), 1)
        self.assertEqual(scenario.stages[0].title, "Withdrawal leach")
        self.assertEqual(scenario.num_cells, 100)
        self.assertEqual(scenario.stages[0].simulation_mode, 1)
        self.assertIs(
            scenario.stages[0].simulation_mode, sansmic.model.SimulationMode.WITHDRAWAL
        )
        self.assertEqual(scenario.stages[0].save_frequency, 240)
        self.assertTrue(scenario.stages[0].set_initial_conditions)
        self.assertEqual(scenario.stages[0].rest_duration, 1440)
        self.assertIn(scenario.advanced.coallescing_wells, [1, None])
        self.assertEqual(scenario.geometry_format, 0)
        self.assertIs(
            scenario.geometry_format, sansmic.model.GeometryFormat.RADIUS_LIST
        )
        self.assertEqual(scenario.stages[0].stop_value, 0.0)
        self.assertIs(
            scenario.stages[0].stop_condition, sansmic.model.StopCondition.DURATION
        )
        self.assertEqual(scenario.cavern_height, 1000.0)
        self.assertEqual(scenario.stages[0].brine_injection_depth, 32.0)
        self.assertEqual(scenario.stages[0].brine_production_depth, 1000.0)
        self.assertEqual(scenario.stages[0].brine_interface_depth, 100.0)
        self.assertEqual(scenario.ullage_standoff, 20.0)
        self.assertEqual(scenario.stages[0].brine_injection_rate, 240000.0)
        self.assertEqual(scenario.stages[0].inner_tbg_inside_diam / 2, 4.925)
        self.assertEqual(scenario.stages[0].inner_tbg_outside_diam / 2, 5.375)
        self.assertEqual(scenario.stages[0].outer_csg_inside_diam / 2, 4.925)
        self.assertEqual(scenario.stages[0].outer_csg_outside_diam / 2, 5.375)
        self.assertEqual(scenario.stages[0].brine_injection_sg, 1.03)
        self.assertEqual(scenario.stages[0].set_cavern_sg, 1.2019)
        self.assertEqual(scenario.stages[0].solver_timestep, 0.1)
        self.assertEqual(scenario.stages[0].injection_duration, 72)
        self.assertEqual(scenario.stages[0].product_injection_rate, 0.0)
        self.assertIsInstance(scenario.geometry_data, dict)
        self.assertIn("radii", scenario.geometry_data)
        self.assertEqual(len(scenario.geometry_data["radii"]), 101)
        self.assertEqual(scenario.floor_depth, 4000.0)
        self.assertIn(scenario.advanced.dissolution_factor, [1.0, None])
        self.assertEqual(scenario.insolubles_ratio, 0.05)

    def test_read_scenario_ordinary_direct(self):
        scenario = sansmic.io.read_scenario(self.ordinary_dir_dat)
        self.assertEqual(len(scenario.stages), 1)
        self.assertEqual(scenario.stages[0].title, "Ordinary leach - direct")
        self.assertEqual(scenario.num_cells, 100)
        self.assertEqual(scenario.stages[0].simulation_mode, 0)
        self.assertIs(
            scenario.stages[0].simulation_mode, sansmic.model.SimulationMode.ORDINARY
        )
        self.assertEqual(scenario.stages[0].save_frequency, 240)
        self.assertTrue(scenario.stages[0].set_initial_conditions)
        self.assertEqual(scenario.stages[0].rest_duration, 1440)
        self.assertIn(scenario.advanced.coallescing_wells, [1, None])
        self.assertEqual(scenario.geometry_format, 0)
        self.assertIs(
            scenario.geometry_format, sansmic.model.GeometryFormat.RADIUS_LIST
        )
        self.assertEqual(scenario.stages[0].stop_value, 4000000.0)
        self.assertIs(
            scenario.stages[0].stop_condition, sansmic.model.StopCondition.VOLUME
        )
        self.assertEqual(scenario.cavern_height, 1000.0)
        self.assertEqual(scenario.stages[0].brine_injection_depth, 100.0)
        self.assertEqual(scenario.stages[0].brine_production_depth, 800.0)
        self.assertEqual(scenario.stages[0].brine_interface_depth, 900.0)
        self.assertEqual(scenario.ullage_standoff, 20.0)
        self.assertEqual(scenario.stages[0].brine_injection_rate, 240000.0)
        self.assertEqual(scenario.stages[0].inner_tbg_inside_diam / 2, 4.925)
        self.assertEqual(scenario.stages[0].inner_tbg_outside_diam / 2, 5.375)
        self.assertEqual(scenario.stages[0].outer_csg_inside_diam / 2, 6.925)
        self.assertEqual(scenario.stages[0].outer_csg_outside_diam / 2, 7.375)
        self.assertEqual(scenario.stages[0].brine_injection_sg, 1.1)
        self.assertEqual(scenario.stages[0].set_cavern_sg, 1.2019)
        self.assertEqual(scenario.stages[0].solver_timestep, 0.1)
        self.assertEqual(scenario.stages[0].injection_duration, 72)
        self.assertEqual(scenario.stages[0].product_injection_rate, 0.0)
        self.assertIsInstance(scenario.geometry_data, dict)
        self.assertIn("radii", scenario.geometry_data)
        self.assertEqual(len(scenario.geometry_data["radii"]), 101)
        self.assertEqual(scenario.floor_depth, 4000.0)
        self.assertIn(scenario.advanced.dissolution_factor, [1.0, None])
        self.assertEqual(scenario.insolubles_ratio, 0.05)

    def test_read_scenario_ordinary_reverse_multistage(self):
        scenario = sansmic.io.read_scenario(self.ordinary_rev_dat)
        self.assertEqual(len(scenario.stages), 3)
        self.assertEqual(scenario.stages[0].title, "Ordinary leach stage 1")
        self.assertEqual(scenario.num_cells, 100)
        self.assertEqual(scenario.stages[0].simulation_mode, 0)
        self.assertIs(
            scenario.stages[0].simulation_mode, sansmic.model.SimulationMode.ORDINARY
        )
        self.assertEqual(scenario.stages[0].save_frequency, 240)
        self.assertTrue(scenario.stages[0].set_initial_conditions)
        self.assertEqual(scenario.stages[0].rest_duration, 1440)
        self.assertIn(scenario.advanced.coallescing_wells, [1, None])
        self.assertEqual(scenario.geometry_format, 0)
        self.assertIs(
            scenario.geometry_format, sansmic.model.GeometryFormat.RADIUS_LIST
        )
        self.assertEqual(scenario.stages[0].stop_value, 0.0)
        self.assertIs(
            scenario.stages[0].stop_condition, sansmic.model.StopCondition.DURATION
        )
        self.assertEqual(scenario.cavern_height, 1000.0)
        self.assertEqual(scenario.stages[0].brine_injection_depth, 800.0)
        self.assertEqual(scenario.stages[0].brine_production_depth, 100.0)
        self.assertEqual(scenario.stages[0].brine_interface_depth, 900.0)
        self.assertEqual(scenario.ullage_standoff, 20.0)
        self.assertEqual(scenario.stages[0].brine_injection_rate, 240000.0)
        self.assertEqual(scenario.stages[0].inner_tbg_inside_diam / 2, 4.925)
        self.assertEqual(scenario.stages[0].inner_tbg_outside_diam / 2, 5.375)
        self.assertEqual(scenario.stages[0].outer_csg_inside_diam / 2, 4.925)
        self.assertEqual(scenario.stages[0].outer_csg_outside_diam / 2, 5.375)
        self.assertEqual(scenario.stages[0].brine_injection_sg, 1.03)
        self.assertEqual(scenario.stages[0].set_cavern_sg, 1.2019)
        self.assertEqual(scenario.stages[0].solver_timestep, 0.1)
        self.assertEqual(scenario.stages[0].injection_duration, 72)
        self.assertEqual(scenario.stages[0].product_injection_rate, 0.0)
        self.assertIsInstance(scenario.geometry_data, dict)
        self.assertIn("radii", scenario.geometry_data)
        self.assertEqual(len(scenario.geometry_data["radii"]), 101)
        self.assertEqual(scenario.floor_depth, 4000.0)
        self.assertIn(scenario.advanced.dissolution_factor, [1.0, None])
        self.assertEqual(scenario.insolubles_ratio, 0.05)
        self.assertEqual(scenario.stages[1].title, "Ordinary leach stage 2")
        self.assertEqual(scenario.stages[1].simulation_mode, 0)
        self.assertIs(
            scenario.stages[1].simulation_mode, sansmic.model.SimulationMode.ORDINARY
        )
        self.assertEqual(scenario.stages[1].save_frequency, 240)
        self.assertFalse(scenario.stages[1].set_initial_conditions)
        self.assertEqual(scenario.stages[1].rest_duration, 1440)
        self.assertEqual(scenario.stages[1].stop_value, 0.0)
        self.assertIs(
            scenario.stages[1].stop_condition, sansmic.model.StopCondition.DURATION
        )
        self.assertEqual(scenario.stages[1].brine_injection_depth, 800.0)
        self.assertEqual(scenario.stages[1].brine_production_depth, 100.0)
        self.assertEqual(scenario.stages[1].brine_interface_depth, 0.0)
        self.assertEqual(scenario.stages[1].brine_injection_rate, 240000.0)
        self.assertEqual(scenario.stages[1].inner_tbg_inside_diam / 2, 4.925)
        self.assertEqual(scenario.stages[1].inner_tbg_outside_diam / 2, 5.375)
        self.assertEqual(scenario.stages[1].outer_csg_inside_diam / 2, 6.925)
        self.assertEqual(scenario.stages[1].outer_csg_outside_diam / 2, 7.375)
        self.assertEqual(scenario.stages[1].brine_injection_sg, 1.03)
        self.assertIsNone(scenario.stages[1].set_cavern_sg)
        self.assertEqual(scenario.stages[1].solver_timestep, 0.1)
        self.assertEqual(scenario.stages[1].injection_duration, 72)
        self.assertEqual(scenario.stages[1].product_injection_rate, 0.0)
        self.assertEqual(scenario.stages[2].title, "Ordinary leach stage 3")
        self.assertEqual(scenario.stages[2].simulation_mode, 0)
        self.assertIs(
            scenario.stages[2].simulation_mode, sansmic.model.SimulationMode.ORDINARY
        )
        self.assertEqual(scenario.stages[2].save_frequency, 240)
        self.assertFalse(scenario.stages[2].set_initial_conditions)
        self.assertEqual(scenario.stages[2].rest_duration, 1440)
        self.assertEqual(scenario.stages[2].stop_value, 0.0)
        self.assertIs(
            scenario.stages[2].stop_condition, sansmic.model.StopCondition.DURATION
        )
        self.assertEqual(scenario.stages[2].brine_injection_depth, 800.0)
        self.assertEqual(scenario.stages[2].brine_production_depth, 100.0)
        self.assertEqual(scenario.stages[2].brine_interface_depth, 0.0)
        self.assertEqual(scenario.stages[2].brine_injection_rate, 240000.0)
        self.assertEqual(scenario.stages[2].inner_tbg_inside_diam / 2, 4.925)
        self.assertEqual(scenario.stages[2].inner_tbg_outside_diam / 2, 5.375)
        self.assertEqual(scenario.stages[2].outer_csg_inside_diam / 2, 8.925)
        self.assertEqual(scenario.stages[2].outer_csg_outside_diam / 2, 9.375)
        self.assertEqual(scenario.stages[2].brine_injection_sg, 1.03)
        self.assertEqual(scenario.stages[2].set_cavern_sg, 1.2019)
        self.assertEqual(scenario.stages[2].solver_timestep, 0.1)
        self.assertEqual(scenario.stages[2].injection_duration, 72)
        self.assertEqual(scenario.stages[2].product_injection_rate, 0.0)

    def test_roundtrip_toml(self):
        scenario = sansmic.io.read_scenario(self.withdrawal_dat)
        scenario.title = ""
        sansmic.io.write_scenario(scenario, self.withdrawal_toml)
        scenario2 = sansmic.io.read_scenario(self.withdrawal_toml)
        scenario2.title = ""
        self.assertEqual(scenario, scenario2)

    def test_roundtrip_yaml(self):
        self.maxDiff = None
        scenario = sansmic.io.read_scenario(self.ordinary_rev_dat)
        scenario.title = ""
        scenario.stages[1].brine_interface_depth = 0
        scenario.stages[2].brine_interface_depth = 0
        scenario.stages[1].set_initial_conditions = None
        scenario.stages[2].set_initial_conditions = None

        sansmic.io.write_scenario(scenario, self.ordinary_rev_yaml, redundant=True)
        scenario2 = sansmic.io.read_scenario(self.ordinary_rev_yaml)
        scenario2.title = ""
        self.assertDictEqual(
            scenario.to_dict(keep_empty=True), scenario2.to_dict(keep_empty=True)
        )

    def test_roundtrip_json(self):
        scenario = sansmic.io.read_scenario(self.ordinary_dir_dat)
        scenario.title = ""
        sansmic.io.write_scenario(scenario, self.ordinary_dir_json)
        scenario2 = sansmic.io.read_scenario(self.ordinary_dir_json)
        scenario2.title = ""
        self.assertEqual(scenario, scenario2)

    def test_roundtrip_spec_format(self):
        scenario = sansmic.io.read_scenario(self.leach_fill_dat, format="dat")
        scenario.title = ""

        sansmic.io.write_scenario(scenario, self.leach_fill_format, format="toml")
        scenario2 = sansmic.io.read_scenario(self.leach_fill_format, format="toml")
        scenario2.title = ""
        self.assertEqual(scenario, scenario2)

        sansmic.io.write_scenario(scenario, self.leach_fill_format, format="json")
        scenario3 = sansmic.io.read_scenario(self.leach_fill_format, format="json")
        scenario3.title = ""
        self.assertEqual(scenario, scenario3)

        sansmic.io.write_scenario(scenario, self.leach_fill_format, format="yaml")
        scenario4 = sansmic.io.read_scenario(self.leach_fill_format, format="yaml")
        scenario4.title = ""
        self.assertEqual(scenario, scenario4)

    @classmethod
    def tearDownClass(cls):
        cls.tempdir.cleanup()

    _WITHDRAWAL_DAT_DATA = """Withdrawal leach
100 1   240 0   0   1440    1   0   0
1000.   32.0    1000.0  100.0   20.0
240000.0
4.925   5.3750  4.925   5.3750
1.03    1.2019
0.1     72
0.0 0.0 0.0
10
50
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
90
80
50
33
10
4
2
2
2
2
2
1.0   0.05   4000.   4000.
END
"""

    _ORDINARY_DIRECT_DATA = """Ordinary leach - direct
100	0	240	0	0	1440	 1	0	4000000.0
1000.   100.0  800.0	900.0   20.0
240000.0
4.925	5.3750	6.925	7.3750
1.1 1.2019
0.1     72
0.0	0.0	0.0
10
50
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
90
80
50
33
10
4
2
2
2
2
2
1.0   0.05   4000.   4000.
END
"""

    _ORDINARY_REVERSE_DATA = """Ordinary leach stage 1
100	0	240	0	0	1440	 1	0	0
1000.   800.0  100.0	900.0   20.0
240000.0
4.925	5.3750	4.925	5.3750
1.03	1.2019
0.1     72
0.0	0.0	0.0
10
50
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
90
80
50
33
10
4
2
2
2
2
2
1.0   0.05   4000.   4000.
Ordinary leach stage 2
100	0	240	1	0	1440	 1	0	0
1000.   800.0  100.0	0.0   20.0
240000.0
4.925	5.3750	6.925	7.3750
1.03	0.0
0.1     72
0.0	0.0	0.0
Ordinary leach stage 3
100	0	240	1	0	1440	 1	0	0
1000.   800.0  100.0	0.0   20.0
240000.0
4.925	5.3750	8.925	9.3750
1.03	1.2019
0.1     72
0.0	0.0	0.0
END
"""

    _LEACH_FILL_DATA = """Leach fill
100	2	240	0	0	1440	 1	0	-200.0
1000.   32.0	350.0  900.0   20.0
240000.0
4.925	5.3750	4.925	5.3750
1.03	1.2019
0.1     72
240000.0	0.0	0.0
10
50
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
100
90
80
50
33
10
4
2
2
2
2
2
1.01   0.05   4000.   4000.
END
"""


if __name__ == "__main__":
    unittest.main()
