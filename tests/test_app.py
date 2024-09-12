# coding: utf-8

import glob
import os
import tempfile
import unittest
from os.path import abspath, dirname, join

import numpy as np
import pandas as pd
import sansmic
import sansmic.app
import sansmic.io
import sansmic.model

testdir = dirname(abspath(str(__file__)))


class TestApplication(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Set up temporary directory
        cls.tempdir = tempfile.TemporaryDirectory()
        cls.tempdirname = cls.tempdir.name
        # cls.tempdirname = abspath('.')

        # Set temporary filenames
        cls.withdrawal_dat = os.path.join(cls.tempdirname, "withdrawal.dat")
        cls.withdrawal_toml = os.path.join(cls.tempdirname, "withdrawal.toml")

        # Write the .DAT formatted data
        with open(cls.withdrawal_dat, "w") as fin:
            fin.write(cls._WITHDRAWAL_DAT)

    def test_run_app(self):
        """Test the 'sansmic' command line application."""
        scenario = sansmic.io.read_scenario(self.withdrawal_dat)
        with scenario.new_simulation() as sim:
            sim.run_sim()
        res0 = sim.results
        res1 = sansmic.app.main(
            args=[
                self.withdrawal_dat,
                "--no-toml",
                "--no-csv",
                "--no-json",
                "--no-hdf",
                "--no-tst",
                "--no-old-out",
            ],
            ret=True,
        )
        res2 = sansmic.app.main(
            args=[
                self.withdrawal_dat,
                "-o",
                join(self.tempdirname, "app-test"),
                "--toml",
                "--csv",
                "--json",
                "--hdf",
                "--tst",
                "--old-out",
            ],
            ret=False,
        )
        res3 = sansmic.app.main(
            args=[
                self.withdrawal_dat,
                "--no-toml",
                "--no-csv",
                "--no-json",
                "--no-hdf",
                "--no-tst",
                "--no-old-out",
                "-v",
            ],
            ret=True,
        )
        res4 = sansmic.app.main(
            args=[
                self.withdrawal_dat,
                "--no-toml",
                "--no-csv",
                "--no-json",
                "--no-hdf",
                "--no-tst",
                "--no-old-out",
                "-v",
                "-v",
            ],
            ret=True,
        )
        res5 = sansmic.app.main(
            args=[
                self.withdrawal_dat,
                "--no-toml",
                "--no-csv",
                "--no-json",
                "--no-hdf",
                "--no-tst",
                "--no-old-out",
                "-q",
            ],
            ret=True,
        )
        self.assertIsNotNone(res1)
        self.assertIsNone(res2)

        self.assertEqual(len(glob.glob(join(self.tempdirname, "app-test.toml"))), 1)
        self.assertEqual(len(glob.glob(join(self.tempdirname, "app-test.log"))), 1)
        self.assertEqual(len(glob.glob(join(self.tempdirname, "app-test*.csv"))), 7)
        self.assertEqual(len(glob.glob(join(self.tempdirname, "app-test.json"))), 1)
        self.assertEqual(len(glob.glob(join(self.tempdirname, "app-test.h5"))), 1)
        self.assertEqual(len(glob.glob(join(self.tempdirname, "app-test.tst"))), 1)
        self.assertEqual(len(glob.glob(join(self.tempdirname, "app-test.out"))), 1)
        self.assertEqual(len(glob.glob(join(self.tempdirname, "app-test*"))), 13)

        self.assertTrue((res0.df_t_1D == res1.df_t_1D).all().all())
        self.assertTrue((res0.df_z_1D == res1.df_z_1D).all().all())
        self.assertTrue((res0.df_t_z_2D == res1.df_t_z_2D).all().all())

        self.assertTrue((res0.df_t_1D == res3.df_t_1D).all().all())
        self.assertTrue((res0.df_z_1D == res3.df_z_1D).all().all())
        self.assertTrue((res0.df_t_z_2D == res3.df_t_z_2D).all().all())

        self.assertTrue((res0.df_t_1D == res4.df_t_1D).all().all())
        self.assertTrue((res0.df_z_1D == res4.df_z_1D).all().all())
        self.assertTrue((res0.df_t_z_2D == res4.df_t_z_2D).all().all())

        self.assertTrue((res0.df_t_1D == res5.df_t_1D).all().all())
        self.assertTrue((res0.df_z_1D == res5.df_z_1D).all().all())
        self.assertTrue((res0.df_t_z_2D == res5.df_t_z_2D).all().all())

    def test_convert_app(self):
        """Test the 'sansmic-convert' command line application."""
        scenario0 = sansmic.io.read_scenario(self.withdrawal_dat)
        scenario0.title = ""
        sansmic.app.convert(args=[self.withdrawal_dat, self.withdrawal_toml])
        scenario1 = sansmic.io.read_scenario(self.withdrawal_toml)
        scenario1.title = ""
        self.assertEqual(scenario0, scenario1)
        sansmic.app.convert(args=[self.withdrawal_dat, self.withdrawal_toml, "--full"])
        scenario2 = sansmic.io.read_scenario(self.withdrawal_toml)
        scenario2.title = ""
        self.maxDiff = None
        print(scenario0.stages[0].product_injection_rate)
        print(scenario1.stages[0].product_injection_rate)
        print(scenario2.stages[0].product_injection_rate)
        self.assertDictEqual(scenario0.to_dict(), scenario2.to_dict())

    @classmethod
    def tearDownClass(cls):
        cls.tempdir.cleanup()

    _WITHDRAWAL_DAT = """Withdrawal leach
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

    _WITHDRAWAL_TOML = """num-cells = 100
geometry-data.radii = [10.0, 50.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 90.0, 80.0, 50.0, 33.0, 10.0, 4.0, 2.0, 2.0, 2.0, 2.0, 2.0]
geometry-format = 'radius-list'
cavern-height = 1000.0
floor-depth = 4000.0
ullage-standoff = 20.0
insolubles-ratio = 0.05
units = 'ft-in-bbl'

[[stages]]
title = 'Withdrawal leach'
simulation-mode = 'withdrawal'
solver-timestep = 0.1
save-frequency = 240
injection-duration = 72.0
rest-duration = 1440.0
inner-tbg-inside-diam = 9.85
inner-tbg-outside-diam = 10.75
outer-csg-inside-diam = 9.85
outer-csg-outside-diam = 10.75
brine-injection-sg = 1.03
brine-injection-depth = 32.0
brine-production-depth = 1000.0
brine-injection-rate = 240000.0
set-cavern-sg = 1.2019
brine-interface-depth = 100.0
"""
