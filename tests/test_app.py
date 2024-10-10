# coding: utf-8

import glob
import os
import subprocess
import tempfile
import unittest
from os.path import abspath, dirname, join

import click
import click.testing
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
        res1 = sansmic.app.run(
            [
                self.withdrawal_dat,
                "--no-csv",
                "--no-toml",
                "--no-json",
                "--no-hdf",
                "--no-tst",
                "--no-oldout",
            ],
            standalone_mode=False,
        )
        res2 = sansmic.app.run(
            [
                self.withdrawal_dat,
                "--prefix",
                join(self.tempdirname, "app-test"),
                "--csv",
                "--toml",
                "--json",
                "--hdf",
                "--tst",
                "--oldout",
            ],
            standalone_mode=False,
        )
        res3 = sansmic.app.run(
            [
                self.withdrawal_dat,
                "--no-csv",
                "--no-hdf",
                "--no-tst",
                "-v",
            ],
            standalone_mode=False,
        )
        res4 = sansmic.app.run(
            [
                self.withdrawal_dat,
                "--no-toml",
                "--no-csv",
                "--no-json",
                "--no-hdf",
                "--no-tst",
                "--no-oldout",
                "-vv",
            ],
            standalone_mode=False,
        )
        res5 = sansmic.app.run(
            [
                self.withdrawal_dat,
                "--no-toml",
                "--no-csv",
                "--no-json",
                "--no-hdf",
                "--no-tst",
                "--no-oldout",
                "-q",
            ],
            standalone_mode=False,
        )

        res6 = sansmic.app.run(
            [
                self.withdrawal_dat,
                "--debug",
            ],
            standalone_mode=False,
        )

        self.assertIsNotNone(res1)

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

        self.assertTrue((res0.df_t_1D == res6.df_t_1D).all().all())
        self.assertTrue((res0.df_z_1D == res6.df_z_1D).all().all())
        self.assertTrue((res0.df_t_z_2D == res6.df_t_z_2D).all().all())

        with self.assertRaises(
            (click.FileError, click.ClickException, FileNotFoundError)
        ):
            sansmic.app.run(
                [
                    "thisWontWork.TOML",
                ],
                standalone_mode=False,
            )

        with self.assertRaises(
            (click.FileError, click.ClickException, FileNotFoundError)
        ):
            sansmic.app.run(
                [
                    self.withdrawal_dat,
                    "-o",
                    os.path.abspath(os.path.join("this", "is", "a", "fake.", "path")),
                ],
                standalone_mode=False,
            )

    def test_license(self):
        """Test the license file echo"""
        runner = click.testing.CliRunner()
        results = runner.invoke(sansmic.app.run, ["--license"])
        self.assertEqual(results.output.strip(), sansmic.__license__.strip())

    def test_copyright(self):
        """Test the license file echo"""
        runner = click.testing.CliRunner()
        results = runner.invoke(sansmic.app.run, ["--copyright"])
        self.assertEqual(results.output.strip(), sansmic.__copyright__.strip())

    def test_version(self):
        """Test the license file echo"""
        runner = click.testing.CliRunner()
        results = runner.invoke(sansmic.app.run, ["--version"])
        self.assertEqual(results.output.strip(), sansmic.__version__.strip())

    def test_convert_app(self):
        """Test the 'sansmic-convert' command line application."""
        scenario0 = sansmic.io.read_scenario(self.withdrawal_dat)
        scenario0.title = ""
        sansmic.app.convert(
            [self.withdrawal_dat, self.withdrawal_toml], standalone_mode=False
        )
        scenario1 = sansmic.io.read_scenario(self.withdrawal_toml)
        scenario1.title = ""
        self.assertEqual(scenario0, scenario1)
        sansmic.app.convert(
            [self.withdrawal_dat, self.withdrawal_toml, "--full"], standalone_mode=False
        )
        scenario2 = sansmic.io.read_scenario(self.withdrawal_toml)
        scenario2.title = ""
        self.maxDiff = None
        self.assertDictEqual(scenario0.to_dict(), scenario2.to_dict())

        with self.assertRaises(
            (click.FileError, click.ClickException, FileNotFoundError)
        ):
            sansmic.app.convert(
                ["thisWontWork.TOML", "neither.will.this."], standalone_mode=False
            )

    def test_sansmic_cmd(self):
        proc = subprocess.run(
            [
                "sansmic",
                self.withdrawal_dat,
                "-o",
                join(self.tempdirname, "sansmic-call-test"),
                "-v",
            ],
            capture_output=True,
            # avoid having to explicitly encode
            text=True,
        )
        data = proc.stdout
        result = proc.returncode
        self.assertTrue("Progress..." in data)
        self.assertTrue("960000.0" in data)
        print(glob.glob(join(self.tempdirname, "sansmic-call-test*")))
        self.assertEqual(
            len(glob.glob(join(self.tempdirname, "sansmic-call-test*"))), 10
        )

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
Stage 2
100 1   240 0   0   120    1   0   0
1000.   32.0    1000.0  100.0   20.0
240000.0
4.925   5.3750  4.925   5.3750
1.03    1.2019
0.1     24
0.0 0.0 0.0
END
"""
