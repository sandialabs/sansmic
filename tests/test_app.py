# coding: utf-8

import os
import unittest
from os.path import abspath, dirname, join

import numpy as np
import pandas as pd
import sansmic, sansmic.io, sansmic.model, sansmic.app
import tempfile

testdir = dirname(abspath(str(__file__)))


class TestApplication(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Set up temporary directory
        cls.tempdir = tempfile.TemporaryDirectory()
        cls.tempdirname = cls.tempdir.name
        # cls.tempdirname = abspath('tests')

        # Set temporary filenames
        cls.withdrawal_dat = os.path.join(cls.tempdirname, "withdrawal.dat")
        cls.withdrawal_toml = os.path.join(cls.tempdirname, "withdrawal.toml")

        # Write the .DAT formatted data
        with open(cls.withdrawal_dat, "w") as fin:
            fin.write(cls._WITHDRAWAL_DAT)

    def test_run_app(self):
        ret = sansmic.app.main(args=[self.withdrawal_dat, "--toml", "--csv"], ret=True)

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
