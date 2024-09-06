# coding: utf-8

import os
import sansmic
import unittest
import numpy as np


class TestModelRun(unittest.TestCase):
    def runTest(self):
        filename = "baseline.dat"
        model = sansmic.read_scenario(filename)
        sansmic.write_scenario(model, "regression.toml")
        with model.new_simulation("regression") as sim:
            sim.run_sim()
        resPy = sim.results
        resF = sansmic.io.read_classic_out_ddl("baseline")
        diff = resF.by_time.iloc[0:, :] - resPy.by_time.iloc[1:, :].reset_index()
        return resF, resPy, diff


if __name__ == "__main__":
    test = TestModelRun()
    resF, resPy, diff = test.runTest()
