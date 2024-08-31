# coding: utf-8

import os
import sansmic
import unittest
import numpy as np


class TestModelRun(unittest.TestCase):
    def runTest(self):
        filename = "old.dat"
        model = sansmic.read_scenario(filename)
        sansmic.write_scenario(model, "new.toml")
        with model.new_simulation("new") as sim:
            sim.run_sim()
        resPy = sim.results
        resF = sansmic.io.read_classic_out_ddl("old")
        diff = resF.summary.iloc[0:, :] - resPy.summary.iloc[1:, :].reset_index()
        return resF, resPy


if __name__ == "__main__":
    test = TestModelRun()
    resF, resPy = test.runTest()
