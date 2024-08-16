# coding: utf-8

import unittest
import sansmic
import sansmic.outputs

class TestModelRun(unittest.TestCase):
    def runTest(self):
        filename = 'ordinary-f90.dat'
        model = sansmic.read_scenario(filename)
        sansmic.write_scenario(model, 'ordinary-c++.toml')
        with model.new_simulation('ordinary-c++') as sim:
            sim.run()
            resPy = sim.get_results()
        resF = sansmic.outputs.read_classic('ordinary-f90')
        res_1 = resPy
        res_0 = resF
        resPy = resPy.get_timeseries_data()
        resF = resF.get_timeseries_data()
        resPy['modtime'] = (((resPy['time']*100)//10) / 10.0)
        resF['modtime'] = (((resF['time']*100)//10) / 10.0)
        resPy.set_index('modtime',inplace=True)
        resF.set_index('modtime',inplace=True)
        bothtimes = list(set(resF.index).intersection(set(resPy.index)))
        bothtimes.sort()
        self.assertLess((resPy.loc[bothtimes,'cavern volume']-resF.loc[bothtimes,'cavern volume']).abs().max(),1e-4*resF['cavern volume'].max(), 'Volumes not within 0.01% through the entrie simulation')


if __name__ == '__main__':
    test = TestModelRun()
    test.runTest()
