# coding: utf-8

import unittest
import sansmic

class TestModelRun(unittest.TestCase):
    def runTest(self):
        filename = 'old.dat'
        model = sansmic.read_scenario(filename)
        sansmic.write_scenario(model, 'new.toml')
        with model.new_simulation('new') as sim:
            sim.run()
            resPy = sim.get_results()
        resF = sansmic.io.read_classic_out_ddl('old')
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
        # self.assertLess((resPy.loc[bothtimes,'cavern volume']-resF.loc[bothtimes,'cavern volume']).abs().max(),1e-4*resF['cavern volume'].max(), 'Volumes not within 0.01% through the entrie simulation')
        return resF, resPy

resF = resPy = None

if __name__ == '__main__':
    test = TestModelRun()
    resF, resPy = test.runTest()
