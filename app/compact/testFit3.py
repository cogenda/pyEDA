from pyEDA.Compact.MOSParamFit import *
from pyEDA.Compact.AuroraData import *
import math
import numpy as np
from matplotlib import pyplot

class MOSFit3(MOS_FitProject):
    def __init__(self):
        super(MOSFit3, self).__init__()

        print 'loading data files...'
        self.loadAuroraFile('data/0p35.dat')
        print self.datasets.keys()
        print 'done.'

    def step1(self, param):
        fit = MOS_IV_Fit(param, {'VTH0':1.0, 'U0':600, 'UA':-1e-9, 'UB':-1e-18})

        dataset = self.datasets['n5x10.gtr']
        curve   = dataset.getCurve('Vgs', {'Vbs':0., 'Vds':0.1}, 'Id')
        fit.addDataSource(dataset.mosID(), curve)

        result,err = fit.doFit(plot=pyplot)

        for k,v in result.iteritems():
            param[k] = v
        return param, result, err

    def step2(self, param):
        fit = MOS_IV_Fit(param, {'K1':0.5, 'K2':0.1, 'U0':param['U0'], 'UA':param['UA'], 'UB':param['UB'], 'UC':0.0})

        dataset = self.datasets['n5x10.gtr']
        fit.addDataSource(dataset.mosID(), dataset.getCurve('Vgs', {'Vbs':0.,   'Vds':0.1}, 'Id'))
        fit.addDataSource(dataset.mosID(), dataset.getCurve('Vgs', {'Vbs':-1.0, 'Vds':0.1}, 'Id'))
        fit.addDataSource(dataset.mosID(), dataset.getCurve('Vgs', {'Vbs':-2.0, 'Vds':0.1}, 'Id'))
        fit.addDataSource(dataset.mosID(), dataset.getCurve('Vgs', {'Vbs':-3.0, 'Vds':0.1}, 'Id'))

        result,err = fit.doFit(plot=pyplot)

        for k,v in result.iteritems():
            param[k] = v
        return param, result, err

#    def step4(self, param):
#        fit = MOS_IV_Fit(param, {'K3':80, 'W0':2.5e-6, 'WINT':0., 'DWG':0.})
#
#        dataset = self.datasets['n5x10.gtr']
#        fit.addDataSource(dataset.mosID(), dataset.getCurve('Vgs', {'Vbs':0.,   'Vds':0.1}, 'Id'))
#
#        dataset = self.datasets['n0p6x10.gtr']
#        fit.addDataSource(dataset.mosID(), dataset.getCurve('Vgs', {'Vbs':0.,   'Vds':0.1}, 'Id'))
#
#        result,err = fit.doFit(plot=pyplot)
#
#        for k,v in result.iteritems():
#            param[k] = v
#        return param, result, err
#

param0 = {
    'TOX':7.2132e-9,
    'NCH':1.7e17,
    'T': 25.0,
    }

proj = MOSFit3()
proj.run(param0)
pyplot.show()

