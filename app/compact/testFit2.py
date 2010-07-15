from pyEDA.Compact.MOSParamFit import *
from pyEDA.Compact.BSPData import *
import math
from matplotlib import pyplot

class MOSFit2(MOS_FitProject):
    def __init__(self):
        super(MOSFit2, self).__init__()

        print 'loading data files...'
        self.loadBSimProFile('NEA40X10', 'data/NEA40X10.DAT')
        self.loadBSimProFile('NEA40X3',  'data/NEA40X3.DAT')
        self.loadBSimProFile('NEA1X10',  'data/NEA1X10.DAT')
        self.loadBSimProFile('NED40X1',  'data/NED40X1.DAT')
        print 'done.'

    def step1(self, param):
        fit = MOS_IV_Fit(param, {'VT0':1.0, 'U0':800, 'THETA':0.1})

        dataset = self.datasets['NEA40X10']
        curve   = dataset.getIVSeries('Id_Vg')[ {'Vbs':0, 'Vds':0.1} ]
        fit.addDataSource(dataset.mosID(), curve)

        result,err = fit.doFit(plot=pyplot)

        for k,v in result.iteritems():
            param[k] = v
        return param, result, err

    def step2(self, param):
        fit = MOS_IV_Fit(param, {'VT0':param['VT0'], 'LD':0., 'THETA':param['THETA']})

        dataset = self.datasets['NEA40X3']
        curve   = dataset.getIVSeries('Id_Vg')[ {'Vbs':0, 'Vds':0.1} ]
        fit.addDataSource(dataset.mosID(), curve)

        result,err = fit.doFit(plot=pyplot)

        param['LD'] = result['LD']
        return param, result, err

    def step3(self, param):
        fit = MOS_IV_Fit(param, {'VT0':param['VT0'], 'WD':0., 'THETA':param['THETA']})

        dataset = self.datasets['NEA1X10']
        curve   = dataset.getIVSeries('Id_Vg')[ {'Vbs':0, 'Vds':0.1} ]
        fit.addDataSource(dataset.mosID(), curve)

        result,err = fit.doFit(plot=pyplot)

        param['WD'] = result['WD']
        return param, result, err

    def step4(self, param):
        fit = MOS_IV_Fit(param, {'VT0':param['VT0'], 'RSDW':0., 'THETA':param['THETA']})

        dataset = self.datasets['NEA40X10']
        curve   = dataset.getIVSeries('Id_Vg')[ {'Vbs':0, 'Vds':0.1} ]
        fit.addDataSource(dataset.mosID(), curve)

        dataset = self.datasets['NEA40X3']
        curve   = dataset.getIVSeries('Id_Vg')[ {'Vbs':0, 'Vds':0.1} ]
        fit.addDataSource(dataset.mosID(), curve)

        result,err = fit.doFit(plot=pyplot)

        param['RSDW'] = result['RSDW']
        param['THETA'] = result['THETA']
        return param, result, err

    def step5(self, param):
        fit = MOS_IV_Fit(param, {'DELTA':0.0, 'THETA':param['THETA']})

        dataset = self.datasets['NEA1X10']
        curve   = dataset.getIVSeries('Id_Vg')[ {'Vbs':0, 'Vds':0.1} ]
        fit.addDataSource(dataset.mosID(), curve)

        result,err = fit.doFit(plot=pyplot)

        param['DELTA'] = result['DELTA']
        return param, result, err

    def step6(self, param):
        fit = MOS_IV_Fit(param, {'NSUB':param['NSUB']})

        dataset = self.datasets['NEA40X10']
        series  = dataset.getIVSeries('Id_Vg')
        fit.addDataSource(dataset.mosID(), series[ {'Vbs':0,    'Vds':0.1} ])
        fit.addDataSource(dataset.mosID(), series[ {'Vbs':-2.5, 'Vds':0.1} ])
        fit.addDataSource(dataset.mosID(), series[ {'Vbs':-5.0, 'Vds':0.1} ])

        result,err = fit.doFit(plot=pyplot)

        param['NSUB'] = result['NSUB']
        return param, result, err

    def step7(self, param):
        fit = MOS_IV_Fit(param, {'NSUB':param['NSUB'], 'XJ':param['XJ']})

        dataset = self.datasets['NEA40X10']
        series  = dataset.getIVSeries('Id_Vg')
        fit.addDataSource(dataset.mosID(), series[ {'Vbs':0,    'Vds':0.1} ])
        fit.addDataSource(dataset.mosID(), series[ {'Vbs':-2.5, 'Vds':0.1} ])
        fit.addDataSource(dataset.mosID(), series[ {'Vbs':-5.0, 'Vds':0.1} ])

        dataset = self.datasets['NEA40X3']
        series  = dataset.getIVSeries('Id_Vg')
        fit.addDataSource(dataset.mosID(), series[ {'Vbs':0,    'Vds':0.1} ])
        fit.addDataSource(dataset.mosID(), series[ {'Vbs':-2.5, 'Vds':0.1} ])
        fit.addDataSource(dataset.mosID(), series[ {'Vbs':-5.0, 'Vds':0.1} ])

        result,err = fit.doFit(plot=pyplot)

        param['NSUB'] = result['NSUB']
        param['XJ'] = result['XJ']
        return param, result, err

    def step9(self, param):
        fit = MOS_IV_Fit(param, {'VMAX':1e5, 'KAPPA':0.2, 'ETA':0.05})

        vgss= [2.034, 3.022, 4.010, 4.9980]
        dataset = self.datasets['NEA40X3']
        series  = dataset.getIVSeries('id_vd')
        for vgs in vgss: 
            fit.addDataSource(dataset.mosID(), series[ {'Vgs':vgs,  'Vbs':0.0} ])
            fit.addDataSource(dataset.mosID(), series[ {'Vgs':vgs,  'Vbs':-5.0} ])

        result,err = fit.doFit(plot=pyplot)

        param['VMAX']  = result['VMAX']
        param['KAPPA'] = result['KAPPA']
        param['ETA']   = result['ETA']
        return param, result, err

param0 = {
    'TOX':20e-9,
    'NSUB':1e16,
    'XJ': 50e-9
    }

proj = MOSFit2()
proj.run(param0)
pyplot.show()

