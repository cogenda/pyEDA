__all__=['MOS_IV_Fit', 'MOS_FitProject']

from pyEDA.PDE.AutoDeriv import *
from pyEDA.Circuit.MOS3 import *
from pyEDA.Compact.BSPData import *

import numpy as np
from scipy.linalg import norm
from scipy.optimize import leastsq
from matplotlib import pyplot

class MOS_IV_FitData(object):
    def __init__(self):
        # curves is a list of (mosID, MOS_IV_Curve)
        self.curves = []

    def __len__(self):
        cnt = 0
        for mosID, curve in self.curves:
            cnt = cnt + len(curve)
        return cnt

    def addCurve(self, mosID, curve):
        self.curves.append( (mosID,curve) )

    def iterData(self):
        for mosID, curve in self.curves:
            IOut = curve.IOut
            for vbias, curr in curve.iterData():
                weight = 1.0
                yield (mosID, IOut, vbias, curr, weight)

    def plotData(self, plt, modelCalc=None):
        '''
        @param modelCalc  a function fn(mosID, IOut, vbias) to calculate mos current
        '''
        for mosID, curve in self.curves:
            dataVScan = curve.dataVScan()
            dataCurr  = curve.dataCurr()
            dataModel = np.zeros(np.shape(dataCurr))

            if modelCalc:
                for i,v in enumerate(dataVScan):
                    vbias = curve.makeVBias(v)
                    dataModel[i] = modelCalc(mosID, curve.IOut, vbias)
                plt.plot(dataVScan, dataCurr, '+', dataVScan, dataModel, '-')
            else:
                plt.plot(dataVScan, dataCurr, '+')


class CachedMOS(object):
    def __init__(self, baseParam, fitParam):
        '''
        @param baseParam   MOS parameters (name:value map)
        @param fitParam    parameters to fit (list of names)
        '''
        self.baseParam = baseParam
        self.fitParam = fitParam
        self.param = None
        self.cache = {}
        self.mos = None

    def initMOS(self, param):
        '''
        @param param        fitting parameters, numpy array, same order as in fitParam
        '''
        if not self.mos == None  and self.fitParam ==None and np.allclose(param, self.param): 
            return  # no need to init

        if isinstance(param, np.float):
            self.param = np.array([param])
        else:
            self.param = param
        self.cache = {}

        MOS_param = self.baseParam.copy()
        for i,name in enumerate(self.fitParam):
            MOS_param[name] = ADVar(self.param[i], i)

        self.mos = MOSLv3(**MOS_param)

    def ids(self, param, vbias):
        self.initMOS(param)
        if self.cache.has_key(vbias):
            return self.cache[vbias]

        Vgs, Vds, Vbs = vbias
        ids = self.mos._Ids(Vgs, Vds, Vbs)
        self.cache[vbias] = ids
        return ids

class MOS_IV_Fit(object):
    def __init__(self, baseParam, fitParam):
        '''
        @param baseParam   MOS parameters (name : value map)
        @param fitParam    parameters to fit (name : init value map)
        '''
        self.baseParam = baseParam
        self.fitParam  = fitParam
        self.mos = {}  # map of mosID:CachedMOS
        self.dataSrc = MOS_IV_FitData()

    def addDataSource(self, mosID, curve):
        '''
        @param mosID    (W,L,T) tuple
        '''
        MOS_param = self.baseParam.copy()
        MOS_param['W'] = mosID[0] * 1e-6
        MOS_param['L'] = mosID[1] * 1e-6
        MOS_param['T'] = mosID[2] + 273.15

        if not self.mos.has_key(mosID):
            self.mos[mosID] = CachedMOS(MOS_param, self.fitParam.keys())
        self.dataSrc.addCurve(mosID, curve)

    def _modelCalc(self, mosID, IOut, param, vbias):
        mos = self.mos[mosID]
        if IOut=='Id':
            return mos.ids(param, vbias)
        else:
            pass # not supported yet
        return None

    def fun(self, param):
        res = np.zeros(len(self.dataSrc))

        i = 0
        for mosID, IOut, vbias,curr,weight in self.dataSrc.iterData():
            curr_m = self._modelCalc(mosID, IOut, param, vbias)

            res[i] = (curr_m - curr) * weight
            i += 1
        return res

    def jac(self, param):
        res = np.zeros([len(self.dataSrc), len(self.fitParam)])

        i = 0
        for mosID, IOut, vbias,curr,weight in self.dataSrc.iterData():
            curr_deriv = self._modelCalc(mosID, IOut, param, vbias).getDeriv()

            for iVar, pDeriv in curr_deriv:
                res[i,iVar] = pDeriv*weight
            i += 1
        return res

    def doFit(self, plot=None):
        guess = self.fitParam.values()

        result, success = leastsq(self.fun, guess, (), self.jac)
        if isinstance(result, np.float):
            result = [result]
        e = self.fun(result)

        if plot:
            self.visualize(result, plot)

        param = {}
        for i,k in enumerate(self.fitParam.keys()):
            param[k] = result[i]
        return param, norm(e)

    def visualize(self, param, plt=pyplot):
        plt.figure()
        self.dataSrc.plotData(plt, 
                              lambda mosID, IOut, vbias: self._modelCalc(mosID, IOut, param, vbias) 
                             )


class MOS_FitProject(object):
    def __init__(self):
        self.datasets = {}

    def loadBSimProFile(self, name, fname):
        ins = BSP_MOSFET_Instance()
        ins.loadBSimProFile(fname)
        self.datasets[name] = ins

    def paramToStr(self, param):
        i=0
        s=''
        for k,v in param.iteritems():
            s += '%8s: %12G' % (k,v)
            i+=1
            if i%4==0:
              s += '\n'
        return s

    def run(self, param0):
        i = 0
        param = param0

        for i in xrange(100):
            nameStep = 'step%d' % i
            if hasattr(self, nameStep):
                fnStep = self.__getattribute__(nameStep)
                print '========== Step %d ===========' % i
                print 'known params:\n%s' % self.paramToStr(param)

                param, result, err = fnStep(param)

                print 'fit error: %g\nfit result:\n%s' % (err, self.paramToStr(result))
                print 'accept params:\n%s' % self.paramToStr(param)
                print '\n'




