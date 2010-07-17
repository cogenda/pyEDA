__all__=['MOS_IV_FitData', 'MOS_IV_Fit', 'MOS_FitProject']

from pyEDA.PDE.AutoDeriv import *
from pyEDA.Circuit.MOS3 import *
from pyEDA.Circuit.BSIM3v3 import *
from pyEDA.Compact.BSPData import *
from pyEDA.Compact.AuroraData import *

import numpy as np
from scipy.linalg import norm
from scipy.optimize import leastsq
from matplotlib import pyplot
import pickle
import os

class MOS_IV_FitData(object):
    def __init__(self):
        # curves is a list of (mosID, MOS_IV_Curve, weight_multiplier)
        self.curves = []
        self.subVth = False
        self.curr_th=1e-6
        self.curr_min=1e-12

    def __len__(self):
        cnt = 0
        for mosID, curve, _ in self.curves:
            cnt = cnt + len(curve)
        return cnt

    def mosIDList(self):
        list = []
        for mosID, _, _ in self.curves:
            list.append(mosID)
        return list

    def addCurve(self, mosID, curve, wMult=1.0):
        '''
        @param mosID    (W,L,T) tuple
        @param curve    MOSFET_IV_Curve object
        @param wMult    weight multiplier
        '''
        self.curves.append( (mosID,curve, wMult) )

    def iterData(self):
        for mosID, curve, wMult in self.curves:
            IOut = curve.IOut
            for vbias, curr in curve.iterData():
                weight = 1.0
                if self.subVth:
                    weight = (self.curr_th+abs(curr))/max(abs(curr), self.curr_min)
                weight *= wMult
                yield (mosID, IOut, vbias, curr, weight)

    def __add__(self, other):
        res = MOS_IV_FitData()

        res.curves.extend(self.curves)
        res.curves.extend(other.curves)
        
        res.subVth   = self.subVth or other.subVth
        res.curr_th  = min(self.curr_th, other.curr_th)
        res.curr_min = max(self.curr_min, other.curr_min)
        return res

    def __rmul__(self, other):
        ''' right-multiplied by a scalar'''
        s = float(other)
        res = MOS_IV_FitData()

        for mosID, curve, wMult in self.curves:
            res.curves.append( (mosID, curve, wMult*s) )
        
        res.subVth   = self.subVth
        res.curr_th  = self.curr_th
        res.curr_min = self.curr_min
        return res

    def plotData(self, plt, modelCalc=None, name=''):
        '''
        @param modelCalc  a function fn(mosID, IOut, vbias) to calculate mos current
        '''
        for mosID, curve, wMult in self.curves:
            dataVScan = curve.dataVScan()
            dataCurr  = curve.dataCurr()*wMult
            dataModel = np.zeros(np.shape(dataCurr))

            if modelCalc:
                for i,v in enumerate(dataVScan):
                    vbias = curve.makeVBias(v)
                    dataModel[i] = modelCalc(mosID, curve.IOut, vbias)*wMult
                if self.subVth:
                    plt.semilogy(dataVScan, dataCurr, '+', dataVScan, dataModel, '-')
                    plt.ylim(ymin=0.1*self.curr_min)
                else:
                    plt.plot(dataVScan, dataCurr, '+', dataVScan, dataModel, '-')
            else:
                if self.subVth:
                    plt.semilogy(dataVScan, dataCurr, '+')
                else:
                    plt.plot(dataVScan, dataCurr, '+')
        plt.title(name)


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
        #if self.mos and self.fitParam and np.allclose(param, self.param):

        if self.mos and self.fitParam and np.array_equal(param, self.param):
            return  # no need to init

        if isinstance(param, np.float):
            self.param = np.array([param])
        else:
            self.param = np.copy(param)
        self.cache = {}

        MOS_param = self.baseParam.copy()
        for i,name in enumerate(self.fitParam):
            MOS_param[name] = ADVar(self.param[i], i)

        #self.mos = MOSLv3(**MOS_param)
        self.mos = MOSBSim3v3(**MOS_param)


    def ids(self, param, vbias):
        self.initMOS(param)
        if self.cache.has_key(vbias):
            return self.cache[vbias]

        Vgs, Vds, Vbs = vbias
        #ids = self.mos._Ids(Vgs, Vds, Vbs)
        #self.cache[vbias] = ids
        #return ids

        Id, Is, Isub  = self.mos._DC_Curr(Vgs,Vds,Vbs)
        self.cache[vbias] = Id
        return Id


class MOS_IV_Fit(object):
    def __init__(self, baseParam, fitParam, name=''):
        '''
        @param baseParam   MOS parameters (name : value map)
        @param fitParam    parameters to fit (name : init value map)
        '''
        self.baseParam = baseParam
        self.fitParam  = fitParam
        self.mos = {}  # map of mosID:CachedMOS
        self.dataSrc = MOS_IV_FitData()
        self.name = name

    def setDataSource(self, src):
        '''
        @param src  MOS_IV_FitData object
        '''
        self.dataSrc = src
        for mosID in src.mosIDList():
            MOS_param={}
            for k,v in self.baseParam.iteritems():
                if isinstance(v, tuple):
                    MOS_param[k] = v[0]
                else:
                    MOS_param[k] = v
    
            W,L,TEMP = [float(x) for x in mosID]
            if W>1e-3:  W*=1e-6
            if L>1e-3:  L*=1e-6
            MOS_param['W'] = W
            MOS_param['L'] = L
            MOS_param['TEMP'] = TEMP

    
            if not self.mos.has_key(mosID):
                self.mos[mosID] = CachedMOS(MOS_param, self.fitParam)


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
        guess=[]
        scale=[] # scaling factor of variables
        for k in self.fitParam:
            v = self.baseParam[k]
            if isinstance(v, tuple):
                invS = min(1e5, max(abs(v[1]), abs(v[2])))
                v=v[0]
            else:
                invS = min(1e5, max(1e-6, abs(float(v))))

            guess.append(v)
            scale.append(1./invS)

        #result, success = leastsq(self.fun, guess, (), self.jac, warning=True, factor=1., diag=scale)
        result, cov_x, infodict, mesg, success = leastsq(self.fun, guess, (), self.jac, warning=True, factor=1., full_output=1)
        print '*****', success, mesg
        if isinstance(result, np.float):
            result = [result]
        e = self.fun(result)

        if plot:
            self.visualize(result, plot)

        param = {}
        for i,k in enumerate(self.fitParam):
            param[k] = result[i]
        return param, norm(e)

    def visualize(self, param, plt=pyplot):
        plt.figure()
        self.dataSrc.plotData(plt, 
                              lambda mosID, IOut, vbias: self._modelCalc(mosID, IOut, param, vbias),
                              self.name
                             )


class MOS_FitProject(object):
    def __init__(self, param0):
        self.datasets = {}
        self.param = param0
        self.name = type(self).__name__

    def loadBSimProFile(self, name, fname):
        ins = BSP_MOSFET_Instance()
        ins.loadBSimProFile(fname)
        self.datasets[name] = ins

    def loadAuroraFile(self, fname):
        file = AuroraFile(fname)
        file.make_instances()
        inss = file.instances
        for name,ins in inss.iteritems():
            self.datasets[name] = ins

    def paramToStr(self, param):
        i=0
        s=''
        keys=param.keys()
        keys.sort()
        for k in keys:
            v = param[k]
            if isinstance(v, tuple):
                v = v[0]
            s += '%8s: %12G' % (k,v)
            i+=1
            if i%4==0:
              s += '\n'
        return s

    def run(self, start=0, stop=1000):
        savefn = lambda x : '%s.save.%d' % (self.name, x)
      
        if start>0:
            for i in xrange(start):
                fn = savefn(i)
                if os.path.exists(fn):
                    print 'loading previous parameters: %s' % fn
                    fload = open(fn)
                    self.param = pickle.load(fload)
                    fload.close()

        print 'known params:\n%s' % self.paramToStr(self.param)
        for i in xrange(start,stop+1):
            nameStep = 'step%d' % i
            if hasattr(self, nameStep):
                fnStep = self.__getattribute__(nameStep)
                print '========== Step %d ===========' % i

                result, err = fnStep()

                print 'fit error: %g\nfit result:\n%s' % (err, self.paramToStr(result))
                #print 'accept params:\n%s' % self.paramToStr(self.param)
                print '\n'

                # save params
                fsave = open(savefn(i), 'w')
                pickle.dump(self.param,fsave)
                fsave.close()

        
        #print 'final params:\n%s' % self.paramToStr(self.param)

    def acceptParam(self, fitResult, keys):
        for k in keys:
            old = self.param[k]
            if isinstance(old, tuple):
                new = (fitResult[k], old[1], old[2])
            else:
                new = fitResult[k]
            self.param[k] = new


