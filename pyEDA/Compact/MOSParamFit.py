__all__=['MOS_IV_FitData', 'MOS_IV_Fit', 'MOS_FitProject']

from pyEDA.PDE.AutoDeriv import *
from pyEDA.Circuit.MOS3 import *
from pyEDA.Circuit.BSIM3v3 import *
from pyEDA.Compact.BSPData import *
from pyEDA.Compact.AuroraData import *

import numpy as np
from scipy.linalg import norm
from scipy.optimize import leastsq
from openopt import NLLSP
from matplotlib import pyplot
import pickle
import os
import string

def _validFileName(str):
    valid_chars = "-_.() %s%s" % (string.ascii_letters, string.digits)
    return ''.join(c for c in str if c in valid_chars).replace(' ','-')


class MOS_IV_FitData(object):
    def __init__(self, name=''):
        # curves is a list of (mosID, MOS_IV_Curve, weight_multiplier)
        self.curves = []
        self.subVth = False
        self.curr_th=1e-6
        self.curr_min=1e-12
        self.name = name

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

    def copy(self):
        res = MOS_IV_FitData()
        res.curves.extend(self.curves)
        res.subVth   = self.subVth
        res.curr_th  = self.curr_th
        res.curr_min = self.curr_min
        res.name     = self.name
        return res

    def __add__(self, other):
        return MOS_IV_FitData_Group() + self + other

    def __rmul__(self, other):
        ''' right-multiplied by a scalar'''
        s = float(other)
        res = MOS_IV_FitData()

        for mosID, curve, wMult in self.curves:
            res.curves.append( (mosID, curve, wMult*s) )
        
        res.subVth   = self.subVth
        res.curr_th  = self.curr_th
        res.curr_min = self.curr_min
        res.name     = self.name
        return res

    def setSubVth(self, subVth):
        self.subVth = subVth

    def plotData(self, modelCalc=None, title='', timeout=-1):
        '''
        @param modelCalc  a function fn(mosID, IOut, vbias) to calculate mos current
        '''
        fig = pyplot.figure()
        for mosID, curve, _ in self.curves:
            dataVScan = curve.dataVScan()
            dataCurr  = curve.dataCurr()
            dataModel = np.zeros(np.shape(dataCurr))

            if modelCalc:
                for i,v in enumerate(dataVScan):
                    vbias = curve.makeVBias(v)
                    dataModel[i] = modelCalc(mosID, curve.IOut, vbias)
                if self.subVth:
                    pyplot.semilogy(dataVScan, dataCurr, '+', dataVScan, dataModel, '-')
                    pyplot.ylim(ymin=0.1*self.curr_min)
                else:
                    pyplot.plot(dataVScan, dataCurr, '+', dataVScan, dataModel, '-')
            else:
                if self.subVth:
                    pyplot.semilogy(dataVScan, dataCurr, '+')
                else:
                    pyplot.plot(dataVScan, dataCurr, '+')

            pyplot.xlabel(curve.VScan)
            pyplot.ylabel(curve.IOut)

        title = '%s : %s' % (title,self.name)
        pyplot.title(title)
        fig.show()
        if timeout>0.0:
            pyplot.waitforbuttonpress(timeout)
        else:
            fig.savefig('output/' + _validFileName(title+'.png'), dpi=300)
            pyplot.close(fig)


class MOS_IV_FitData_Group(object):
    def __init__(self):
        self.srcs = []

    def __len__(self):
        cnt = 0
        for src in self.srcs:
            cnt = cnt + len(src)
        return cnt

    def copy(self):
        res = MOS_IV_FitData_Group()
        for src in self.srcs:
            res = res + src.copy()
        return res

    def mosIDList(self):
        list = []
        for src in self.srcs:
            list.extend(src.mosIDList())
        return list


    def __add__(self, other):
        res = MOS_IV_FitData_Group()
        if isinstance(other, MOS_IV_FitData_Group):
            res.srcs.extend(self.srcs)
            res.srcs.extend(other.srcs)
        elif isinstance(other, MOS_IV_FitData):
            res.srcs.extend(self.srcs)
            res.srcs.append(other)
        
        return res

    def __rmul__(self, other):
        ''' right-multiplied by a scalar'''
        s = float(other)
        res = MOS_IV_FitData_Group()
        for src in self.srcs:
            res.srcs.append(s*src)
        return res

    def setSubVth(self, subVth):
        for src in self.srcs:
            src.setSubVth(subVth)

    def iterData(self):
        for src in self.srcs:
            for t in src.iterData():
                yield t

    def plotData(self, modelCalc=None, name='', timeout=-1):
        for src in self.srcs:
            src.plotData(modelCalc, name, timeout)


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
        self.solver='nlp:scipy_tnc'

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

    def doFit(self):
        guess=[]
        scale=[] # scaling factor of variables
        ub=[] # upper bound
        lb=[] # lower bound
        for k in self.fitParam:
            v = self.baseParam[k]
            if isinstance(v, tuple):
                vmin = v[1]
                vmax = v[2]
                invS = abs(vmax-vmin)
                v=v[0]
            else:
                vmin = -1e100
                vmax = 1e100
                invS = min(1e5, max(1e-6, abs(float(v))))

            guess.append(v)
            lb.append(vmin)
            ub.append(vmax)
            scale.append(1./invS)

        #result, success = leastsq(self.fun, guess, (), self.jac, warning=True, factor=1., diag=scale)
        #result, cov_x, infodict, mesg, success = leastsq(self.fun, guess, (), self.jac, warning=True, factor=1., full_output=1)
        #print '*****', success, infodict['nfev'], mesg

        LSP = NLLSP(self.fun, guess, df = self.jac, ub=ub, lb=lb, scale=scale, xtol = 1e-12, ftol = 1e-12, gtol=1e-12, maxFunEvals = 100)
        r = LSP.solve(self.solver)
        result = r.xf
        if isinstance(result, np.float):
            result = [result]
        e = self.fun(result)

        param = {}
        for i,k in enumerate(self.fitParam):
            param[k] = result[i]
        return param, norm(e)

    def visualize(self, param, timeout=-1):
        paramVec = []
        for k in self.fitParam:
            paramVec.append(param[k])

        self.dataSrc.plotData(lambda mosID, IOut, vbias: self._modelCalc(mosID, IOut, paramVec, vbias),
                              self.name, timeout
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
        savefn = lambda x : 'output/%s.save.%d' % (self.name, x)
      
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

        
        print 'final params:\n%s' % self.paramToStr(self.param)

    def acceptParam(self, fitResult, keys):
        for k in keys:
            old = self.param[k]
            if isinstance(old, tuple):
                new = (fitResult[k], old[1], old[2])
            else:
                new = fitResult[k]
            self.param[k] = new


