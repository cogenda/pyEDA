__all__=['IV_Curve', 'MOSFET_IV_Curve', 'MOSFET_Instance']

import numpy as np

class IV_Curve(object):
    Terminals = [ ('V1', 'I1'), ('V2', 'I2') ]
    def __init__(self, VScan, VConsts, IOut, data=None):
        '''
        VScan is the name of the scanning voltage
        VConsts is a map of terminal voltage name:value
        IOut is the name of the output current
        data is a 2D numarray, each row for a data point in the voltage scan
        '''
        self.VScan = VScan
        self.VConsts = VConsts
        self.IOut = IOut
        self.data = data

        self._vbiasBuf = [0.0] * len(self.Terminals)
        i, idxVScan = (0, -1)
        for vname, iname in self.Terminals:
            if vname==self.VScan:
                self._idxVScan = i
            else:
                self._vbiasBuf[i] = self.VConsts[vname]
            i = i+1

    def __len__(self):
        return np.shape(self.data)[0]
    
    def dataVScan(self):
        return self.data[:,0]

    def dataCurr(self):
        return self.data[:,1]

    def makeVBias(self, scanVal):
        self._vbiasBuf[self._idxVScan] = scanVal
        return tuple(self._vbiasBuf)

    def iterData(self):
        l = np.shape(self.data)[0]

        for i in range(l):
            volt, curr = self.data[i]
            vbias = self.makeVBias(volt)
            yield (vbias, curr)

class CurveSeries(object):
    def __init__(self, label):
        self.label = label
        self.curves = []
        self.VScan = None
        self.IOut  = None

    def __len__(self):
        return len(self.curves)

    def __getitem__(self, key):
        '''
        get a curve from the series,

        @param key      index of the curve or a filter function.
                        In the second case, each curve c is passed to the filter fn(c),
                        if fn() yields True, c is returned.
        '''
        if isinstance(key, int):
            idx = key
            if idx>=0 and idx<len(self):
                return self.curves[idx]
            else:
                raise IndexError
        elif isinstance(key, dict):
            cond = key
            for c in self.curves:
                flg = True
                for n,v in cond.iteritems():
                    if abs(c.VConsts[n]-v)>1e-8:
                        flg = False
                if flg:
                    return c
        elif callable(key):
            fn = key
            for c in self.curves:
                if fn(c):
                    return c
        return None

    def __setitem__(self, key, curve):
        '''
        '''
        if not isinstance(key, int):
            raise TypeError
        if not isinstance(curve, IV_Curve):
            raise TypeError

        if self.VScan and self.IOut:
            if not (curve.VScan==self.VScan and curve.IOut==self.IOut):
                raise ValueError
        else:
            self.VScan = curve.VScan
            self.IOut  = curve.IOut

        if key==len(self):
            self.curves.append(curve)
        elif key>=0 and key<len(self):
            self.curves[key] = curve
        else:
            raise IndexError

    def append(self, curve):
        self[len(self)] = curve


class MOSFET_IV_Curve(IV_Curve):
    Terminals = [ ('Vgs', 'Ig'),
                  ('Vds', 'Id'),
                  ('Vbs', 'Isub') ]

class MOSFET_Instance(object):
    def __init__(self, W=None, L=None, T=None):
        self.W = W
        self.L = L
        self.T = T
        self.IVData = {}

    def mosID(self):
        return (self.W, self.L, self.T)

    def addIVSeries(self, label, series):
        self.IVData[label] = series

    def getIVSeries(self, label):
        return self.IVData[label]

    def listIVSeries(self):
        return self.IVData.keys()

    def loadBSimProFile(self, fname):
        from ParserBSP import makeInstance
        makeInstance(self, fname)

if __name__=='__main__':
    mos = MOSFET_Instance()
    mos.loadBSimProFile('data/NEA1X10.DAT')
    
    for sLabel in mos.listIVSeries():
        print 'Series:', sLabel
        series = mos.getIVSeries(sLabel)
        for i in xrange(len(series)):
            curve = series[i]

            print 'Curve --------------',len(curve)
            for vbias, curr in curve.iterData():
                print vbias, curr


