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


class MOSFET_IV_Curve(IV_Curve):
    Terminals = [ ('Vgs', 'Ig'),
                  ('Vds', 'Id'),
                  ('Vbs', 'Isub') ]

class MOSFET_Instance(object):
    def __init__(self, W=None, L=None, T=None):
        self.W = W
        self.L = L
        self.T = T

    def mosID(self):
        return (self.W, self.L, self.T)

