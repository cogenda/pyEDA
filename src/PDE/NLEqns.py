'''
Created on Jun 12, 2009

@author: hash
'''
__all__ = ['NLEqnState', 'NLEqns', 'OPADD', 'OPSET']

import numpy as np
import scipy
from scipy import sparse
from scipy import linalg
from scipy.sparse import linalg as splinalg
from scipy.sparse.linalg import dsolve

from AutoDeriv import *

OPADD = 0
OPSET = 1

class NLEqnState(object):
    def __init__(self, n=0):
        self.N = n
        self.x = None
        self.J = None
        self.dx = None
        self.clock = None
        self.NTStep = 2 # past time step to keep 
        self.ptime = None
        self.px = None
        
        if n>0:
            self._prepareData()

    def _prepareData(self):
        self.x = scipy.zeros(self.N)
        self.b = scipy.zeros(self.N)
        self.J = sparse.lil_matrix((self.N, self.N))
        self.dx = None
        self.clock = 0
        self.ptime = []
        self.px = []
        
    def size(self):
        return self.nn
    
    def setSize(self, n):
        if not n>0:
            raise ValueError
        self.N = n
        _prepareData()

    def getVar(self, idx):
        if not idx<self.N:
            raise IndexError
        return ADVar(self.x[idx], idx)
    
    def getVars(self, indices):
        vars = []
        for idx in indices:
            vars.append(self.getVar(idx))
        return vars
    
    def getTimeDeriv(self, idx):
        if self.clock==0.0:
            return 0.0
        if not idx<self.N:
            raise IndexError
        if not len(self.ptime)>0:
            raise IndexError
        
        x1 = self.px[0]
        dt1 = self.clock - self.ptime[0]
        if dt1<1e-15:
            return 0
        else:
            return (self.getVar(idx)-x1[idx])/dt1

    def getTimeDerivs(self, indices):
        ret = []
        for idx in indices:
            ret.append(self.getTimeDeriv(idx))
        return ret
   
    def saveTimeStep(self):
        self.ptime.insert(0, self.clock)
        self.px.insert(0, self.x)
        
        nstep = len(self.ptime) 
        if nstep > self.NTStep:
            del self.ptime[nstep-1]
            del self.px[nstep-1]
        
    def advanceClock(self, dt):
         self.clock += dt
    
    def setVar(self, idx, v=0):
        if not idx<self.N:
            raise IndexError
        self.x[idx] = v
    
    def setVec(self, vec):
        self.x = vec
    
    def setFunJac(self, idx, advar, op=OPADD):
        if not idx<self.N:
            raise IndexError
        if not isinstance(advar, ADVar):
            raise TypeError
                
        if op==OPADD:
            self.b[idx] += advar.val
        else:
            self.b[idx] = advar.val
            
        deriv = advar.getDeriv()
        for ix,dx in deriv:
            if op==OPADD:
                self.J[idx,ix] += dx
            else:   
                self.J[idx,ix] = dx

    def resetEqn(self, idx):
        self.b[idx] = 0.0
        self.J[idx,:] = 0.0
        
    def connectVar(self, idx1, idx2, advar=None, op=OPADD):
        self.J[idx1,:] += self.J[idx2,:]
        self.b[idx1] += self.b[idx2] 
        
        self.J[idx2,:] = 0.0
        self.J[idx2,idx1] = -1.0
        self.J[idx2,idx2] = 1.0
        self.b[idx2] = 0.0

    def clearFunJac(self):
        self.b = scipy.zeros(self.N)
        self.J = sparse.lil_matrix((self.N, self.N))
        
class NLEqns(object):
    def __init__(self):
        self.state = None

    def calcFunJac(self):
        pass
    
    def initGuess(self):
        pass
        
    def checkConv(self):
        if self.state.N==0:
            return (True,0)
        
        res = self.state.b
        norm = linalg.norm(res)/self.state.N
        return (norm<1e-8, norm)

    def dampStep(self, dx):
        return dx
    
    def solve(self):
        maxiter=10
        trace = True
        for iter in xrange(0,maxiter):
            self.state.clearFunJac()
            self.calcFunJac()
            flagConv, err = self.checkConv()
            if trace:
                print iter, err
            if flagConv:
                break

            dx = np.negative(dsolve.spsolve(self.state.J.tocsr(), self.state.b))
            self.state.dx = self.dampStep(dx)

            self.state.x = np.add(self.state.x, self.state.dx)

        #x = splinalg.bicg(eqns.J, eqns.b, None, 1e-6)
