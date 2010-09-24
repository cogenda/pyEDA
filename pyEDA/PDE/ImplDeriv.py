'''
Automatic Differetiation for Implicitly Defined Functions


'''

__all__ = ['ImplDeriv']

from AutoDeriv import *
from NLEqns import *
from scipy import linalg
import numpy as np

class ImplDeriv(NLEqns):
    def __init__(self, m, p):
        '''m dependent variables (y) are implicitly defined as functions 
           of p truely independent variables (x).
           [ y_0, y_1, y_2, ... , y_m-1,  x_0, x_1, ... , x_p-1 ] 
           '''
        super(ImplDeriv,self).__init__()
        self.sizeM = m
        self.sizeP = p
        self.iVars = None  # independent variables
        n = m + p
        #if n > 20:  # we only plan to handle small implicit problems
        #    raise ValueError
        self.state = NLEqnState(n)
        self.Jinv = None

    def setIndepVars(self, vars):
        '''
        set the values of the independent variables.
        @param: vars     list of ADVar
        '''
        if not isinstance(vars, list):
            raise TypeError
        if not len(vars)==self.sizeP:
            raise ValueError
        if not isinstance(vars[0], ADVar):
            raise TypeError

        self.iVars = vars

    def getDeriv(self, vars=None):
        '''
        Note: call solve() before calling this function.
        @param: vars list of dependent variable indeces
        @return: dependent variables as list of ADVar
        '''
        res = []

        if self.Jinv==None:
            self.Jinv = linalg.inv(self.state.J.todense())
        Jinv=self.Jinv

        if vars==None:
            rows=range(self.sizeM)
        elif isinstance(vars, list): 
            rows=vars
        elif isinstance(vars, int): 
            rows=[vars]
        else:
            raise TypeError

        for i in rows:
            v = ADVar()
            for j in xrange(self.sizeP):
                v += Jinv[i,j] * self.iVars[j]
            v.setVal(self.state.x[i])
            res.append(v)
        if len(res)==1:
            return res[0]
        return res

    def initGuess(self):
        for i in xrange(self.sizeP):
            p = self.sizeM + i
            self.state.setVar(p, self.iVars[i].getVal())

    def calcFunJac(self):
        for i in xrange(self.sizeP):
            p = self.sizeM + i
            self.state.setFunJac(i, self.state.getVar(p) - self.iVars[i].getVal())

    def solve(self):
        self.Jinv = None
        super(ImplDeriv,self).solve()

