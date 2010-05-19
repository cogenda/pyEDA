'''
Automatic Differetiation for Implicitly Defined Functions


'''

__all__ = ['ImplDeriv']

from PDE.AutoDeriv import *
from PDE.NLEqns import *
from scipy import linalg

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
        if n > 20:  # we only plan to handle small implicit problems
            raise ValueError
        self.state = NLEqnState(n)

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

    def getDeriv(self):
        '''
        call solve() before calling this function.
        @return: dependent variables as list of ADVar
        '''
        res = self.sizeM * [ADVar()]
        Jinv = linalg.inv(self.state.J.todense())
        for i in xrange(self.sizeM):
            for j in xrange(self.sizeP):
                res[i] += Jinv[i,j] * self.iVars[j]
            res[i].setVal(self.state.x[i])
        return res

    def calcFunJac(self):
        for i in xrange(self.sizeP):
            p = self.sizeM + i
            self.state.setFunJac(i, self.state.getVar(p) - self.iVars[i].getVal())


