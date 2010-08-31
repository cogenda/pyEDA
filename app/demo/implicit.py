'''
Created on Jun 12, 2009

@author: hash
'''
from pyEDA.PDE.AutoDeriv import *
from pyEDA.PDE.ImplDeriv import *
import numpy as np
import math

class PolarCoord(ImplDeriv):
    def __init__(self):
        # we have 4 variables (rho, theta, x, y).
        # the first 2 are dependent variables
        # the last 2 are truly independent variables
        super(PolarCoord,self).__init__(2, 2)

    def initGuess(self):
        super(PolarCoord, self).initGuess()
        self.state.setVar(0, 2.)
        self.state.setVar(1, -1.)

    def calcFunJac(self):
        super(PolarCoord,self).calcFunJac() # must first call parent's version of this method
        rho,theta,x,y = self.state.getVars(range(4))
        self.state.setFunJac(self.sizeP+0, x*x+y*y-rho*rho)
        self.state.setFunJac(self.sizeP+1, rho * cos(theta) -x)

eqn = PolarCoord()
x = ADVar(1, 0)
y = ADVar(-1, 2)
eqn.setIndepVars([x,y]) # set the values of the truly independent variables
eqn.initGuess()
eqn.solve()
print '***********'
res = eqn.getDeriv()
for v in res:
  print v


