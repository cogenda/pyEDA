'''
Created on Jun 12, 2009

@author: hash
'''
from PDE.AutoDeriv import *
from PDE.ImplDeriv import *
import numpy as np
import math

class PolarCoord(ImplDeriv):
    def __init__(self):
        super(PolarCoord,self).__init__(2, 2)

    def initGuess(self):
        self.state.setVec([2, -1, 1, -1])

    def calcFunJac(self):
        super(PolarCoord,self).calcFunJac()
        rho,theta,x,y = self.state.getVars(range(4))
        self.state.setFunJac(self.sizeP+0, x*x+y*y-rho*rho)
        self.state.setFunJac(self.sizeP+1, rho * cos(theta) -x)

eqn = PolarCoord()
x = ADVar(1, 0)
y = ADVar(1, 2)
eqn.setIndepVars([x,y])
eqn.initGuess()
eqn.solve()
print '***********'
res = eqn.getDeriv()
for v in res:
  print v


