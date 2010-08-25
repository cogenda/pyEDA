from pyEDA.PDE.NLEqns import *
from pyEDA.PDE.AutoDeriv import *
import numpy as np
import scipy

class SimpleEqns(NLEqns):
    ''' Solution class
    '''

    def __init__(self):
        super(SimpleEqns,self).__init__()
        self.state = NLEqnState(2) # prepare solution state data for 2 unknowns
	print '''Setting up the equations
	----------------------------------
	    F(x,y)=0

                {  f0 = x*x + y*y - 4.0
	    F = {
                {  f1 = x*x - y*y - 1.0
	----------------------------------'''

    def initGuess(self):
        # initial guess
        self.state.x = np.array([1, 1])

    def calcFunJac(self):
        x, y = self.state.getVars([0,1])
        
        # compute function, the derivative information is contained in f0 and f1
        f0 = x*x + y*y - 4.0
        f1 = x*x - y*y - 1.0
        # print f0, f1      # to see the variable and derivatives

        # set function value as well as jacobian matrix elements
        self.state.setFunJac(0, f0)
        self.state.setFunJac(1, f1)

# main program

eqn = SimpleEqns()
eqn.initGuess()
eqn.solve()

print ''
print ' ------------ Result ---------------'
print eqn.state.x
