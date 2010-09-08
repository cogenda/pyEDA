'''
Nonlinear PDE Framework

A framework for solving nonlinear PDE of the form F(x)=0.
The class NLEqns should not be instantiated directly.
Instead, one extends NLEqns, and implements the C{calcFunJac()}
function that calculate the PDE function ad Jacobian.
One may want to also implement the following methods to 
influence various aspects of the solver.
    - calcFunJac()
    - initGuess()
    - checkConv()
    - dampStep()
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
    '''
    Current state of the PDE solver.
        - C{x}  solution vector of the current iteration
        - C{b}  residue of the current iteration
        - C{J}  Jacobian matrix of the current iteration
        - C{dx} solution update to be applied on the current iteration
        - C{clock}  simulation system time
        - C{NTStep} number of past time steps to be saved
        - C{ptime}  clock time of past saved steps
        - C{px} solution of past saved steps
        - C{integ} integration
    '''

    def __init__(self, n=0):
        '''
        constructor.
        @param: n   number of equations/variables
        '''
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
        '''
        initialize data structure
        '''
        self.x = scipy.zeros(self.N)

        self.dx = None
        self.clock = 0
        self.ptime = []
        self.px = []
        self.integ = scipy.zeros(self.N)

        self.clearFunJac()
        
    def size(self):
        '''Return the number of equations'''
        return self.N
    
    def setSize(self, n):
        '''Change the number of equations'''
        if not n>0:
            raise ValueError
        self.N = n
        _prepareData()

    def getVar(self, idx):
        '''
        Get the variable idx from the current iteration.

        @param: idx  variable index 
        @return: the AD variable at index idx.
        '''
        if not idx<self.N:
            raise IndexError
        return ADVar(self.x[idx], idx)
    
    def getVars(self, indices):
        '''
        Get a list of AD variables from the current iteration.

        @param: indices     list of variable indices
        '''
        vars = []
        for idx in indices:
            vars.append(self.getVar(idx))
        return vars
    
    def getTimeDeriv(self, idx):
        '''
        Get the time derivative of a variable, which is an AD variable

        @param:  idx     variable index
        '''
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
        '''
        Get the time derivative of a list of variables
        '''
        ret = []
        for idx in indices:
            ret.append(self.getTimeDeriv(idx))
        return ret
    
    def getTimeInteg(self, idx):
        '''
        Get the integration of a variable from time zero,
        or from the last time resetTimeInteg() is called.
        
        @param: idx   variable index 
        '''
        if self.clock==0.0:
            return 0.0
        if not idx<self.N:
            raise IndexError
        if not len(self.ptime)>0:
            raise IndexError
        
        x1 = self.px[0]
        dt1 = self.clock - self.ptime[0]
        return self.integ[idx] + 0.5*(self.getVar(idx)+x1[idx])*dt1
    
    def getTimeIntegs(self, indices):
        '''
        Get the integration of a list of variables
        '''
        ret = []
        for idx in indices:
            ret.append(self.getTimeDeriv(idx))
        return ret
    
    def resetTimeInteg(self):
        self.integ = scipy.zeros(self.N)
    
    def saveTimeStep(self):
        '''
        Save the current iteration as the accepted solution.

        This will be used for time-derivative-calculation in the future.
        '''
        # update integration
        if len(self.ptime)>0:
            dt = self.clock - self.ptime[0]
            self.integ += 0.5*(self.px[0]+self.x)*dt

        # update saved solution
        self.ptime.insert(0, self.clock)
        self.px.insert(0, self.x)
        
        # delete excessive old solution
        nstep = len(self.ptime) 
        if nstep > self.NTStep:
            del self.ptime[nstep-1]
            del self.px[nstep-1]
            
        
    def advanceClock(self, dt):
        ''' Advance the system clock by dt. '''
        self.clock += dt
    
    def setVar(self, idx, v=0):
        '''Set a variable in the current iteration solution

        @param: idx     variable index
        @param: v       variable value (normal scalar)
        '''
        if not idx<self.N:
            raise IndexError
        self.x[idx] = v
    
    def setVec(self, vec):
        '''Set the entire current solution vector'''
        if isinstance(vec, np.ndarray):
            if len(vec) == self.N:
                self.x = vec
            else:
                raise ValueError
        elif isinstance(vec, list):
            if len(vec) == self.N:
                self.x = np.array(vec)
            else:
                raise ValueError
        else:
            raise TypeError
            
    
    def setFunJac(self, idx, advar, op=OPADD):
        '''
        Evaluate the function of an equation.

        Set (or add to) the current iteration residue b.
        Also set (or add to) the Jacobian matrix.

        @param: idx     Equation(variable) index.
        @param: advar   The value of the function for this equation.
                        The numerical value of advar contributes to b at row idx,
                        whereas the partial derivatives contributes to the Jacobian at row idx.
        @param: op      if C{OPADD}, advar will be added to b and J at idx.
                        if C{OPSET}, advar will replace any previous b and J at idx.
        '''
        if not idx<self.N:
            raise IndexError
        if not isinstance(advar, ADVar):
            raise TypeError
                
        deriv = advar.getDeriv()
        if op==OPADD:
            self.b[idx] += advar.getVal()
            for ix,dx in deriv:
                self.Jrows[idx].append((ix,dx))
        else: #OPSET
            raise Exception # not implemented
            
    def resetEqn(self, idx):
        '''
        Clear the residue vector b and Jacobian matrix J at row idx.
        '''
        self.b[idx] = 0.0
        self.Jrows[idx] = []
        
    def connectVar(self, idx1, idx2, advar=None, op=OPADD):
        '''
        Connect the variable at idx1 and idx2.

        Three operations:
            - add row idx2 to row idx1
            - clear row idx2
            - J[idx2, idx1] = 1, J[idx2, idx2] = -1, b[idx2]=0

        This effectively force the two variables equal.
        '''
        self.b[idx1] += self.b[idx2] 
        self.b[idx2] = 0.0

        # add row idx2 to row idx1
        self.Jrows[idx1].extend(self.Jrows[idx2])
        # clear row idx2
        self.Jrows[idx2] = []

        # J[idx2, idx1]= 1
        self.Jrows[idx2].append((idx1,1.0))

        # J[idx2, idx2]=-1
        self.Jrows[idx2].append((idx2,-1.0))
 
    def clearFunJac(self):
        '''
        Clear all entries in Jacobian and residue.
        '''
        self.b = scipy.zeros(self.N)
        self.J = None # Jacobian Matrix not assembled
        self.Jrows = []
        for i in xrange(self.N):
            self.Jrows.append([]) # list of (col, val)


    def assembleJac(self):
        '''
        Assemble jacobian matrix. add duplicate items
        '''
        ii=[]
        jj=[]
        vv=[]
        for i,row in enumerate(self.Jrows):
            for j,v in row:
                ii.append(i)
                jj.append(j)
                vv.append(v)
        self.J=sparse.coo_matrix((vv,(ii,jj)), (self.N,self.N)).tocsr()
        
class NLEqns(object):
    def __init__(self):
        '''
        Constructor.

        Reimplementation of this method should setup mesh and other data.
        '''
        self.state = None

    def calcFunJac(self):
        '''
        Calculate the residue vector and Jacobian matrix using the current solution x 
        (in self.state).

        This method must be reimplemented for the NL PDE solver to function.
        '''
        pass
    
    def initGuess(self):
        '''
        Set in self.state.x the initial guess for the solution.

        Reimplement this to improve convergence.
        '''
        pass
        
    def checkConv(self):
        '''
        Check if the residue is small enough so we can accept the current iteration solution.

        Reimplement this if more advanced convergence checking is needed.
        '''
        if self.state.N==0:
            return (True,0)
        
        res = self.state.b
        norm = linalg.norm(res)/self.state.N
        return (norm<1e-8, norm)

    def dampStep(self, dx):
        '''
        Damp the solution update step.

        Reimplement this to enforce non-negative variable, perform damping, line search etc.
        '''
        return dx
    
    def solve(self):
        '''
        Solve the PDE.
        '''
        maxiter=50
        trace = True
        for iter in xrange(0,maxiter):
            self.state.clearFunJac()
            self.calcFunJac()
            self.state.assembleJac()

            flagConv, err = self.checkConv()
            if trace:
                print iter, err
            if flagConv:
                break

            dx = np.negative(dsolve.spsolve(self.state.J, self.state.b))
            self.state.dx = self.dampStep(dx)

            self.state.x = np.add(self.state.x, self.state.dx)

        #x = splinalg.bicg(eqns.J, eqns.b, None, 1e-6)
