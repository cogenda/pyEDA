'''
Created on Jun 13, 2009

@author: hash
'''
from pyEDA.PDE.AutoDeriv import *
from pyEDA.PDE.NLEqns import *
from pyEDA.Mesh.Mesh1D import *

import numpy as np
import scipy

import PhysUnit as Unit

class MOS(Mesh1D):
    def __init__(self):
        NOX = 4
        NSi = 100
        Nd = -1e17*pow(Unit.cm,-3.0)

        xx1 = np.linspace(-10e-6*Unit.cm, 0.0, NOX+1)
        xx2 = np.linspace(0.0, 1.0e-4*Unit.cm, NSi+1)
        xx = np.unique(np.concatenate((xx1,xx2)))
        Mesh1D.__init__(self, xx, 
                        rgns=[(0, NOX, 'oxide', 1), (NOX,NOX+NSi, 'silicon', 1)],
                        bnds=[(0,'gate'),(NOX+NSi,'substrate')])
        self.setFieldByFunc(1, 'C', lambda x : Nd)
        
class MOSEqn (NLEqns):
    def __init__(self, device):
        NLEqns.__init__(self)
        self.device = device
        self.state = NLEqnState(device.eqnCnt)
        self.bcs = [0.0, 0.0]
    
    def setVg(self, Vg):
        self.device.boundaries[0].voltage = Vg
        self.device.boundaries[1].voltage = 0.0
        
    def calcFunJac(self):
        epsOx=Unit.eps0*3.9
        epsSi=Unit.eps0*11.7
        VT = Unit.V*0.0258

        for r in [0,1]:
            for elem in self.device.regions[r].elems:
                varA = elem.cells[0].vars[0]
                varB = elem.cells[1].vars[0]
                
                h = elem.len
                Va = self.state.getVar(varA)
                Vb = self.state.getVar(varB)
                
                if r==0:
                    eps=epsOx
                else:
                    eps=epsSi

                self.state.setFunJac(varA, (Va-Vb)/h*eps)
                self.state.setFunJac(varB, (Vb-Va)/h*eps)
                
            if r==1:
                for cell in self.device.regions[r].cells:
                    C = cell.fields['C']
                    V = self.state.getVar(cell.vars[0])
                    rho = Unit.e*C * (1-exp(-V/VT))
                    self.state.setFunJac(cell.vars[0], - rho*cell.volume())
            
        for interface in self.device.interfaces:
            (c1,c2),dummy = interface.cellPairs[0]
            self.state.connectVar(c1.vars[0], c2.vars[0])
        
        for bnd in self.device.boundaries:
            cell,dummy = bnd.cells[0]
            var = cell.vars[0]
            V = self.state.getVar(var)

            self.state.resetEqn(var)
            self.state.setFunJac(var, V-bnd.voltage)

if __name__=='__main__':
    mos = MOS()
    mosEqn = MOSEqn(mos)

    mosEqn.setVg(5)
    mosEqn.solve()

    #import Gnuplot
    #g = Gnuplot.Gnuplot()
    #g('set data style linespoints')
    
    vi = mos.getVarIdx(1, 0)
    
    print mosEqn.state.x[vi]

    #result = mos.extractVar(mosEqn.state.x, None, 0)
    #d = Gnuplot.Data(mos.xx, result, title='potential in silicon');
    #g.plot(d) 
