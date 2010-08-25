'''
Created on Jun 18, 2009

@author: hash
'''

from pyEDA.PDE.AutoDeriv import *
from pyEDA.PDE.NLEqns import *
from pyEDA.Mesh.Mesh1D import *
from pyEDA.FVMEqn.FVMEqn import *
import pyEDA.Device.PhysUnit as Unit
from pyEDA.Device.DDEqns import *

import numpy as np
import scipy
import math

class Silicon(SemiconductorMaterial):
    def __init__(self):
        super(Silicon, self).__init__()

class Diode(Mesh1D):
    def __init__(self):
        NN = 200

        xx = np.linspace(-1e-4*Unit.cm, 1e-4*Unit.cm, NN+1)
        Mesh1D.__init__(self, xx, 
                        rgns=[(0, NN, 'silicon')],
                        bnds=[(0,'anode'),(NN,'cathode')])

        self.setRegionMaterial('silicon', Silicon())

        def doping(x):
            if x<0:
                return -1e17*pow(Unit.cm,-3)
            elif x>0:
                return 1e17*pow(Unit.cm,-3)
            else:
                return 0.0
            
        self.setFieldByFunc(0, 'C', doping)

class DiodeEqns(FVMEqns):
    def __init__(self, device):
        super(DiodeEqns, self).__init__(device)

        self.eqnSi = SemiconductorRegionEqn()
        self.setRegionEqn('silicon', self.eqnSi)

        self.bcAnode = OhmicBoundaryEqn()
        self.bcCathode = OhmicBoundaryEqn()
        self.setBoundaryEqn('anode', self.bcAnode)
        self.setBoundaryEqn('cathode', self.bcCathode)
        self.setupEqns()

if __name__ == '__main__':
    diode = Diode()
    diodeEqn = DiodeEqns(diode)
    diodeEqn.initGuess()
    diodeEqn.bcAnode.setVoltage(0.0)
    diodeEqn.solve()

    diodeEqn.bcAnode.setVoltage(0.3)
    diodeEqn.solve()
    
    vi = diode.getVarIdx('silicon', 0) # potential in Substrate
    print '-----------------------'
    print 'Potential in silicon:'
    print diodeEqn.state.x[vi]

    cmc=pow(Unit.cm,-3)
    ni = diode.getVarIdx('silicon', 1) # electron conc in Substrate
    print '-----------------------'
    print 'Elec. conc. in silicon:'
    print diodeEqn.state.x[ni]/cmc

    pi = diode.getVarIdx('silicon', 2) # hole conc in Substrate
    print '-----------------------'
    print 'Hole conc. in silicon:'
    print diodeEqn.state.x[pi]/cmc

