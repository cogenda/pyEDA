'''
Created on Jun 13, 2009

@author: hash
'''
from pyEDA.PDE.AutoDeriv import *
from pyEDA.PDE.NLEqns import *
from pyEDA.Mesh.Mesh1D import *
from pyEDA.FVMEqn.FVMEqn import *
from pyEDA.Device.DDEqns import *

import numpy as np
import scipy

import pyEDA.Device.PhysUnit as Unit

class SiO2(InsulatorMaterial):
    def __init__(self):
        super(SiO2, self).__init__()

class Silicon(SemiconductorMaterial):
    def __init__(self):
        super(Silicon, self).__init__()

class MOS(Mesh1D):
    def __init__(self):
        NOX = 4
        NSi = 100
        Nd = -1e17*pow(Unit.cm,-3.0)

        xx1 = np.linspace(-10e-6*Unit.cm, 0.0, NOX+1)
        xx2 = np.linspace(0.0, 1.0e-4*Unit.cm, NSi+1)
        xx = np.unique(np.concatenate((xx1,xx2)))
        Mesh1D.__init__(self, xx, 
                        rgns=[(0, NOX, 'Gox'), (NOX,NOX+NSi, 'Substrate')],
                        bnds=[(0,'Gate'),(NOX+NSi,'Substrate')])
        self.setFieldByFunc(1, 'C', lambda x : Nd)
        
        self.setRegionMaterial('Gox', SiO2())
        self.setRegionMaterial('Substrate', Silicon())

class MOSEqn (FVMEqns):
    def __init__(self, device):
        super(MOSEqn, self).__init__(device)


        self.OxideEqn  = InsulatorRegionEqn()

        self.SiliconEqn  = SemiconductorRegionEquEqn()
    
    	self.setRegionEqn('Gox', self.OxideEqn)
        self.setRegionEqn('Substrate', self.SiliconEqn)

        self.bcGate = GateBoundaryEqn()
        self.bcSubstrate = OhmicBoundaryEquEqn()
        self.setBoundaryEqn('Gate', self.bcGate)
        self.setBoundaryEqn('Substrate', self.bcSubstrate)

        self.setInterfaceEqn('Gox', 'Substrate', SimpleIFEqn())

        self.setupEqns()

    def setVg(self, Vg):
        self.bcGate.setVoltage(Vg)

    def setVb(self, Vb):
        self.bcSubstrate.setVoltage(Vb)
        
if __name__=='__main__':
    mos = MOS()
    mosEqn = MOSEqn(mos)
    mosEqn.initGuess()

    mosEqn.setVg(10.0)
    mosEqn.solve()

    vi = mos.getVarIdx('Substrate', 0) # potential in Substrate
    print '-----------------------'
    print 'Potential in Substrate:'
    print mosEqn.state.x[vi]

