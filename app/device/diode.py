'''
Created on Jun 18, 2009

@author: hash
'''

from pyEDA.PDE.AutoDeriv import *
from pyEDA.PDE.NLEqns import *
from pyEDA.Mesh.Mesh1D import *
from pyEDA.FVMEqn.FVMEqn import *
import pyEDA.Device.PhysUnit as Unit

import numpy as np
import scipy
import math

class Silicon(object):
    def __init__(self):
        self.ni = 1.45e10*pow(Unit.cm,-3)
        self.mun = 1000*Unit.cm*Unit.cm/Unit.V/Unit.s
        self.mup = 400*Unit.cm*Unit.cm/Unit.V/Unit.s
        self.tau = 1e-7*Unit.s
        self.eps = 11.7*Unit.eps0

class DDRegionEqn(RegionEqn):
    def __init__(self):
        super(DDRegionEqn,self).__init__()
        
    def eqnPerCell(self):
        return 3
    
    def cellEqn(self, state, cell):
        mtl = cell.region.material
        ni,mun,mup,tau,eps = (mtl.ni, mtl.mun, mtl.mup, mtl.tau, mtl.eps)
        vol = cell.volume()
        vars = cell.vars
            
        dummy, n, p = state.getVars(vars)
        rho = Unit.e*(p-n+cell.fields['C'])
        Rsrh = (n*p-ni*ni)/(tau*(n+ni)+tau*(p+ni))
            
        state.setFunJac(vars[0], -rho*vol)
        state.setFunJac(vars[1], -Rsrh*vol)
        state.setFunJac(vars[2], -Rsrh*vol)
    
    def elemEqn(self, state, elem):
        mtl = elem.region.material
        ni,mun,mup,tau,eps = (mtl.ni, mtl.mun, mtl.mup, mtl.tau, mtl.eps)
        VT = 0.0258*Unit.V
        
        varsA = elem.cells[0].vars
        varsB = elem.cells[1].vars
        
        h = elem.len
        Va, na, pa = state.getVars(varsA)
        Vb, nb, pb = state.getVars(varsB)
        
        nmid = na*aux2((Va-Vb)/2.0/VT) + nb*aux2((Vb-Va)/2.0/VT)
        dndx = aux1((Va-Vb)/2.0/VT) * (nb-na)/h
        pmid = pa*aux2((Vb-Va)/2.0/VT) + pb*aux2((Va-Vb)/2.0/VT)
        dpdx = aux1((Va-Vb)/2.0/VT) * (pb-pa)/h
    
        Fn =  mun * (-nmid* (Va-Vb)/h - VT*dndx)
        Fp =  mup * ( pmid* (Va-Vb)/h - VT*dpdx)
    
        state.setFunJac(varsA[0], (Va-Vb)/h*eps)
        state.setFunJac(varsA[1], -Fn)
        state.setFunJac(varsA[2], -Fp)
    
        state.setFunJac(varsB[0], (Vb-Va)/h*eps)
        state.setFunJac(varsB[1], Fn)
        state.setFunJac(varsB[2], Fp)

    def initGuess(self, state, cell):
        ni = cell.region.material.ni
        VT = 0.0258*Unit.V

        C = cell.fields['C']
        n=0.5*(C+math.sqrt(C*C+4*ni*ni))
        p=0.5*(-C+math.sqrt(C*C+4*ni*ni))
        phi=VT*math.log(n/ni)

        state.setVar(cell.vars[0],phi)
        state.setVar(cell.vars[1], n)
        state.setVar(cell.vars[2], p)        

class OhmicBoundaryEqn(BoundaryEqn):
    def __init__(self):
        super(OhmicBoundaryEqn,self).__init__()
        self.voltage = 0.0
        
    def setVoltage(self, V):
        self.voltage = V

    def cellEqn(self, state, cell, dt=0):
        VT = 0.0258*Unit.V
        ni = cell.region.material.ni
        C = cell.fields['C']
        V, n, p = state.getVars(cell.vars)
        
        ni = cell.region.material.ni
        nb=0.5*(C+math.sqrt(C*C+4*ni*ni))
        pb=0.5*(-C+math.sqrt(C*C+4*ni*ni))
        Vb=self.voltage + VT*math.log(nb/ni)

        state.resetEqn(cell.vars[0])
        state.resetEqn(cell.vars[1])
        state.resetEqn(cell.vars[2])
        state.setFunJac(cell.vars[0], V-Vb)
        state.setFunJac(cell.vars[1], n-nb)
        state.setFunJac(cell.vars[2], p-pb)

class Diode(Mesh1D):
    def __init__(self):
        NN = 200

        xx = np.linspace(-1e-4*Unit.cm, 1e-4*Unit.cm, NN+1)
        Mesh1D.__init__(self, xx, 
                        rgns=[(0, NN, 'silicon')],
                        bnds=[(0,'anode'),(NN,'cathode')])

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
        self.setRegionEqn('silicon', DDRegionEqn())
        self.device.regions[0].material = Silicon()
        self.bcAnode = OhmicBoundaryEqn()
        self.bcCathode = OhmicBoundaryEqn()
        self.setBoundaryEqn('anode', self.bcAnode)
        self.setBoundaryEqn('cathode', self.bcCathode)
        self.setupEqns()

if __name__ == '__main__':
    diode = Diode()
    diodeEqns = DiodeEqns(diode)
    diodeEqns.initGuess()
    diodeEqns.bcAnode.setVoltage(0.0)
    diodeEqns.solve()

    diodeEqns.bcAnode.setVoltage(0.3)
    diodeEqns.solve()
    
