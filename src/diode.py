'''
Created on Jun 15, 2009

@author: hash
'''

from PDE.AutoDeriv import *
from PDE.NLEqns import *
from Mesh.Mesh1D import *
import PhysUnit as Unit

import numpy as np
import scipy
import math


class Diode(Mesh1D):
    def __init__(self):
        NN = 200

        xx = np.linspace(-1e-4*Unit.cm, 1e-4*Unit.cm, NN+1)
        Mesh1D.__init__(self, xx, 
                        rgns=[(0, NN, 'silicon', 3)],
                        bnds=[(0,'anode'),(NN,'cathode')])

        def doping(x):
            if x<0:
                return -1e17*pow(Unit.cm,-3)
            elif x>0:
                return 1e17*pow(Unit.cm,-3)
            else:
                return 0.0
            
        self.setFieldByFunc(0, 'C', doping)

class Material:
    def __init__(self):
        self.ni = 1.45e10*pow(Unit.cm,-3)
        self.mun = 1000*Unit.cm*Unit.cm/Unit.V/Unit.s
        self.mup = 400*Unit.cm*Unit.cm/Unit.V/Unit.s
        self.tau = 1e-7*Unit.s
        self.eps = 11.7*Unit.eps0
        
class DiodeEqn(NLEqns):
    def __init__(self, device, material):
        NLEqns.__init__(self)
        self.device = device
        self.material = material
        self.state = NLEqnState(device.eqnCnt)

    def setVd(self, Vd):
        self.device.boundaries[0].voltage = Vd
        self.device.boundaries[1].voltage = 0.0
        
    def calcOhmic(self, b):
        VT = 0.0258*Unit.V
        cell,dummy = b.cells[0]

        C = cell.fields['C']
        V = b.voltage
        
        ni = self.material.ni
        n=0.5*(C+math.sqrt(C*C+4*ni*ni))
        p=0.5*(-C+math.sqrt(C*C+4*ni*ni))
        phi=V + VT*math.log(n/ni)
        return [phi,n,p]
            
    def initGuess(self):
        ni = self.material.ni
        VT = 0.0258*Unit.V
        for c in self.device.regions[0].cells:
            C = c.fields['C']
            n=0.5*(C+math.sqrt(C*C+4*ni*ni))
            p=0.5*(-C+math.sqrt(C*C+4*ni*ni))
            phi=VT*math.log(n/ni)

            self.state.setVar(c.vars[0],phi)
            self.state.setVar(c.vars[1], n)
            self.state.setVar(c.vars[2], p)        
        
    def calcFunJac(self):
        VT = 0.0258*Unit.V
        mun=self.material.mun
        mup=self.material.mup
        ni = self.material.ni
        tau = self.material.tau
        epsSi = self.material.eps

        for elem in self.device.regions[0].elems:
            varsA = elem.cells[0].vars
            varsB = elem.cells[1].vars
            
            h = elem.len
            Va, na, pa = self.state.getVars(varsA)
            Vb, nb, pb = self.state.getVars(varsB)
            
            nmid = na*aux2((Va-Vb)/2.0/VT) + nb*aux2((Vb-Va)/2.0/VT)
            dndx = aux1((Va-Vb)/2.0/VT) * (nb-na)/h
            pmid = pa*aux2((Vb-Va)/2.0/VT) + pb*aux2((Va-Vb)/2.0/VT)
            dpdx = aux1((Va-Vb)/2.0/VT) * (pb-pa)/h

            Fn =  mun * (-nmid* (Va-Vb)/h - VT*dndx)
            Fp =  mup * ( pmid* (Va-Vb)/h - VT*dpdx)

            self.state.setFunJac(varsA[0], (Va-Vb)/h*epsSi)
            self.state.setFunJac(varsA[1], -Fn)
            self.state.setFunJac(varsA[2], -Fp)

            self.state.setFunJac(varsB[0], (Vb-Va)/h*epsSi)
            self.state.setFunJac(varsB[1], Fn)
            self.state.setFunJac(varsB[2], Fp)
                
        for cell in self.device.regions[0].cells:
            vol = cell.volume()
            vars = cell.vars
            
            dummy, n, p = self.state.getVars(vars)
            rho = Unit.e*(p-n+cell.fields['C'])
            Rsrh = (n*p-ni*ni)/(tau*(n+ni)+tau*(p+ni))
            
            self.state.setFunJac(vars[0], -rho*vol)
            self.state.setFunJac(vars[1], -Rsrh*vol)
            self.state.setFunJac(vars[2], -Rsrh*vol)
            
        for bnd in self.device.boundaries:
            cell,dummy = bnd.cells[0]
            vol = cell.volume()
            vars = cell.vars
            
            V, n, p = self.state.getVars(vars)
            Vbnd, nbnd, pbnd = self.calcOhmic(bnd)
            
            self.state.resetEqn(vars[0])
            self.state.resetEqn(vars[1])
            self.state.resetEqn(vars[2])
            self.state.setFunJac(vars[0], V-Vbnd)
            self.state.setFunJac(vars[1], n-nbnd)
            self.state.setFunJac(vars[2], p-pbnd)

if __name__ == '__main__':
    diode = Diode()
    diodeEqn = DiodeEqn(diode, Material())
    diodeEqn.initGuess()

    diodeEqn.setVd(0.1)
    diodeEqn.solve()

    vi = diode.getVarIdx(0, 0)
    ni = diode.getVarIdx(0, 1)
    pi = diode.getVarIdx(0, 2)
#    print diodeEqn.state.x[vi]
#    print diodeEqn.state.x[ni]
#    print diodeEqn.state.x[pi]
    