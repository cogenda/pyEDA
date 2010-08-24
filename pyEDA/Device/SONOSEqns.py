'''
Created on Jul 5, 2009

@author: hash
'''

from pyEDA.PDE.AutoDeriv import *
from pyEDA.PDE.NLEqns import *
from pyEDA.Mesh.Mesh1D import *
from pyEDA.FVMEqn.FVMEqn import *
import pyEDA.Device.PhysUnit as Unit
import math
import numpy as np

class OxideMaterial(object):
    def __init__(self):
        self.affinity = 0.97*Unit.V
        self.eps = 3.9*Unit.eps0
        self.Eg = 9.0*Unit.V
        self.me_tnl = 0.5*Unit.me
        self.mh_tnl = 0.5*Unit.me
        
class TrappingMaterial(object):
    def __init__(self):
        self.VT = 300*Unit.kb*Unit.K
        self.eps = 7.5*Unit.eps0
        self.affinity = 1.4*Unit.V
        self.Eg = 5.5*Unit.V
        self.me_tnl = 0.5*Unit.me
        self.mh_tnl = 0.5*Unit.me

        self.vsat = 1e7*Unit.cm/Unit.s
        self.Esat = 1e6*Unit.V/Unit.cm
        self.D = 1.0*Unit.cm*Unit.cm/Unit.s

        self.Snp = 8.6e-13*Unit.cm*Unit.cm
        self.Spn = 2.7e-13*Unit.cm*Unit.cm
        self.Sno = 0.1*self.Snp
        self.Spo = 0.1*self.Spn

class SubstrateMaterial(object):
    def __init__(self, doping):
        self.VT = 0.0258*Unit.V
        self.C = doping
        self.eps = 11.7*Unit.eps0
        self.affinity = 4.17*Unit.V
        self.Eg = 1.12*Unit.V
        self.ni = 1.45e10*pow(Unit.cm, -3.0)
        self.me_dos = 0.5*Unit.me
        
        self.n0=0.5*(self.C+math.sqrt(self.C*self.C+4*self.ni*self.ni))
        self.p0=0.5*(-self.C+math.sqrt(self.C*self.C+4*self.ni*self.ni))
        if self.C>0:
            self.phib = + self.VT * math.log(self.n0/self.ni) # n-type
        else:  
            self.phib = - self.VT * math.log(self.p0/self.ni) # p-type
        self.phib = self.phib - 0.5*self.Eg - self.affinity

class GateBoundaryEqn(BoundaryEqn):
    def __init__(self):
        self.V = 0
        self.workfunc = 4.17*Unit.V
    
    def setVoltage(self, V):
        self.V = V
        
    def voltage(self):
        return self.V

    def cellEqn(self, state, cell):
        phi = state.getVar(cell.vars[0])

        state.resetEqn(cell.vars[0])
        state.setFunJac(cell.vars[0], phi-(self.V-self.workfunc))

class MOSSubstrateEqn(BoundaryEqn):
    def __init__(self):
        self.material = None

    def setMaterial(self, material):
        self.material = material
    
    def cellEqn(self, state, cell):
        material = self.material
        phi = state.getVar(cell.vars[0]) - material.phib
        f1 = material.p0 * ( exp(-phi/material.VT) + phi/material.VT - 1)
        f2 = material.n0 * ( exp(phi/material.VT) - phi/material.VT - 1)
        F = material.VT * Unit.e / material.eps * ( f1 + f2 )
        if phi>0:
            E = pow(F,0.5)
        else:
            E = -pow(F,0.5)
        
        state.setFunJac(cell.vars[0], E*material.eps)

class InsulatorRegionEqn(RegionEqn):
    def __init__(self):
        super(InsulatorRegionEqn, self).__init__()
        self.material = None
    
    def eqnPerCell(self):
        return 1

    def setMaterial(self, material):
        self.material = material

    def elemEqn(self, state, elem):
        material = self.material
        varA = elem.cells[0].vars[0]
        varB = elem.cells[1].vars[0]
        h = elem.len
        
        Va = state.getVar(varA)
        Vb = state.getVar(varB)
        state.setFunJac(varA, (Va-Vb)/h*material.eps)
        state.setFunJac(varB, (Vb-Va)/h*material.eps)
        
    def damp(self, state, cell, dx):
        vars = cell.vars
        V = state.getVars(vars)
        dV = dx[vars[0]]
        if dV>0.1:
            dV = 0.1*(1+math.log(dV/0.1))
        elif dV<-0.1:
            dV = -0.1*(1+math.log(-dV/0.1))
        dx[vars[0]] = dV

    def initGuess(self, state, cell):
        state.setVar(cell.vars[0], -5.0*Unit.V)

class TrappingRegionEqn(RegionEqn):
    def __init__(self):
        super(TrappingRegionEqn, self).__init__()
        self.material = None
        
    def setMaterial(self, material):
        self.material = material

    def velocity(self, E):
        tmp = E/self.material.Esat
        return self.material.vsat*(tmp/(1.0+abs(tmp)))

    def eqnPerCell(self):
        return 5
    
    def cellEqn(self, state, cell):
        material = self.material
        vars = cell.vars
        V, n, p, nT, pT = state.getVars(vars)
        dVdt, dndt, dpdt, dnTdt, dpTdt = state.getTimeDerivs(vars)
        NT = cell.fields['NT']
        oT = NT-nT-pT
        
        rho = Unit.e*(p-n+pT-nT)
        Cpo = material.Spo * material.vsat * p * oT  # hole capture at neutral trap
        Cpn = material.Spn * material.vsat * p * nT  # hole capture at negative trap
        Cno = material.Sno * material.vsat * n * oT  # electron capture at neutral trap
        Cnp = material.Snp * material.vsat * n * pT  # electron capture at positive trap
        
        state.setFunJac(vars[0], - rho*cell.volume())
        state.setFunJac(vars[1], -(Cnp+Cno+dndt)*cell.volume())
        state.setFunJac(vars[2], -(Cpn+Cpo+dpdt)*cell.volume())
        state.setFunJac(vars[3], (Cno-Cpn-dnTdt)*cell.volume())
        state.setFunJac(vars[4], (Cpo-Cnp-dpTdt)*cell.volume())
    
    def elemEqn(self, state, elem):
        material = self.material
        varsA = elem.cells[0].vars
        varsB = elem.cells[1].vars
        h = elem.len
        
        Va, na, pa, nTa, pTa = state.getVars(varsA)
        Vb, nb, pb, nTb, pTb = state.getVars(varsB)
        
        Fn = 0.5*(na+nb) * self.velocity((Vb-Va)/h) + material.D * (na-nb)/h
        Fp = 0.5*(pa+pb) * self.velocity((Va-Vb)/h) + material.D * (pa-pb)/h
        
        state.setFunJac(varsA[0], (Va-Vb)/h*material.eps)
        state.setFunJac(varsA[1], -Fn)
        state.setFunJac(varsA[2], -Fp)

        state.setFunJac(varsB[0], (Vb-Va)/h*material.eps)
        state.setFunJac(varsB[1], Fn)
        state.setFunJac(varsB[2], Fp)

    def initGuess(self, state, cell):
        NT = cell.fields['NT']
        state.setVar(cell.vars[0], -5.0*Unit.V)
        state.setVar(cell.vars[1], 1e-10*NT)
        state.setVar(cell.vars[2], 1e-10*NT)
        state.setVar(cell.vars[3], 1e-10*NT)
        state.setVar(cell.vars[4], 1e-10*NT)
        
    def damp(self, state, cell, dx):
        vars = cell.vars
        V, n, p, nT, pT = state.getVars(vars)
        dV, dn, dp, dnT, dpT = dx[vars[0]], dx[vars[1]], dx[vars[2]], dx[vars[3]], dx[vars[4]] 
        
        if dV>0.1:
            dV = 0.1*(1+math.log(dV/0.1))
        elif dV<-0.1:
            dV = -0.1*(1+math.log(-dV/0.1))
        if n+dn<0:
            dn = -0.99*n.getVal()
        if p+dp<0:
            dp = -0.99*p.getVal()
        if nT+dnT<0:
            dnT = -0.99*nT.getVal()
        if pT+dpT<0:
            dpT = -0.99*pT.getVal()
        dx[vars[0]], dx[vars[1]], dx[vars[2]], dx[vars[3]], dx[vars[4]] = (dV, dn, dp, dnT, dpT)


class SimpleIFEqn(InterfaceEqn):
    def __init__(self):
        super(SimpleIFEqn, self).__init__()
    
    def cellPairEqn(self, state, cell1, cell2):
        state.connectVar(cell1.vars[0], cell2.vars[0])


