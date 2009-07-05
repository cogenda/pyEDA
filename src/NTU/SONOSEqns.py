'''
Created on Jul 5, 2009

@author: hash
'''

from PDE.AutoDeriv import *
from PDE.NLEqns import *
from Mesh.Mesh1D import *
from FVMEqn.FVMEqn import *
import PhysUnit as Unit

class TrappingMaterial(object):
    def __init__(self):
        self.VT = 300*Unit.kb*Unit.K
        self.eps = 4.5*Unit.eps0
        self.affinity = 2.3*Unit.V
        self.Eg = 6*Unit.V

        self.vsat = 1e7*Unit.cm/Unit.s
        self.Esat = 1e6*Unit.V/Unit.cm
        self.D = 1.0*Unit.cm*Unit.cm/Unit.s

        self.Snp = 8.6e-13*Unit.cm*Unit.cm
        self.Spn = 2.7e-13*Unit.cm*Unit.cm
        self.Sno = 0.1*self.Snp
        self.Spo = 0.1*self.Spn

class TrappingRegionEqn(RegionEqn):
    def __init__(self):
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
        
        if n+dn<0:
            dn = -0.99*n.getVal()
        if p+dp<0:
            dp = -0.99*p.getVal()
        if nT+dnT<0:
            dnT = -0.99*nT.getVal()
        if pT+dpT<0:
            dpT = -0.99*pT.getVal()
        dx[vars[0]], dx[vars[1]], dx[vars[2]], dx[vars[3]], dx[vars[4]] = (dV, dn, dp, dnT, dpT)
