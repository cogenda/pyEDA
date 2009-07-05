'''
Created on Jun 16, 2009

@author: hash
'''

from PDE.AutoDeriv import *
from PDE.NLEqns import *
from Mesh.Mesh1D import *
import PhysUnit as Unit

import numpy as np
import scipy
import math

class RegionEqn:
    def eqnPerCell(self):
        return 0
    
    def cellEqn(self, state, cell, dt=0):
        pass
    
    def elemEqn(self, state, elem, dt=0):
        pass
    
    def initGuess(self, state, cell, dt=0):
        pass
    
    def saveTimeStep(self, dt):
        pass

class BoundaryEqn:
    def cellEqn(self, state, cell, dt=0):
        pass

    def saveTimeStep(self, dt):
        pass
    
class InterfaceEqn:
    def cellPairEqn(self, state, cell1, cell2, dt=0):
        pass

    def saveTimeStep(self, dt):
        pass

class SubstrateBoundaryEqn(BoundaryEqn):
    def __init__(self):
        self.V = 0
        self.VT = 0.0258*Unit.V
        self.affinity = 4.17*Unit.V
        self.Eg = 1.12*Unit.V
        self.ni = 1.45e10*pow(Unit.cm,-3.0)
    
    def setVoltage(self, V):
        self.V = V

    def cellEqn(self, state, cell, dt=0):
        C = cell.fields['C']
        phi = state.getVar(cell.vars[0])

        n=0.5*(C+math.sqrt(C*C+4*self.ni*self.ni))
        phiB = self.V + self.VT*math.log(n/self.ni) - self.affinity - self.Eg/2.0 
        
        state.resetEqn(cell.vars[0])
        state.setFunJac(cell.vars[0], phi-phiB)
        
class GateBoundaryEqn(BoundaryEqn):
    def __init__(self):
        self.V = 0
        self.workfunc = 4.17*Unit.V
    
    def setVoltage(self, V):
        self.V = V

    def cellEqn(self, state, cell, dt=0):
        phi = state.getVar(cell.vars[0])

        state.resetEqn(cell.vars[0])
        state.setFunJac(cell.vars[0], phi-(self.V-self.workfunc))
        
        
class InsulatorRegionEqn(RegionEqn):
    def __init__(self):
        self.eps = 3.9*Unit.eps0
        self.affinity = 1.0*Unit.V
        self.Eg = 9.0*Unit.V

    def eqnPerCell(self):
        return 1

    def elemEqn(self, state, elem, dt=0):
        varA = elem.cells[0].vars[0]
        varB = elem.cells[1].vars[0]
        h = elem.len
        
        Va = state.getVar(varA)
        Vb = state.getVar(varB)
        state.setFunJac(varA, (Va-Vb)/h*self.eps)
        state.setFunJac(varB, (Vb-Va)/h*self.eps)
        
    def initGuess(self, state, cell, dt=0):
        state.setVar(cell.vars[0], -5.0*Unit.V)
    
class SiliconRegionEqn(RegionEqn):
    def __init__(self):
        self.VT = 0.0258*Unit.V
        self.eps = 11.7*Unit.eps0
        self.affinity = 4.17*Unit.V
        self.Eg = 1.12*Unit.V
        self.ni = 1.45e10*pow(Unit.cm,-3.0)

    def eqnPerCell(self):
        return 1
    
    def cellEqn(self, state, cell, dt=0):
        C = cell.fields['C']
        Ei = state.getVar(cell.vars[0]) + self.affinity + self.Eg/2.0
        
        n = self.ni * exp(Ei/self.VT)
        p = self.ni * exp(-Ei/self.VT)
        rho = Unit.e* (C+p-n)
        state.setFunJac(cell.vars[0], - rho*cell.volume())
    
    def elemEqn(self, state, elem, dt=0):
        varA = elem.cells[0].vars[0]
        varB = elem.cells[1].vars[0]
        h = elem.len
        
        Va = state.getVar(varA)
        Vb = state.getVar(varB)
        state.setFunJac(varA, (Va-Vb)/h*self.eps)
        state.setFunJac(varB, (Vb-Va)/h*self.eps)

    def initGuess(self, state, cell):
        C = cell.fields['C']
        n=0.5*(C+math.sqrt(C*C+4*self.ni*self.ni))
        phi = self.VT*math.log(n/self.ni) - self.affinity - self.Eg/2.0 
        state.setVar(cell.vars[0], phi)

class TrappingRegionEqn(RegionEqn):
    def __init__(self):
        self.VT = 0.0258*Unit.V
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

    def velocity(self, E):
        tmp = E/self.Esat
        return self.vsat*(tmp/(1.0+abs(tmp)))

    def eqnPerCell(self):
        return 5
    
    def cellEqn(self, state, cell, dt=0):
        vars = cell.vars
        V, n, p, nT, pT = state.getVars(vars)
        NT = cell.fields['NT']
        oT = NT-nT-pT 
        
        rho = Unit.e*(p-n+pT-nT)
        Cpo = self.Spo * self.vsat * p * oT  # hole capture at neutral trap
        Cpn = self.Spn * self.vsat * p * nT  # hole capture at negative trap
        Cno = self.Sno * self.vsat * n * oT  # electron capture at neutral trap
        Cnp = self.Snp * self.vsat * n * pT  # electron capture at positive trap
        
        state.setFunJac(vars[0], - rho*cell.volume())
        state.setFunJac(vars[1], -(Cnp+Cno)*cell.volume())
        state.setFunJac(vars[2], -(Cpn+Cpo)*cell.volume())
        state.setFunJac(vars[3], (Cnp+Cno)*cell.volume())
        state.setFunJac(vars[4], (Cpn+Cpo)*cell.volume())
    
    def elemEqn(self, state, elem, dt=0):
        varsA = elem.cells[0].vars
        varsB = elem.cells[1].vars
        h = elem.len
        
        Va, na, pa, nTa, pTa = state.getVars(varsA)
        Vb, nb, pb, nTb, pTb = state.getVars(varsB)
        
        Fn = 0.5*(na+nb) * self.velocity((Vb-Va)/h) + self.D * (na-nb)/h
        Fp = 0.5*(pa+pb) * self.velocity((Va-Vb)/h) + self.D * (pa-pb)/h
        
        state.setFunJac(varsA[0], (Va-Vb)/h*self.eps)
        state.setFunJac(varsA[1], -Fn)
        state.setFunJac(varsA[2], -Fp)

        state.setFunJac(varsB[0], (Vb-Va)/h*self.eps)
        state.setFunJac(varsB[1], Fn)
        state.setFunJac(varsB[2], Fp)

    def initGuess(self, state, cell):
        NT = cell.fields['NT']
        state.setVar(cell.vars[0], -5.0*Unit.V)
        state.setVar(cell.vars[1], 1e-2*NT)
        state.setVar(cell.vars[2], 1e-2*NT)
        state.setVar(cell.vars[3], 1.1e-2*NT)
        state.setVar(cell.vars[4], 0.9e-2*NT)

class OxideSiliconIFEqn(InterfaceEqn):
    def __init__(self):
      pass
    
    def cellPairEqn(self, state, cell1, cell2):
        state.connectVar(cell1.vars[0], cell2.vars[0])

class OxideTrappingIFEqn(InterfaceEqn):
    def __init__(self):
        self.injection = 0.0
        self.J = 1e-3*Unit.A*pow(Unit.cm,-2.0)
        
    def cellPairEqn(self, state, cell1, cell2):
        state.connectVar(cell1.vars[0], cell2.vars[0])

        if cell1.region.name=='trapping':
            cell = cell1
        elif cell2.region.name=='trapping':
            cell = cell2
        
        vars = cell.vars
        state.setFunJac(vars[1], ADVar(self.injection*self.J/Unit.e))

class SONOS(Mesh1D):
    def __init__(self):
        NOX = 20
        NSi = 40
        Nd = -1e17*pow(Unit.cm,-3.0)
        NT = 1e18*pow(Unit.cm,-3.0)

        xx1 = np.linspace(-20e-7*Unit.cm, 0.0, NOX+1)
        xx2 = np.linspace(0.0, 2.0e-5*Unit.cm, NSi+1)
        xx = np.unique(np.concatenate((xx1,xx2)))
        Mesh1D.__init__(self, xx, 
                        rgns=[(0, 10, 'blocking'),
                              (10,15, 'trapping'),
                              (15,NOX, 'oxide'), 
                              (NOX,NOX+NSi, 'silicon')],
                        bnds=[(0,'gate'),(NOX+NSi,'substrate')])
        self.setFieldByFunc(1, 'NT', lambda x : NT)
        self.setFieldByFunc(3, 'C', lambda x : Nd)

class SONOSEqns(NLEqns):
    def __init__(self, device):
        NLEqns.__init__(self)

        self.device = device
        self.regionEqns=[]
        self.interfaceEqns=[]
        self.boundaryEqns=[]

        eqnCnt=0
        for region in device.regions:
            if region.name=='oxide'  or region.name=='blocking':
                eqn =InsulatorRegionEqn()
            elif region.name=='silicon':
                eqn = SiliconRegionEqn()
            elif region.name=='trapping':
                eqn = TrappingRegionEqn()
            else:
                raise ValueError

            for cell in region.cells:
                eqnPerCell = eqn.eqnPerCell()
                cell.vars = xrange(eqnCnt, eqnCnt+eqnPerCell)
                eqnCnt += eqnPerCell

            self.regionEqns.append(eqn)

        for boundary in device.boundaries:
            if boundary.name=='gate':
                self.boundaryEqns.append( GateBoundaryEqn() )
            elif boundary.name=='substrate':
                self.boundaryEqns.append( SubstrateBoundaryEqn() )
            else:
                raise ValueError
                
        for i,interface in enumerate(device.interfaces):
            if matchIFName(interface.name, 'oxide', 'silicon'):
                self.interfaceEqns.append( OxideSiliconIFEqn() )
            elif matchIFName(interface.name, 'oxide', 'trapping'):
                self.ifc1 = OxideTrappingIFEqn()
                self.ifc1.injection=1.0
                self.interfaceEqns.append( self.ifc1 )
            elif matchIFName(interface.name, 'blocking', 'trapping'):
                self.ifc2 = OxideTrappingIFEqn()
                self.ifc2.injection=-1.0
                self.interfaceEqns.append( self.ifc2 )
            else:
                raise ValueError
        
        self.eqnCnt = eqnCnt
        self.state = NLEqnState(self.eqnCnt)

    def initGuess(self):
        for r,region in enumerate(self.device.regions):
            initGuess = self.regionEqns[r].initGuess
            for cell in region.cells:
                initGuess(self.state, cell)
        
    def calcFunJac(self):
        for r,region in enumerate(self.device.regions):
            elemEqn = self.regionEqns[r].elemEqn
            cellEqn = self.regionEqns[r].cellEqn
            for elem in region.elems:
                elemEqn(self.state, elem)
            
            for cell in region.cells:
                cellEqn(self.state, cell)
        
        for i,interface in enumerate(self.device.interfaces):
            (c1,c2),dummy = interface.cellPairs[0]
            self.interfaceEqns[i].cellPairEqn(self.state, c1, c2)
        
        for b,boundary in enumerate(self.device.boundaries):
            cell,dummy = boundary.cells[0]
            self.boundaryEqns[b].cellEqn(self.state, cell)
            
        print self.state.J
            
def matchIFName(ifname, r1name, r2name):
    if ifname==r1name+'|'+r2name or ifname==r2name+'|'+r1name:
        return True
    else:
        return False


if __name__ == '__main__':
    sonos = SONOS()
    sonosEqns = SONOSEqns(sonos)
    sonosEqns.initGuess()
    
    sonosEqns.boundaryEqns[0].setVoltage(0)
    sonosEqns.solve()
#    sonosEqns.boundaryEqns[0].setVoltage(5)
#    sonosEqns.solve()
#    sonosEqns.boundaryEqns[0].setVoltage(10)
#    sonosEqns.solve()
#    sonosEqns.ifc1.injection=1.0
#    sonosEqns.ifc2.injection=-1.0
#    sonosEqns.boundaryEqns[0].setVoltage(10)
#    sonosEqns.solve()
    
    vi = sonos.getVarIdx(0, 0)    
    print sonosEqns.state.x[vi]
    vi = sonos.getVarIdx(1, 0)    
    print sonosEqns.state.x[vi]
    vi = sonos.getVarIdx(2, 0)    
    print sonosEqns.state.x[vi]
    vi = sonos.getVarIdx(3, 0)    
    print sonosEqns.state.x[vi]
