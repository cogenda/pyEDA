'''
Created on Jun 18, 2009

@author: hash
'''

from pyEDA.PDE.AutoDeriv import *
from pyEDA.PDE.NLEqns import *
from pyEDA.Mesh.Mesh1D import *
from pyEDA.FVMEqn.FVMEqn import *
import pyEDA.Device.PhysUnit as Unit
from pyEDA.Device.SONOSEqns import *

import numpy as np
import scipy
import math

class ToxMaterial(OxideMaterial):
    def __init__(self):
        super(ToxMaterial, self).__init__()
        self.eps = 4.2*Unit.eps0
        self.Eg = 8.1*Unit.V

class BlockMaterial(OxideMaterial):
    def __init__(self):
        super(BlockMaterial, self).__init__()

class Si3N4Material(TrappingMaterial):
    def __init__(self):
        super(Si3N4Material, self).__init__()

class DrainIFEqn(InterfaceEqn):
    def __init__(self, elec=False, hole=False):
        super(DrainIFEqn, self).__init__()
        self.material=None
        self.elecDrain = elec
        self.holeDrain = hole
        
    def setMaterial(self, material):
        self.material = material

    def setDrainMode(self, elec=False, hole=False):
        self.elecDrain=elec
        self.holeDrain=hole
        
    def cellPairEqn(self, state, cell1, cell2):
        state.connectVar(cell1.vars[0], cell2.vars[0])
        
        cell=None
        if (cell1.region.name=='Si3N4'):
            cell=cell1
        else:
            cell=cell2
        material = self.material
        
        vars = cell.vars
        V, n, p, nT, pT = state.getVars(vars)
        
        if self.elecDrain:
            Fn = -material.vsat * n
            state.setFunJac(vars[1], Fn)

        if self.holeDrain:
            Fp = -material.vsat * p
            state.setFunJac(vars[2], Fp)



class NonLocalTunneling(object):
    def __init__(self, mesh):
        self.materials = {}
        self.injectMaterial = None
        self.mesh = mesh
        self.inject = None
        self.dir = None
        self.band = 'Ec' # or 'Ev'
        
    def setInjectBoundary(self, boundary_or_interface, material, dir):
        if isinstance(boundary_or_interface, FVM.Boundary) or isinstance(boundary_or_interface, FVM.Interface):
            self.inject = boundary_or_interface
            self.injectMaterial = material
            self.dir = dir
        else:
            raise TypeError
    
    def setMaterials(self, materials):
        self.materials = materials
    
    def __call__(self, state):
        iter = ElemIterator1D(self.mesh)

        polarity = 1.0
        if self.band=='Ev':
            polarity = -1.0

        if isinstance(self.injectMaterial, SubstrateMaterial):
            c,dummy = self.inject.cells[0]
            Einj0 = float(state.getVar(c.vars[0])) + self.injectMaterial.affinity
            if self.band == 'Ev':
                Einj0 += self.injectMaterial.Eg
            Ef = 0.0
        elif isinstance(self.injectMaterial, GateMaterial):
            c,dummy = self.inject.cells[0]
            Ef = float(state.getVar(c.vars[0])) + self.injectMaterial.workfunc
            Einj0 = Ef+0.1*polarity
        
            
        Einj = Einj0
        Erange = np.linspace(Einj0, Einj0-0.2*polarity, 11)
        dE = abs(Erange[1]-Erange[0])
        VT = 300*Unit.kb*Unit.K/Unit.e
        Jtot = 0
        for Einj in Erange:
            TE=0
            termElem = None
            termPos = None
            
            iter.setElem(self.inject, self.dir)
            for elem in iter:
                if self.dir>0:
                    c1,c2 = elem.cells
                else:
                    c2,c1 = elem.cells
                material = self.materials[elem.region.name]
                if self.band=='Ec':
                    me = material.me_tnl
                else:
                    me = material.mh_tnl

                prefac = 2.0/3.0 * math.sqrt(2*me*Unit.e)/Unit.hbar
                E1 = state.getVar(c1.vars[0]) + material.affinity
                E2 = state.getVar(c2.vars[0]) + material.affinity
                if self.band=='Ev':
                    E1=E1+material.Eg
                    E2=E2+material.Eg
                    
                phi1 = (Einj-E1)*polarity
                phi2 = (Einj-E2)*polarity
                if phi1<=0:
                    termElem = elem
                    termPos = 0
                    break
                if phi2<=0:
                    dx = phi1/(phi1-phi2)*abs(c1.node.pos - c2.node.pos)
                    TE += prefac * pow(phi1, 0.5) * dx
                    termElem = elem
                    termPos = phi1/(phi1-phi2)
                    break
                else:
                    dx = abs(c1.node.pos-c2.node.pos)
                    if (abs(E1-E2)<1e-8):
                        TE += prefac * pow(phi1, 0.5) * dx
                    else:
                        TE += prefac * (pow(phi1,1.5)-pow(phi2,1.5))/(phi1-phi2) * dx

            if not termElem==None and TE<30: 
                if not termElem.region.name =='Si3N4':
                    for elem in iter:
                        if elem.region.name =='Si3N4':
                            termElem = elem
                            termPos=0.0
                            break
                    
                if self.dir>0:
                    c1,c2 = elem.cells
                else:
                    c2,c1 = elem.cells
                dJ = exp(-2.0*TE) * math.log(1+ exp((Einj-Ef)/VT)) * dE
            
                prefac = self.injectMaterial.me_dos * Unit.kb * 300*Unit.K / 2.0 / (math.pi*math.pi) / pow(Unit.hbar,3.0)
                dJ = dJ * prefac
                
                if self.band=='Ec':
                    eqn1, eqn2 = c1.vars[1], c2.vars[1]
                else:
                    eqn1, eqn2 = c1.vars[2], c2.vars[2]
                    
                state.setFunJac(eqn1, dJ*(1-termPos))
                state.setFunJac(eqn2, dJ*(1-termPos))
                
                Jtot+=dJ
        print "......Jtot = %8e (A/cm^2)" % (Jtot/Unit.A * Unit.cm *Unit.cm)

class Trapping(Mesh1D):
    def __init__(self):
        NBlock = 6
        NTrap = 6
        NTunnel = 3
        NN = NBlock+NTrap+NTunnel

        xx = np.linspace(-15e-7*Unit.cm, 0.0, NN+1)
        Mesh1D.__init__(self, xx, 
                        rgns=[(0, NBlock, 'Block'), 
                              (NBlock, NBlock+NTrap, 'Si3N4'),
                              (NBlock+NTrap, NN, 'Tunnel')],
                        bnds=[(0,'anode'),(NN,'cathode')])
        
        NT = 1e19*pow(Unit.cm, -3.0)
        self.setFieldByFunc('Si3N4', 'NT', lambda x: NT)
    
class TrappingSolver(FVMEqns):
    def __init__(self, device):
        super(TrappingSolver, self).__init__(device)
        
        self.materials = {'Si3N4': Si3N4Material(),
                          'Block': BlockMaterial(),
                          'Tunnel':ToxMaterial(),
                          'Substrate':SubstrateMaterial(-1e17*pow(Unit.cm, -3.0))}

        self.BlockEqn  = InsulatorRegionEqn()
        self.BlockEqn.setMaterial(self.materials['Block'])

        self.Si3N4Eqn = TrappingRegionEqn()
        self.Si3N4Eqn.setMaterial(self.materials['Si3N4'])

        self.ToxEqn = InsulatorRegionEqn()
        self.ToxEqn.setMaterial(self.materials['Tunnel'])

        self.setRegionEqn('Block', self.BlockEqn)
        self.setRegionEqn('Si3N4', self.Si3N4Eqn)
        self.setRegionEqn('Tunnel', self.ToxEqn)

        self.bcAnode   = GateBoundaryEqn()
        self.bcCathode = MOSSubstrateEqn()
        self.bcCathode.setMaterial(self.materials['Substrate'])  # p-type
        self.setBoundaryEqn('anode', self.bcAnode)
        self.setBoundaryEqn('cathode', self.bcCathode)

        self.topIFEqn = DrainIFEqn(True, False)
        self.topIFEqn.setMaterial(self.materials['Si3N4'])
        self.bottomIFEqn = DrainIFEqn(False, True)
        self.bottomIFEqn.setMaterial(self.materials['Si3N4'])
        self.setInterfaceEqn('Block', 'Si3N4', self.topIFEqn)
        self.setInterfaceEqn('Tunnel', 'Si3N4', self.bottomIFEqn)
        
        self.tunneling = NonLocalTunneling(self.device)
        self.tunneling.setInjectBoundary(self.device.getBoundary('cathode'), self.materials['Substrate'], -1)
        self.tunneling.setMaterials(self.materials)
        self.addCustomEqn(self.tunneling)
        
        self.setupEqns()
    
    def setMode(self,mode):
        if mode=='program':
            eqns.tunneling.band = 'Ec'
            self.topIFEqn.setDrainMode(True, False)  # top interface drains electron
            self.bottomIFEqn.setDrainMode(False, True) # bottom interface drains holes
        elif mode=='erase':
            eqns.tunneling.band = 'Ev'
            self.topIFEqn.setDrainMode(False, True)
            self.bottomIFEqn.setDrainMode(True, False)

if __name__ == '__main__':
    device = Trapping()
    eqns = TrappingSolver(device)
    eqns.initGuess()
    
    eqns.bcAnode.setVoltage(0.0)
    eqns.solve()

    #----- voltage source
    def VSource(t):
        Vmax = -20
        trise = 1e-6*Unit.s
        thigh = 10e-6*Unit.s
        tfall = 1e-6*Unit.s
        if t<=trise:
            return t/trise * Vmax
        t-=trise
        
        if t<=thigh:
            return Vmax
        t-=thigh
        
        if t<=tfall:
            return (1.0-t/tfall) * Vmax
        
        return 0.0
    #-----
        

    for ti in xrange(1,21):
        eqns.state.saveTimeStep()
        eqns.state.advanceClock(1e-7*Unit.s)
        vg = VSource(eqns.state.clock)
        print '--------- time:%8g (s),   Vg:%6g (V) ---------' % (eqns.state.clock/Unit.s, vg/Unit.V) 
        eqns.bcAnode.setVoltage(vg)
        if vg>=0:
            eqns.setMode('program')
        else:
            eqns.setMode('erase')

        eqns.solve()
        
        region = device.getRegion('Si3N4')
        Qn=0
        Qp=0
        for cell in region.cells:
            if cell.region.name=='Si3N4':
                Qn += float(eqns.state.getVar(cell.vars[3])) * cell.volume()
                Qp += float(eqns.state.getVar(cell.vars[4])) * cell.volume()
        print '---------- trapped elec: %8g (cm^-2) -------' % (Qn*Unit.cm*Unit.cm)
        print '---------- trapped hole: %8g (cm^-2) -------\n\n' % (Qp*Unit.cm*Unit.cm)
        
    
    iV = device.getVarIdx('Block', 0)
    print eqns.state.x[iV]
    iV = device.getVarIdx('Si3N4', 0)
    print eqns.state.x[iV]
    iV = device.getVarIdx('Tunnel', 0)
    print eqns.state.x[iV]
    
#    iV = device.getVarIdx(0, 0);
#    print eqns.state.x[iV]
#    ielec = device.getVarIdx(0, 1);
#    print eqns.state.x[ielec]
#    ihole = device.getVarIdx(0, 2);
#    print eqns.state.x[ihole]
#    inT = device.getVarIdx(0, 3);
#    print eqns.state.x[ielec]
#    ipT = device.getVarIdx(0, 4);
#    print eqns.state.x[ihole]
        
