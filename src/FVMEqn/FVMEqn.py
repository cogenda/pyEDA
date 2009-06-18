'''
Created on Jun 18, 2009

@author: hash
'''
__all__=['RegionEqn', 'BoundaryEqn', 'InterfaceEqn', 'FVMEqns']

from PDE.NLEqns import *

class RegionEqn(object):
    def __init__(self):
        pass
    
    def eqnPerCell(self):
        return 0
    
    def cellEqn(self, state, cell, dt=0):
        pass
    
    def elemEqn(self, state, elem, dt=0):
        pass
    
    def initGuess(self, state, cell):
        pass
    
    def saveTimeStep(self, dt):
        pass

class BoundaryEqn(object):
    def __init__(self):
        pass

    def cellEqn(self, state, cell, dt=0):
        pass

    def saveTimeStep(self, dt):
        pass
    
class InterfaceEqn(object):
    def __init__(self):
        pass

    def cellPairEqn(self, state, cell1, cell2, dt=0):
        pass

    def saveTimeStep(self, dt):
        pass

class FVMEqns(NLEqns):
    def __init__(self, device):
        super(FVMEqns, self).__init__()

        self.device = device
        self.regionEqns=[None]*len(device.regions)
        self.interfaceEqns=[None]*len(device.interfaces)
        self.boundaryEqns=[None]*len(device.boundaries)
        self.eqnCnt=0

    def setupEqns(self):
        eqnCnt=0
        for r,region in enumerate(self.device.regions):
            eqn = self.regionEqns[r]
            for cell in region.cells:
                eqnPerCell = eqn.eqnPerCell()
                cell.vars = xrange(eqnCnt, eqnCnt+eqnPerCell)
                eqnCnt += eqnPerCell
        self.eqnCnt = eqnCnt
        self.state = NLEqnState(self.eqnCnt)

    def setRegionEqn(self, regionName, eqn):
        if not isinstance(eqn, RegionEqn):
            raise TypeError
        
        for r,region in enumerate(self.device.regions):
            if region.name==regionName:
                self.regionEqns[r]=eqn
                return
        # region name not found...
        raise ValueError

    def setInterfaceEqn(self, r1Name, r2Name, eqn):
        if not isinstance(eqn, InterfaceEqn):
            raise TypeError

        def matchIFName(ifname, r1name, r2name):
            if ifname==r1name+'|'+r2name or ifname==r2name+'|'+r1name:
                return True
            else:
                return False

        for i,interface in enumerate(self.device.interfaces):
            if matchIFName(interface.name, r1Name, r2Name):
                self.interfaceEqns[i]=eqn
                return
        # region name not found...
        raise ValueError

    def setBoundaryEqn(self, bndName, eqn):
        if not isinstance(eqn, BoundaryEqn):
            raise TypeError
        
        for b,boundary in enumerate(self.device.boundaries):
            if boundary.name==bndName:
                self.boundaryEqns[b]=eqn
                return
        # boundary name not found
        raise ValueError
    
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
  