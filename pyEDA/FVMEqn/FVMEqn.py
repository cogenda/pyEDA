'''
Finite Volume PDE framework.

Equations are classified as Region equation, boundary equation and interface equation.
'''
__all__=['RegionEqn', 'BoundaryEqn', 'InterfaceEqn', 'FVMEqns']

from pyEDA.PDE.NLEqns import *

class RegionEqn(object):
    '''
    Region equation.
    '''
    def __init__(self):
        self.region=None
    
    def eqnPerCell(self):
        ''' Return the number of equations per cell in this region'''
        return 0
    
    def cellEqn(self, state, cell):
        ''' Returns an ADVar of the evaluated equation associated with cell.  '''
        pass
    
    def elemEqn(self, state, elem):
        ''' Returns an ADVar of the evaluated equation associated with the element. '''
        pass
    
    def initGuess(self, state, cell):
        ''' initial guess '''
        pass
    
    def damp(self, state, cell, dx):
        ''' damping '''
        pass
    
class BoundaryEqn(object):
    ''' Boundary equation '''
    def __init__(self):
        self.region=None

    def cellEqn(self, state, cell):
        ''' Returns an ADVar of the evaluated equation associated with boundary cell.  '''
        pass

    
class InterfaceEqn(object):
    ''' Interface Equation '''
    def __init__(self):
        self.region1=None
        self.region2=None

    def cellPairEqn(self, state, cell1, cell2):
        ''' Returns an ADVar of the evaluated equation associated with cell pair cell1, cell2 across the interface.  '''
        pass


class FVMEqns(NLEqns):
    ''' FVM equation solver '''
    def __init__(self, device):
        '''
        constructor

        @param:  device     device mesh of the type Mesh1D
        '''
        super(FVMEqns, self).__init__()

        self.device = device
        self.regionEqns=[None]*len(device.regions)
        self.interfaceEqns=[None]*len(device.interfaces)
        self.boundaryEqns=[None]*len(device.boundaries)
        self.customEqns=[]
        self.eqnCnt=0

    def setupEqns(self):
        '''
        Setup the equations, count the number of equations, init data structure.
        '''
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
        '''
        Attach an equation to the region.

        @param:  regionName  Name of the region
        @param:  eqn         RegionEqn
        '''
        if not isinstance(eqn, RegionEqn):
            raise TypeError
        
        for r,region in enumerate(self.device.regions):
            if region.name==regionName:
                self.regionEqns[r]=eqn
                eqn.region=region
                return
        # region name not found...
        raise ValueError

    def setInterfaceEqn(self, r1Name, r2Name, eqn):
        '''
        Attach an equation to the interface between two regions.

        @param:  r1Name  Name of the first region
        @param:  r2Name  Name of the second region
        @param:  eqn         InterfaceEqn
        '''
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
                eqn.region1=interface.region1
                eqn.region2=interface.region2
                return
        # region name not found...
        raise ValueError

    def setBoundaryEqn(self, bndName, eqn):
        '''
        Attach an equation to the boundary

        @param:  bndName   Name of the Boundary
        @param:  eqn       BoundaryEqn
        '''
        if not isinstance(eqn, BoundaryEqn):
            raise TypeError
        
        for b,boundary in enumerate(self.device.boundaries):
            if boundary.name==bndName:
                self.boundaryEqns[b]=eqn
                eqn.region=boundary.region
                return
        # boundary name not found
        raise ValueError
    
    def addCustomEqn(self, eqn):
        ''' Add a custom equation.  '''
        self.customEqns.append(eqn)
        
    def initGuess(self):
        ''' call the initial guess method of all region cells, elements,
            boundary cells and interface cell pairs.
        '''
        for r,region in enumerate(self.device.regions):
            initGuess = self.regionEqns[r].initGuess
            for cell in region.cells:
                initGuess(self.state, cell)

    def calcFunJac(self):
        ''' 
        Evaluate the function and jacobian at all region cells, elements,
        boundary cells and interface cell pairs.
        Also evaluate the custom equations.
        '''
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
            
        for eqn in self.customEqns:
            eqn(self.state)
  
    def dampStep(self, dx):
        '''
        call the region equations' damp() method.
        '''
        for r,region in enumerate(self.device.regions):
            damp = self.regionEqns[r].damp
            for cell in region.cells:
                damp(self.state, cell, dx)
#        print dx
        return dx
