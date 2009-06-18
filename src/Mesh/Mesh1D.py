'''
Created on Jun 15, 2009

@author: hash
'''
__all__=['FVM', 'Mesh1D']

import FVM1D as FVM
import scipy
import numpy as np

class Mesh1D(object):
    def __init__(self, xx, rgns, bnds):
        # xx : coordinates of nodes
        # rgns : [ (istart, iend, name), (...), ...]
        # bnds : [ (nodeNo, name), (...), ... ]
        
        bnds.sort( lambda bnd1, bnd2 : cmp(bnd1[0], bnd2[0]) )
        rgns.sort( lambda r1, r2: cmp(r1[0], r2[0]) )

        self.regions = []
        self.boundaries = []
        self.interfaces=[]
        
        nodes=[]
        for x in xx:
            node = FVM.Node(x)
            nodes.append(node)
        
        eqnCnt=0
        rgnCnt=len(rgns)
        for r in xrange(0,rgnCnt):
            istart,iend,name = rgns[r]
            region = FVM.Region(name)
            for i in xrange(istart, iend+1):
                cell = FVM.Cell(nodes[i])
                cell.fields = {}
                cell.region = region
                region.cells.append(cell)
    
            j=0
            for i in xrange(istart, iend):
                c1 = region.cells[j]
                c2 = region.cells[j+1]
                elem = FVM.Elem1D([c1,c2])
                elem.region = region
                region.elems.append(elem)
                j+=1

            self.regions.append(region)
        
        for r in xrange(0, rgnCnt-1):
            rgn1 = self.regions[r]
            rgn2 = self.regions[r+1]
            interface = FVM.Interface(rgn1.name+'|'+rgn2.name)
            c1 = rgn1.cells[len(rgn1.cells)-1]
            c2 = rgn2.cells[0]
            interface.addCellPair(c1, c2)
            self.interfaces.append(interface)
            
        for node,name in bnds:
            boundary = FVM.Boundary(name)
            #boundary nodes can only belong one region, hence one cell
            boundary.addCell(nodes[node].cells[0]) 
            self.boundaries.append(boundary)
        
    def setFieldByFunc(self, rIdx, name, func):
        if not rIdx<len(self.regions):
            raise IndexError
        region = self.regions[rIdx]
        
        for c in region.cells:
            c.fields[name] = func(c.node.pos)
                
    def getField(self, rIdx, name, cIdx=None):
        if not rIdx<len(self.regions):
            raise IndexError
        region = self.regions[rIdx]
        
        if not cIdx==None:
            if not cIdx<len(region.cells):
                raise IndexError
            return region.cells[cIdx].fields[name]
        else:
            vec = np.zeros(len(region.cells))
            for i,c in enumerate(region.cells):
                vec[i] = c.fields[name]
            return vec
    
    def getVarIdx(self, rIdx, var):
        if not rIdx<len(self.regions):
            raise IndexError
        region = self.regions[rIdx]
        vec = np.zeros(len(region.cells), dtype=int)
        for i,c in enumerate(region.cells):
            vec[i] = c.vars[var]
        return vec

if __name__=='__main__':
    xx = np.linspace(0,1,11)
    
    mesh = Mesh1D(xx, 
                  [(0,5,'haha',2), (5,10,'xixi',2)],
                  [(0,'A'),(10,'B')] )

    for r in mesh.regions:
        print r
    for b in mesh.boundaries:
        print b
    for i in mesh.interfaces:
        print i
        
    mesh.setFieldByFunc(0, 'electron', lambda x: 100.0)
    mesh.setFieldByFunc(0, 'hole', lambda x: 0.01)
    mesh.setFieldByFunc(1, 'pairs', lambda x: 3.14)

    print mesh.getField(0, 'hole')
    print mesh.getVarIdx(1, 0)
    