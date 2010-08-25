''' One-dimensional finite-volume mesh.  '''
__all__=['FVM', 'Mesh1D', 'ElemIterator1D']

import FVM1D as FVM
import scipy
import numpy as np

class ElemIterator1D(object):
    '''
    Iterator to walk through the mesh's elements.
    '''
    def __init__(self, mesh, dir=1):
        '''
        constructor
        
        @param: mesh    Mesh1D
        @param: dir    -1 to walk to the left, +1 to walk to the right
        '''
        self.mesh = mesh
        self.dir  = dir
        
        if dir>0: # from left to right
            self.curr_r = 0
            self.curr_e = 0
        else:
            self.curr_r = len(mesh.regions)-1
            region = mesh.regions[self.curr_r]
            self.curr_e = len(region.elems)-1
        
    def setElem(self, elem_cell_bnd_or_interface, dir):
        '''
        Set the current element using an element, cell, boundary or interface.
        Also set walking direction.
        '''
        self.dir = dir
        if isinstance(elem_cell_bnd_or_interface, FVM.Elem1D):
            elem = elem_cell_bnd_or_interface
            for r,region in enumerate(self.mesh.regions):
                for i,e in enumerate(region.elems):
                    if e==elem:
                        if dir<0 and i>0:
                            self.curr_r = r
                            self.curr_e = i-1
                            return True
                        elif dir>0 and i+1<len(region.elems):
                            self.curr_r = r
                            self.curr_e = i+1
                            return True
                        else:
                            return False
            return False
        elif isinstance(elem_cell_bnd_or_interface, FVM.Cell):
            cell = elem_cell_bnd_or_interface
            for r,region in enumerate(self.mesh.regions):
                for i,e in enumerate(region.elems):
                    if (dir>0 and e.cells[0]==cell) or (dir<0 and e.cells[1]==cell):
                        self.curr_r = r
                        self.curr_e = i
                        return True
            return False
        elif isinstance(elem_cell_bnd_or_interface, FVM.Interface):
            interface = elem_cell_bnd_or_interface
            (c1,c2),dummy = interface.cellPairs[0]
            for i in self.mesh.interfaces:
                if interface==i:
                    if dir<0:
                        return self.setElem(c1, dir) # going left
                    else:
                        return self.setElem(c2, dir) # going right
            return False
        elif isinstance(elem_cell_bnd_or_interface, FVM.Boundary):
            bnd = elem_cell_bnd_or_interface
            cell, dummy = bnd.cells[0]
            for b in self.mesh.boundaries:
                if bnd==b:
                    return self.setElem(cell, dir)
            return False        

    def __iter__(self):
        return self
    
    def next(self):
        '''
        Return the current element, and advance to the next in the preset direction.
        '''
        mesh=self.mesh
        dir =self.dir
        if self.curr_r<0 or self.curr_r >= len(mesh.regions):
            raise StopIteration
        region = mesh.regions[self.curr_r]
        if self.curr_e<0 or self.curr_e >= len(region.elems):
            raise StopIteration
        elem = region.elems[self.curr_e]

        if dir>0:
            if self.curr_e+1 < len(region.elems):
                self.curr_e +=1
            else:
                self.curr_r +=1
                self.curr_e = 0
        else:
            if self.curr_e > 0:
                self.curr_e -=1
            elif self.curr_r > 0:
                self.curr_r -=1
                self.curr_e = len(mesh.regions[self.curr_r].elems)-1
            else:
                self.curr_r=-1

        return elem
        
class Mesh1D(object):
    def __init__(self, xx, rgns, bnds):
        ''' 
        Construct the mesh structure, setup cells, elements, boundaries and interfaces.
        @param: xx : coordinates of nodes (list or array)
        @param: rgns list of regions, each region expressed as tuple [ (istart, iend, name), (...), ...]
                example::
                    [ (0, 10, 'region1'), (10, 20, 'region2') ]
        @param: bnds list of boundaries, [ (nodeNo, name), (...), ... ]
        '''
        
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
            interface = FVM.Interface(rgn1, rgn2, rgn1.name+'|'+rgn2.name)
            c1 = rgn1.cells[len(rgn1.cells)-1]
            c2 = rgn2.cells[0]
            interface.addCellPair(c1, c2)
            self.interfaces.append(interface)
            
        for node,name in bnds:
            #boundary nodes can only belong one region, hence one cell
            cell=nodes[node].cells[0]
            boundary = FVM.Boundary(cell.region, name)
            boundary.addCell(cell) 
            self.boundaries.append(boundary)
                
    
    def getRegion(self, rName_or_rIdx):
        '''
        Return the region specified by region name of index.
        '''
        if isinstance(rName_or_rIdx, str):
            for r in self.regions:
                if r.name==rName_or_rIdx:
                    return r
            return None
        else:
            rIdx = int(rName_or_rIdx)
            if not rIdx<len(self.regions):
                raise IndexError
            return self.regions[rIdx]

    def setRegionMaterial(self, rName_or_rIdx, material):
        '''
        Short-hand method for setting the material of a region
        '''
        self.getRegion(rName_or_rIdx).material = material
    
    def getBoundary(self, bName_or_bIdx):
        '''
        Return the boundary specified by region name of index.
        '''
        if isinstance(bName_or_bIdx, str):
            for b in self.boundaries:
                if b.name==bName_or_bIdx:
                    return b
            return None
        else:
            bIdx = int(bName_or_bIdx)
            if not bIdx<len(self.boundaries):
                raise IndexError
            return self.boundaries[bIdx]

    def getInterface(self, r1, r2):
        '''
        Return the interface specified by the two region names
        '''
        if not (isinstance(r1, FVM.Interface) and isinstance(r2, FVM.Interface)):
            raise TypeError

        for i in self.interfaces:
            (c1,c2), dummy = i.cellPairs[0]
            if c1.region==r1 and c2.region==r2:
                return i
            elif c1.region==r2 and c2.region==r1:
                return i
        return None

    def setFieldByFunc(self, rName_or_rIdx, name, func):
        '''
        Set the field using a custom function.

        @param: rName_or_rIdx   region name or index
        @param: name            field name
        @param: func            func(pos), function will be supplied by node coordinates, 
                                should return field value.
        '''
        region=self.getRegion(rName_or_rIdx)        
        for c in region.cells:
            c.fields[name] = func(c.node.pos)
                
    def getField(self, rName_or_rIdx, name, cIdx=None):
        '''
        Return the field value at a node, or for all node
        @param: rName_or_rIdx   region name or index
        @param: name            field name
        @param: cIdx            cell index, None if want values at all cells.
        '''
        region=self.getRegion(rName_or_rIdx)        
        if not cIdx==None:
            if not cIdx<len(region.cells):
                raise IndexError
            return region.cells[cIdx].fields[name]
        else:
            vec = np.zeros(len(region.cells))
            for i,c in enumerate(region.cells):
                vec[i] = c.fields[name]
            return vec
    
    def getVarIdx(self, rName_or_rIdx, var):
        '''
        Return the variable index for a type of variables.
        @param: rName_or_rIdx   region name or index
        @param: variable offset with-in the cell
        '''
        region=self.getRegion(rName_or_rIdx)        
        vec = np.zeros(len(region.cells), dtype=int)
        for i,c in enumerate(region.cells):
            vec[i] = c.vars[var]
        return vec

if __name__=='__main__':
    xx = np.linspace(0,1,11)
    
    mesh = Mesh1D(xx, 
                  [(0,5,'haha'), (5,10,'xixi')],
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
    #print mesh.getVarIdx(1, 0)
    print "----------------"
    iter = ElemIterator1D(mesh)
    for e in iter:
        print e
        
    print "----------------"
    iter = ElemIterator1D(mesh)
    iter.setElem(mesh.interfaces[0],-1)
    for e in iter:
        print e

    print "----------------"
    iter = ElemIterator1D(mesh)
    iter.setElem(mesh.regions[0].elems[2],1)
    for e in iter:
        print e
    
