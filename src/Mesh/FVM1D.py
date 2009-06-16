'''
Created on Jun 15, 2009

@author: hash
'''

__all__ = ['Point', 'Node', 'Edge', 'Cell', 'Elem1D', 'Region', 'Boundary', 'Interface']

from FVM import *

class Elem1D(Elem, Edge):
    def __init__(self, cells):
        Edge.__init__(self, cells[0].node, cells[1].node)
        Elem.__init__(self, cells)
        self.__vol = self.len
        self.__pVol = [self.__vol/2.0, self.__vol/2.0]
        cells[0].vol += self.pVolCell(0)
        cells[1].vol += self.pVolCell(1)
        
    def volume(self):
        return self.__vol

    def areaEdge(self, eIdx):
        return 1.0
    
    def pVolCell(self, cIdx):
        return self.__pVol[cIdx]

    def __str__(self):
        return "Elem1D btw "+str(self.cells[0])+" and "+str(self.cells[1])
        
if __name__=='__main__':
    c1 = Cell(Node(0.0))
    c2 = Cell(Node(1.0))
    e = Elem1D([c1,c2])
    print e.volume()
    print e.pVolCell(0)
    print e.pVolCell(1) 