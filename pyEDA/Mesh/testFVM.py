'''
Created on Jun 15, 2009

@author: hash
'''

from FVM1D import *

if __name__=='__main__':
    c1 = Cell(Node(0.0))
    c2 = Cell(Node(1.0))
    e = Elem1D([c1,c2])
    print e.volume()
    print e.pVolCell(0)
    print e.pVolCell(1) 