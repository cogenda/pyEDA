'''
Created on Jun 15, 2009

@author: hash
'''
__all__ = ['Point', 'Node', 'Edge', 'Cell', 'Elem', 'Region', 'Boundary', 'Interface']

class Point(object):
    def __init__(self, pos):
        self.pos = pos
        
    def __str__(self):
        return "point("+str(pos)+")"

class Node(Point):
    def __init__(self, pos):
        Point.__init__(self, pos)
        self.cells = []        

    def __str__(self):
        return "node("+str(self.pos)+")"

class Edge(object):
    def __init__(self, n1, n2):
        self.n1 = n1
        self.n2 = n2
        self.len = abs(n2.pos - n1.pos)
        
    def gradient(self, v1, v2):
        return (v2-v1)/self.len

    def __str__(self):
        return "edge btw "+str(n1)+" and "+str(n2)

############

class Cell(object):
    def __init__(self, node):
        self.node = node
        node.cells.append(self)
        self.region = None
        self.boundaries = []
        self.interfaces = []
        self.vars = []
        self.vol = 0.0
        
    def volume(self): # pure virtual
        return self.vol
    
    def __str__(self):
        return "cell at "+str(self.node) 
    
class Elem(object):
    def __init__(self, cells):
        self.cells = cells
        self.region = None
        
    def volume(self): #pure virtual
        return 0.0

    def areaEdge(self, eIdx):
        return 0.0    #pure virtual
    
    def pVolCell(self, cIdx):
        return 0.0    #pure virtual
    
    def gradient(self, vars):
        return 0.0
        
class Region(object):
    def __init__(self, name=None, material=None):
        self.name = name
        self.material = material
        self.cells = []
        self.elems = []
        
    def __str__(self):
        s = "region "+self.name
        s+="\n cells:"
        for c in self.cells:
            s+=str(c) +","
        s+="\n elems:"
        for e in self.elems:
            s+=str(e) +","
        return s
    
class Boundary(object):
    def __init__(self, region, name=None):
        self.name=name
        self.region=region
        self.cells=[]
        
    def addCell(self, cell, area=1.0):
        self.cells.append((cell,area))

    def __str__(self):
        s = "boundary "+self.name
        s+="\n cells:"
        for c,a in self.cells:
            s+=str(c) +","
        return s
        
class Interface(object):
    def __init__(self, region1, region2, name=None):
        self.name=name
        self.region1 = region1
        self.region2 = region2
        self.cellPairs=[]
        
    def addCellPair(self, c1, c2, area=1.0):
        self.cellPairs.append(((c1,c2),area))

    def __str__(self):
        s = "interface "+self.name
        s+="\n cells:"
        for cp,a in self.cellPairs:
            s+=str(cp[0]) +"|"+ str(cp[1]) + ","
        return s
        
