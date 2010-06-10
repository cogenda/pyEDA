__all__=['CircuitEqns', 'CircuitElem']

from pyEDA.PDE.NLEqns import *
from Elements import *
import scipy

class CircuitEqns(NLEqns):
    DCMode = 1
    TRMode = 2
    def __init__(self):
        super(CircuitEqns, self).__init__()
        self.mode = self.DCMode
        self.nodes = {} # map {name, idx}
        self.elements = []
        self.varCount = 0
        
    def addElemToCircuit(self, elem, nodes):
        if not isinstance(elem, CircuitElem):
            raise TypeError
        
        if not len(nodes) == elem.terminalCount():
            raise ValueError
        
        vars = []
        for nodeName in nodes:
            nodeKey = str(nodeName)
            if self.nodes.has_key(nodeKey):
                vars.append(self.nodes[nodeKey])
            else:
                idx = self.varCount
                vars.append(idx)
                self.varCount += 1
                self.nodes[nodeKey] = idx
        for av in elem.AuxVars:
            idx = self.varCount
            vars.append(idx)
            self.varCount += 1
        
        self.elements.append(elem)
        elem.connectToVars(vars)
        
    def setupEqns(self):
        self.state = NLEqnState(self.varCount)
    
    def calcFunJac(self):
        for elem in self.elements:
            elem.calcFunJac(self.state)
        
        vGnd = self.nodes['0']
        self.state.resetEqn(vGnd)
        self.state.setFunJac(vGnd, self.state.getVar(vGnd))
        
