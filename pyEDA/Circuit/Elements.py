__all__=['CircuitElem', 'Resistor', 'Capacitor', 'VSource', 'Diode']

from pyEDA.PDE.NLEqns import *
from pyEDA.PDE.AutoDeriv import *
from pyEDA.PDE.ImplDeriv import *
import math

class CircuitElem(object):
    SymbolPrefix = '?'
    Terminals = []
    AuxVars = []
    
    def __init__(self):
        super(CircuitElem, self).__init__()
        self.varIdx = [0] * (len(self.Terminals)+len(self.AuxVars))
    
    def terminalCount(self):
        return len(self.Terminals)
    
    def varCount(self):
        return len(self.Terminals)+len(self.AuxVars)
    
    def connectToVars(self, varIdx):
        ''' Connect the element to circuit equations
            
            @param: varIdx    list of equation ID numbers that the element is connected to
        '''
        if not len(varIdx) == self.varCount():
            raise ValueError
        
        self.varIdx = varIdx
    
    def calcFunJac(self, state):
        pass


class Resistor(CircuitElem):
    SymbolPrefix = 'R'
    Terminals = ['1', '2']
    
    def __init__(self, res):
        super(Resistor, self).__init__()
        self.res = float(res)
        
    def calcFunJac(self, state):
        V1, V2 = state.getVars(self.varIdx)
        
        state.setFunJac(self.varIdx[0], (V2-V1)/self.res)
        state.setFunJac(self.varIdx[1], (V1-V2)/self.res)
        
class VSource(CircuitElem):
    SymbolPrefix = 'V'
    Terminals = ['+', '-']
    AuxVars = ['i']
    
    def __init__(self, volt):
        super(VSource, self).__init__()
        self.volt = volt

    def setVolt(self, volt):
        self.volt = float(volt)
        
    def calcFunJac(self, state):
        V1, V2, curr = state.getVars(self.varIdx)
        
        state.setFunJac(self.varIdx[0], curr)
        state.setFunJac(self.varIdx[1], -curr)
        state.setFunJac(self.varIdx[2], V1-V2-self.volt)

class Capacitor(CircuitElem):
    SymbolPrefix = 'C'
    Terminals = ['1', '2']

    def __init__(self, cap):
        super(Capacitor, self).__init__()
        self.cap = cap
        
    def calcFunJac(self, state):
        V1, V2 = state.getVars(self.varIdx)
        dv1dt, dv2dt = state.getTimeDerivs(self.varIdx)
        
        state.setFunJac(self.varIdx[0], self.cap*(dv2dt-dv1dt))
        state.setFunJac(self.varIdx[1], self.cap*(dv1dt-dv2dt))
       
class Diode(CircuitElem):
    SymbolPrefix = 'D'
    Terminals = ['+', '-']

    class _DiodeDC(ImplDeriv):
        def __init__(self, Js, Rs):
            super(Diode._DiodeDC,self).__init__(1, 1)
            self.Js = Js
            self.Rs = Rs

        def initGuess(self):
            self.state.setVec([2, -1])

        def calcFunJac(self):
            super(Diode._DiodeDC,self).calcFunJac()
            id,v = self.state.getVars(range(2))
            id1 =  self.Js * ( exp((v-id*self.Rs)/0.0258)  - 1.0 )
            self.state.setFunJac(self.sizeP+0, id-id1)


    def __init__(self, Js=1e-10, Rs=1.0):
        super(Diode, self).__init__()
        self.Js = Js
        self.Rs = Rs

    def _approxSol(self, v):
        if v > 0.6:
            id0 = (v-0.6)/self.Rs
        else:
            id0 = self.Js * ( exp(v/0.0258) - 1.0 )
        return id0

    def calcFunJac(self, state):
        V1, V2 = state.getVars(self.varIdx)

        vd0 = V1.getVal()-V2.getVal()
        id0 = self._approxSol(vd0)

        dcEqn = self._DiodeDC(self.Js, self.Rs)
        dcEqn.setIndepVars([V1-V2])
        dcEqn.state.setVec([id0, vd0])
        dcEqn.solve()

        res = dcEqn.getDeriv()
        curr = res[0]

        state.setFunJac(self.varIdx[0], -curr)
        state.setFunJac(self.varIdx[1], curr)


