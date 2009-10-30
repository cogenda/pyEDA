__all__=['CircuitElem', 'Resistor', 'Capacitor', 'VSource', 'Diode']

from PDE.NLEqns import *
from PDE.AutoDeriv import *
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

    def __init__(self, Js=1e-10, Rs=0.01):
        super(Diode, self).__init__()
        self.Js = Js
        self.Rs = Rs

    def calcFunJac(self, state):
        V1, V2 = state.getVars(self.varIdx)

        v1, v2 = (V1.val, V2.val)
        if v1-v2>0.6:
            vd=ADVar(0.6,0)
        else:
            vd=ADVar(v1-v2,0)

        err = 1e100
        while err>1e-10:
            cc1 = self.Js * (exp(vd/0.025)-1)
            cc2 = (v1-v2-vd)/self.Rs
            err=abs(cc1-cc2)
            deriv = err.getDeriv(0)
            vd -= err.val/deriv

        Vd = V1-V2
        Vd.setVal(vd.val)
        curr = self.Js * (exp(Vd/0.025)-1)

        state.setFunJac(self.varIdx[0], -curr)
        state.setFunJac(self.varIdx[1], curr)


