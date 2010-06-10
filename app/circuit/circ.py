from pyEDA.Circuit.Circuit import *
from pyEDA.Circuit.Elements import *

circuit = CircuitEqns()
    
R1 = Resistor(1.0)  # 1ohm resistor
R2 = Resistor(1.0)  # 1ohm resistor
V1 = VSource(0.0)   # V source
D1 = Diode(1e-10)   # diode

#circuit.addElemToCircuit(R1, [0,1])
circuit.addElemToCircuit(D1, [1,0])
circuit.addElemToCircuit(R2, [1,2])
circuit.addElemToCircuit(V1, [2,0])

#            R2
#  2________-----_______ 1
#   |       -----       | 
#   |                   |
#   | +                | | 
#  (V)                 | | R1
#   | -                |_|
#   |                   |
#   |___________________|
#                        0

circuit.setupEqns()
circuit.solve()

for i in xrange(0,5):
    v = i+1
    print 'ramp up to %f V' % v
    V1.setVolt(v)
    circuit.solve()

for node in xrange(0,3):
    i = circuit.nodes[str(node)]
    print 'voltage at node %d \t %g' % (node, circuit.state.x[i]) 
    

