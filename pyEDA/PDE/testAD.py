'''
Created on Jun 12, 2009

@author: hash
'''
import unittest
from AutoDeriv import *
import math

class Test(unittest.TestCase):

    def testCmp(self):
        a = ADVar(1.0, 0)
        b = ADVar(1.0, 1) 
        c = ADVar(2.0, 0)
        
        self.assert_( a==b )
        self.assert_( not a==c ) 
        self.assert_( not b==c )
        
        self.assert_( a<=b )
        self.assert_( a<=c )
        self.assert_( a<c )
        
    def testCmpDeriv(self):
        a = ADVar(1.0, 0)
        b = ADVar(1.0, 1)
        c = ADVar(1.0, 0)
        d = ADVar(1.0, 2)
        e = ADVar(-1.0, 1)
        
        self.assert_(a.derivEq(c))
        self.assert_((a+b).derivEq(b+c))
        self.assert_((a+d+b).derivEq(b+d+c))
        self.assert_(a.derivApproxEq(c))
        self.assert_((a+b).derivApproxEq(b+c))
        self.assert_((a+d+b).derivApproxEq(b+d+c))
        
        self.assertFalse(a.derivEq(b))
        self.assertFalse(e.derivEq(b))
        self.assertFalse(a.derivApproxEq(b))
        self.assertFalse(e.derivApproxEq(b))
        
        
    def testAdd(self):
        a = ADVar(0.0, 0)
        b = ADVar(1.0, 1)
        c = ADVar(3.0, 3)
        d = ADVar(5.0, 5)
        e = ADVar(101, 101)
        f = ADVar(1, 1)
        
        s1 = a+b+c+d+e+f
        s2 = f+e+d+c+b+a
        self.assert_(s1.derivEq(s2))
        self.assert_(s1==111)
        
        dl = s1.getDeriv()
        self.assert_(len(dl)==5)
        self.assert_(s1.getDeriv(1)==2.0)
        
        bp1 = 1+b
        self.assert_(bp1==2.0)
        self.assert_(len(bp1.getDeriv())==1)
        self.assert_(bp1.getDeriv(1)==1.0)
        
        self.assert_((a+4).derivEq((a+2)+2))
        self.assert_((a+4).derivApproxEq((a+2)+2))
        
    def testSub(self):
        a = ADVar(0.0, 0)
        b = ADVar(1.0, 1)
        c = ADVar(3.0, 3)
        d = ADVar(5.0, 5)
        e = ADVar(101, 101)
        f = ADVar(1, 1)

        s1 = e-a-b-c-d
        s2 = e-d-c-a-b
        s3 = e-a-f-c-d
        s4 = +e+(-a)+(-b)+(-c)+(-d)
        
        self.assert_(s1.derivEq(s2))
        self.assert_(s1.derivEq(s3))
        self.assert_(s1.derivEq(s4))
        self.assert_(e.derivEq(s1+a+b+c+d))
        
        self.assert_((e+10).derivEq((e-1)+11))
        self.assert_((e+10).derivApproxEq((e-1)+11))
        
    def testMul(self):
        a = ADVar(4.0, 0)
        b = ADVar(1.0, 1)
        c = ADVar(3.0, 3)
        d = ADVar(5.0, 5)
        e = ADVar(101, 101)
        f = ADVar(1, 1)

        p1 = a*b*c*d*e*f
        p2 = f*e*d*c*b*a
        self.assert_(p1.derivApproxEq(p2))
        
        self.assert_(a*c==12.0)
        self.assert_((a*c).getDeriv(0)==3.0)
        self.assert_((a*c).getDeriv(3)==4.0)
        
        self.assert_((3*c).derivApproxEq(c*3.0))
        self.assert_((3*c).derivApproxEq(c+c+c))        
        
    def testDiv(self):
        a = ADVar(4.0, 0)
        b = ADVar(1.0, 1)
        c = ADVar(3.0, 3)
        d = ADVar(5.0, 5)
        e = ADVar(101, 101)
        f = ADVar(1, 1)
        
        q1 = b/a
        q2 = c/a
        q3 = (b+c)/a
        self.assert_(q1==0.25)
        self.assertAlmostEqual(q1.getDeriv(0),-1.0/16)
        self.assertAlmostEqual(q1.getDeriv(1),1.0/4)

        q1 = b/d
        q2 = c/d
        q3 = (b+c)/d
        self.assert_(q3.derivApproxEq(q1+q2))
        
        q1 = e/d/c/f
        q2 = a/d/c/f
        q3 = (a+e)/(d*c*f)
        self.assert_(q3.derivApproxEq(q1+q2))

    def testPow(self):
        b = ADVar(1.0, 1)
        c = ADVar(3.0, 3)
        d = ADVar(5.0, 5)
        a = b*c*d

        p1 = pow(a,6.6)

        p2 = pow(b,6.6)
        p3 = pow(c,6.6)
        p4 = pow(d,6.6)
        
        self.assertAlmostEqual(p1, p2*p3*p4)
        self.assertNotAlmostEqual(p1, p2*p3)
        
        
    def testExp(self):
        a = ADVar(1.0, 0)
        b = ADVar(1.0, 1)
        c = ADVar(1.0, 2)
        
        e = exp(a+b+c)
        f = exp(a)*exp(b)*exp(c)
        self.assert_(e.derivApproxEq(f))
        self.assertAlmostEqual(e, math.exp(3))
        self.assertAlmostEqual(e.getDeriv(0), math.exp(3))
        self.assertAlmostEqual(e.getDeriv(1), math.exp(3))
        self.assertAlmostEqual(e.getDeriv(2), math.exp(3))

    def testSin(self):
        a = ADVar(0.25*math.pi, 0)
        b = ADVar(0.5*math.pi, 1)
        c = ADVar(math.pi/3.0, 2)
        d = sin(a) + sin(b) + sin(c)
        
        self.assertAlmostEqual(d, math.sin(a)+math.sin(b)+math.sin(c))
        self.assertAlmostEqual(d.getDeriv(0), math.cos(a))
        self.assertAlmostEqual(d.getDeriv(1), math.cos(b))
        self.assertAlmostEqual(d.getDeriv(2), math.cos(c))

    def testCos(self):
        a = ADVar(0.25*math.pi, 0)
        b = ADVar(0.5*math.pi, 1)
        c = ADVar(math.pi/3.0, 2)
        d = cos(a) + cos(b) + cos(c)
        
        self.assertAlmostEqual(d, math.cos(a)+math.cos(b)+math.cos(c))
        self.assertAlmostEqual(d.getDeriv(0), -math.sin(a))
        self.assertAlmostEqual(d.getDeriv(1), -math.sin(b))
        self.assertAlmostEqual(d.getDeriv(2), -math.sin(c))

    def testLog(self):
        a = ADVar(1.0, 0)
        b = ADVar(2.0, 1)
        c = ADVar(3.0, 2)
        
        e = log(a*b*c)
        f = log(a)+log(b)+log(c)
        self.assert_(e.derivApproxEq(f))
 
    def testSqrt(self):
        a = ADVar(4.0, 0)
        b = ADVar(2.0, 0)
        
        e = sqrt(a)
        self.assertAlmostEqual(e, 2.0)
        self.assertAlmostEqual(e.getDeriv(0), 0.25)

    def testAux1(self):
        x = ADVar(2.4, 0) - ADVar(-1.1, 1)
        
        a1 = float(x)/math.sinh(float(x))
        a2 = aux1(x)
        self.assertAlmostEqual(a1, float(a2))
        
        self.assertEqual(aux1(ADVar(0,0)).getVal(),1.0)
        self.assertEqual(aux1(ADVar(0,0)).getDeriv(0),0)
        
    def testAux2(self):
        x = ADVar(2.4, 0) - ADVar(-1.1, 1)
        a1 = 1.0/(1.0+math.exp(float(x)))
        a2 = aux2(x)

        self.assertAlmostEqual(a1, float(a2))
        
        self.assertEqual(aux2(ADVar(0,0)).getVal(), 0.5)
        self.assertEqual(aux2(ADVar(0,0)).getDeriv(0), -0.25)

    def testErf(self):
        x = ADVar(0.0, 0)
        self.assertAlmostEqual(erf(x).getVal(),  0.0)
        self.assertAlmostEqual(erfc(x).getVal(), 1.0)

        x = ADVar(0.1, 0)
        self.assertAlmostEqual(erf(x).getVal(),  0.1124629)
        self.assertAlmostEqual(erfc(x).getVal(), 0.8875371)

        x = ADVar(1.234, 0)
        y = erf(x)+erfc(x)
        self.assertAlmostEqual(y.getVal(), 1.0)
        self.assertAlmostEqual(y.getDeriv(0), 0.0)
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testCmp']
    unittest.main()
