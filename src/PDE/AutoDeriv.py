'''
Created on Jun 11, 2009

@author: hash
'''
__all__=['ADVar', 'exp', 'log', 'aux1', 'aux2']

import sys
import math

def _calcDeriv (a, b, func):
    deriv = []
    if not (isinstance(a, ADVar) and isinstance(b, ADVar)):
        return NotImplemented

    i,imax = (0,len(a.deriv))
    j,jmax = (0,len(b.deriv))
    while(i<imax or j<jmax):
        if i<imax:
            ix,dx = a.deriv[i]
        else:
            ix = sys.maxint
        if j<jmax:
            iy,dy = b.deriv[j]
        else:
            iy = sys.maxint

        ii = min(ix,iy)
        if ix<iy:
            deriv.append((ii,func(a.val, b.val, dx,0.0)))
            i=i+1
        elif ix>iy:
            deriv.append((ii,func(a.val, b.val, 0.0,dy)))
            j=j+1
        else:
            deriv.append((ii,func(a.val, b.val, dx,dy)))
            i=i+1
            j=j+1
            
    return deriv

class ADVar(object):
    def __init__(self, val=0, idx=-1):
        self.val = float(val)
        self.deriv = []
        if idx>=0:
            self.deriv.append((idx,1.0))
        
    def __str__(self):
        return 'value:' + str(self.val) + ' deriv: ' + str(self.deriv)

    def __cmp__(self, other):
        if isinstance(other, ADVar):
            return cmp(self.val, other.val)
        else:
            return cmp(self.val, other)
        
    def __int__(self):
        return int(self.val)

    def __long__(self):
        return long(self.val)
    
    def __float__(self):
        return self.val

    def getVal(self):
        return self.val
    
    def setVal(self, v):
        self.val = v
        
    def getDeriv(self, idx=-1):
        if idx<0:
            return self.deriv
        else:
            for i,v in self.deriv:
                if i==idx:
                    return v
            return 0.0
        
    def derivEq(self, other):
        if not self.val == other.val:
            return False
        
        tmp = self-other
        for i,dx in tmp.deriv:
            if not dx==0:
                return False
        return True

    def derivApproxEq(self, other, tol=1e-6):
        absTol = max(abs(self.val), abs(other.val))
        if not abs(self.val-other.val) < absTol:
            return False
        
        tmp = self-other
        for i,dx in tmp.deriv:
            if not abs(dx)<absTol:
                return False
        return True
        
    def __add__(self, other):
        r = ADVar()
        if not isinstance(other, ADVar):
            other = float(other)
            r.val = self.val + other
            r.deriv = self.deriv
            return r
        
        r.val = self.val + other.val
        r.deriv = _calcDeriv(self, other, lambda x,y,dx,dy : dx+dy)
        return r
    
    def __radd__(self,other):
        ''' other + self, other is a scalar'''
        r = ADVar()
        other = float(other)
        r.val = other + self.val
        r.deriv = self.deriv
        return r

    def __sub__(self, other):
        r = ADVar()
        if not isinstance(other, ADVar):
            other = float(other)
            r.val = self.val - other
            r.deriv = self.deriv
            return r
        
        r.val = self.val - other.val
        r.deriv = _calcDeriv(self, other, lambda x,y,dx,dy : dx-dy)
        return r
    
    def __rsub__(self, other):
        ''' other - self, other is a scalar'''
        r = ADVar()
        other = float(other)
        r.val = other - self.val
        for i,dx in self.deriv:
            r.deriv.append((i,-dx))
        return r
    
    def __mul__(self, other):
        r = ADVar()
        if not isinstance(other, ADVar):
            other = float(other)
            r.val = self.val * other
            for i,dx in self.deriv:
                r.deriv.append((i,dx*other))
            return r
        
        r.val = self.val * other.val
        r.deriv = _calcDeriv(self, other, lambda x,y,dx,dy : y*dx + x*dy)
        return r
        
    def __rmul__(self, other):
        ''' other * self, other is a scalar'''
        r = ADVar()
        other = float(other)
        r.val = other * self.val
        for i,dx in self.deriv:
            r.deriv.append((i,other*dx))
        return r
        
    def __div__(self, other):
        r = ADVar()
        if not isinstance(other, ADVar):
            other = float(other)
            r.val = self.val / other
            for i,dx in self.deriv:
                r.deriv.append((i,dx/other))
            return r
        
        r.val = self.val / other.val
        r.deriv = _calcDeriv(self, other, lambda x,y,dx,dy : (y*dx - x*dy)/(y*y) )
        return r
        
    def __rdiv__(self, other):
        ''' other * self, other is a scalar'''
        r = ADVar()
        other = float(other)
        r.val = other / self.val
        for i,dx in self.deriv:
            r.deriv.append( (i, -other*dx/(self.val*self.val)) )
        return r                

    def __pow__(self, other, modulo=0):
        if not modulo==0:
            return NotImplemented
        
        r = ADVar()
        if not isinstance(other, ADVar):
            other = float(other)
            r.val = pow(self.val, other)
            tmp = other * pow(self.val, other-1.0) 
            for i,dx in self.deriv:
                r.deriv.append((i,tmp*dx))
            return r
        
        return NotImplemented
    
    def __abs__(self):
        r = ADVar()
        r.val = abs(self.val)
        
        sign = 0.0
        if self.val>0:
            sign = 1.0
        elif self.val<0:
            sign = -1.0
        
        if not sign==0.0:
            for i,dx in self.deriv:
                r.deriv.append( (i,sign*dx) )
        else:
            for i,dx in self.deriv:
                sign = 0.0
                if dx>0:
                    sign = 1.0
                elif dx<0:
                    sign = -1.0
                r.deriv.append( (i,sign*dx) )

        return r

    def __pos__(self):
        r = ADVar()
        r.val = self.val
        r.deriv = self.deriv
        return r
    
    def __neg__(self):
        r = ADVar()
        r.val = -self.val
        for i, dx in self.deriv:
            r.deriv.append( (i,-dx) )
        return r
    

def exp(x):
    if not isinstance(x, ADVar):
        return math.exp(x)
    
    r = ADVar()
    r.val = math.exp(x)
    for i,dx in x.deriv:
        r.deriv.append( (i, r.val*dx) )
        
    return r

def log(x):
    if not isinstance(x, ADVar):
        return math.log(x)
    
    r = ADVar()
    r.val = math.log(x)
    
    for i,dx in x.deriv:
        if x.val>0 or x.val==0 and dx>=0 :
            r.deriv.append( (i,dx/x.val) )
        else:
            raise ArithmeticError
    return r


_BP0_AUX1  = -8.121635672643270e-03
_BP1_AUX1  =  8.121635672643270e-03
_BP0_DAUX1 = -2.710594333272793e-03
_BP1_DAUX1 =  2.710594333272793e-03
_BP2_DAUX1 = -2.918080926071920e+02
_BP3_DAUX1 =  2.918080926071920e+02
_BP4_DAUX1 = -7.451332191019419e+02
_BP5_DAUX1 =  7.451332191019419e+02
_BP0_AUX2  = -3.673680056967704e+01
_BP1_AUX2  =  3.680808162809191e+01
_BP2_AUX2  =  7.451332191019419e+02
_BP0_DAUX2 = -7.451332191019419e+02
_BP1_DAUX2 = -3.673680056967704e+01
_BP2_DAUX2 =  3.680808162809191e+01
_BP3_DAUX2 =  7.451332191019419e+02
_BP0_MISC  =  7.097827128183643e+02

def aux1(x):
    #
    #                                       x
    #                         Aux1(x) =  -------
    #                                     sinh(x)

    #                    d           sinh(x) - x*cosh(x)
    #                    --Aux1(x) = -------------------
    #                    dx              (sinh(x))^2

    if isinstance(x, ADVar):
        y = x.val
    else:
        y=float(x);
        
    if x < -_BP0_MISC:
        y = -_BP0_MISC
    elif x > _BP0_MISC:
        y =  _BP0_MISC

    if y <= _BP0_AUX1:
        z = y / math.sinh(y)
    elif y <= _BP1_AUX1:
        z = 1 - y*y/6.0*(1.0 - 7.0*y*y/60.0)
    else:
        z = y / math.sinh(y)
        
    if not isinstance(x, ADVar):
        return z
    
    r = ADVar()
    r.val = z
    
    if y <= _BP4_DAUX1:
        pd = 0.0
    elif y <= _BP2_DAUX1:
        pd = -(1+y)*math.exp(y)
    elif y <= _BP0_DAUX1:
        tmp = math.sinh(y); 
        pd = (tmp - y*math.cosh(y))/(tmp*tmp)
    elif y <= _BP1_DAUX1:
        pd = -y/3.0*(1.0 - 7.0*y*y/30.0)
    elif y <= _BP3_DAUX1:
        tmp = math.sinh(y)
        pd = (z - y*math.cosh(y))/(z*z)
    elif y <= _BP5_DAUX1:
        pd =  (1-y)*math.exp(-y)
    else:
        pd = 0.0

    for i,dx in x.deriv:
        r.deriv.append( (i, pd*dx) )

    return r

def aux2(x):
    #                                     1
    #                        Aux2(x) = -------
    #                                  1 + e^x
    #
    #                      d             - e^x
    #                      --Aux2(x) = -----------
    #                      dx          (1 + e^x)^2

    if isinstance(x, ADVar):
        y = x.val
    else:
        y = float(x)
        
    if y <= _BP0_AUX2:
        z = 1.0
    elif y <= _BP1_AUX2:
        z = 1.0 / (1.0 + math.exp(y))
    elif x <= _BP2_AUX2:
        z = math.exp(-y)
    else:
        z = 0.0

    if not isinstance(x, ADVar):
        return z
    
    r = ADVar()
    r.val = z
    if y <= _BP0_DAUX2:
        pd = 0.0
    elif y <= _BP1_DAUX2:
        pd = -math.exp(y)
    elif y <= _BP2_DAUX2:
        tmp1 = math.exp(y)
        tmp2 = tmp1 + 1.0;
        pd = -tmp1/(tmp2*tmp2)
    elif y <= _BP3_DAUX2:
        pd = -math.exp(-x)
    else:
        pd = 0.0

    for i,dx in x.deriv:
        r.deriv.append( (i, pd*dx) )

    return r


if __name__ == '__main__':
    a = ADVar(1.,0)
    b = ADVar(2.,0)
    c = ADVar(3.,1)
    print a
    print b
    print c
    print (a>0)
    print a.getDeriv()
    print a.getDeriv(0)
    print a.getDeriv(1)
    print (a+b)
    print (a+c)
    print (1+c)
