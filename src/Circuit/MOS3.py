__all__=['MOSLv3']

from PDE.AutoDeriv import *
from Elements import *
import math
import string

q0      = 1.6021918e-19    # C
k       = 1.3806266e-23    # J/K
eps0    = 8.85421487e-12   # F/m
epsSi   = 11.7 * eps0
epsOx   = 3.9 * eps0

class MOSLv3(CircuitElem):
    SymbolPrefix = 'M'
    Terminals = ['D', 'G', 'S', 'B']

    def __init__(self, **kwargs):
        kwArgs={}
        for k,v in kwargs.iteritems():
            kwArgs[string.upper(k)] = v

        self.L      = kwArgs.get('L',       2e-6)   # device length         m
        self.W      = kwArgs.get('W',       50e-6)  # device width          m
        self.LD     = kwArgs.get('LD',      0.0)    # diffusion length      m
        self.WD     = kwArgs.get('WD',      0.0)    # diffusion width       m
        self.TOX    = kwArgs.get('TOX',     1e-7)   # oxide thickness       m
        self.PHI    = kwArgs.get('PHI',     0.6)    # surface inversion potential V
        self.NSUB   = kwArgs.get('NSUB',    0.0)    # substrate doping      cm^-3
        self.XJ     = kwArgs.get('XJ',      0.0)    # junction depth        m
        self.VT0    = kwArgs.get('VT0',     0)      # ref threshold voltage V
        self.GAMMA  = kwArgs.get('GAMMA',   0.0)    # bulk threshold param  V^0.5
        self.DELTA  = kwArgs.get('DELTA',   0.0)    # width effect on Vth   -
        self.ETA    = kwArgs.get('ETA',     0.0)    # static feedback on threshold voltage
        self.U0     = kwArgs.get('U0',      600.0)  # surface mobility      cm^2/V/s
        self.KAPPA  = kwArgs.get('KAPPA',   0.2)    # saturation field factor m
        self.THETA  = kwArgs.get('THETA',   0.0)    # mobility modulation   1/V
        self.VMAX   = kwArgs.get('VMAX',    1e5)    # saturation velocity   m/s
        self.NFS    = kwArgs.get('NFS',     0.0)    # fast surface state    cm^-2
        #self.KP     = kwArgs.get('KP',      2.1e-5) # transconductance      A/V^2
        #self.TPG    = kwArgs.get('TPG',     0.0)    # gate material type    -
        #self.NSS    = kwArgs.get('NSS',     0.0)    # surface state density cm^-2
        self.T      = kwArgs.get('T',       300.15) # temperature

    def _Ids(self, VGS, VDS, VBS):
        # oxide capacitance
        Cox = epsOx / self.TOX

        # effective channel length and width
        Leff = self.L - 2*self.LD
        Weff = self.W - 2*self.WD

        # depletion layer width coefficient
        Xd   = (2.0*epsSi/q0/(self.NSUB*1e6))**0.5

        # built-in voltage
        Vbi  = self.VT0 - self.GAMMA * self.PHI**0.5

        if VBS<=0.0:
            SqVbs = (self.PHI - VBS)**0.5
        else:
            SqVbs = self.PHI**0.5 / ( 1.0 + 0.5/self.PHI * VBS * (1.0 + 0.75/self.PHI * VBS) )

        # short-channel effect correction factor
        c0 =  0.0631353
        c1 =  0.8013292
        c2 = -0.01110777
        T1 = self.XJ * ( c0 + c1 * Xd*SqVbs + c2 * (Xd*SqVbs)**2 )
        Fs = 1.0 - (self.LD + T1)/Leff * ( 1.0 - (Xd*SqVbs/(self.XJ+Xd*SqVbs))**2 )**0.5

        # narrow-width effect correction factor
        Fn = math.pi * epsSi * self.DELTA / 2.0 / Cox / Weff

        # DIBL
        sigma = 8.14e-22 * self.ETA / Leff**3 / Cox

        # Threshold voltage
        Vth = Vbi - sigma * VDS + self.GAMMA * SqVbs * Fs + Fn * SqVbs**2

        # sub-threshold operation
        Xn = 1.0 + q0 * self.NFS*1e4 / Cox + Fn/2.0 + self.GAMMA/2.0 * Fs/SqVbs

        # Modified Vth
        if self.NFS>0:
            Von = Vth + k*self.T/q0 * Xn
        else:
            Von = Vth

        # sub-threshold gate voltage
        if VGS > Von:
            Vgsx = VGS
        else:
            Vgsx = Von

        # bulk charge
        Fb = self.GAMMA/4.0 * Fs/SqVbs + 2.0*Fn

        # surface mobility
        Us = self.U0*1e-4 / ( 1.0 + self.THETA * (Vgsx-Vth) )

        # saturation voltage
        Vsat  = (Vgsx - Vth)/(1 + Fb)
        Vc    = self.VMAX * Leff / Us
        Vdsat = Vsat + Vc - (Vsat**2 + Vc**2)**0.5
        
        if VDS < Vdsat:
            Vdsx = VDS
        else:
            Vdsx = Vdsat

        # velocity saturation
        Fdrain = 1.0 / ( 1.0 + Us*Vdsx/Leff/self.VMAX )

        Ueff = Us * Fdrain
        beta = Weff/Leff * Ueff * Cox

        Ids = beta * (Vgsx - Vth - (1+Fb)/2.0*Vdsx ) * Vdsx  # saturation region

        # channel length modulation
        #Ep = self.VMAX / Us / (1.0-Fdrain)
        Ep = Vc * (Vc+Vdsat) / Leff / Vdsat
        dl = Xd * ( Xd**2 * Ep**2 /4.0  +  self.KAPPA*(VDS-Vdsat) ) ** 0.5 - Ep * Xd**2 /2.0

        if dl < 0.5*Leff:
            lfact = Leff/(Leff-dl)
        else:
            lfact = 4.0*dl/Leff

        Ids = Ids * lfact
        if VGS<Von:
            Ids = Ids * exp(q0/k/self.T * (VGS-Von)/Xn)
        
        return Ids

    def calcFunJac(self, state):
        VD, VG, VS, VB = state.getVars(self.varIdx)
        VDS = VD - VS
        VGS = VG - VS
        VBS = VB - VS

        Ids = self._Ids(VGS, VDS, VBS)


if __name__=='__main__':
    mos = MOSLv3(
            L    =1e-6, 
            W    =1e-6,
            TOX  =2e-8,
            NSUB =1e16,
            XJ   =1e-7,
            VT0  =1.0,
            NFS  =1e10,
            GAMMA=0.5
           )

    #Vgs = 5.0
    #Vbs = 0.0
    #for i in xrange(50):
    #    Vds = 0.1*i
    #    Ids = mos._Ids(Vgs,Vds,Vbs)
    #    print Vds, Ids

    Vds = 0.1
    Vbs = 0.0
    for i in xrange(50):
        Vgs = 0.1*i
        Ids = mos._Ids(Vgs,Vds,Vbs)
        print Vgs, Ids



