__all__=['MOSLv3']

from pyEDA.PDE.AutoDeriv import *
from Elements import *
import math
import string

q0      = 1.6021918e-19    # C
kb      = 1.3806266e-23    # J/K
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
        self.T      = kwArgs.get('T',       300.15) # temperature
        self.LD     = kwArgs.get('LD',      0.0)    # diffusion length      m
        self.WD     = kwArgs.get('WD',      0.0)    # diffusion width       m
        self.TOX    = kwArgs.get('TOX',     1e-7)   # oxide thickness       m

        if kwArgs.has_key('NSUB'):
            self.NSUB   = kwArgs.get('NSUB')        # substrate doping      cm^-3

            T = self.T
            Eg = 1.16 - 7.02e-4 * T*T/(T+1108)

            n_i = 1.45e10 * (self.T/300)**1.5 * exp(q0*Eg/2.0/kb * (1./300. - 1./self.T))
            self.PHI    = 2*kb*T/q0*log(self.NSUB/n_i)
            Cox = epsOx / self.TOX
            self.GAMMA  = (2*q0*epsSi*self.NSUB*1e6)**0.5 / Cox
        else:
            self.NSUB   = 0.0
            self.PHI    = kwArgs.get('PHI',     0.6)# surface inversion potential V
            self.GAMMA  = kwArgs.get('GAMMA',   0.0)# bulk threshold param  V^0.5

        self.XJ     = kwArgs.get('XJ',      0.0)    # junction depth        m
        self.VT0    = kwArgs.get('VT0',     0)      # ref threshold voltage V
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

        self.RSDW   = kwArgs.get('RSDW',    0.0)    # source/drain resistance of unit width ohm.m


        if self.XJ<1e-8:
            self.XJ = abs(self.XJ)+1e-9

        self.RSDW = abs(self.RSDW)
        self.NFS   = abs(self.NFS)
        self.ETA   = abs(self.ETA)
        self.KAPPA = abs(self.KAPPA)

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
            Von = Vth + kb*self.T/q0 * Xn
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
        if Vgsx>Vth:
            Fdrain = 1.0 / ( 1.0 + Vdsx/Vc )
        else:
            Fdrain = 1.0

        Ueff = Us * Fdrain
        beta = Weff/Leff * Ueff * Cox

        Ids = beta * (Vgsx - Vth - (1+Fb)/2.0*Vdsx ) * Vdsx  # saturation region

        # channel length modulation
        #Ep = self.VMAX / Us / (1.0-Fdrain)
        if VDS>Vdsat and Vgsx>Vth:
            Ep = Vc * (Vc+Vdsat) / Leff / Vdsat
            dl = Xd * ( Xd**2 * Ep**2 /4.0  +  self.KAPPA*(VDS-Vdsat) ) ** 0.5 - Ep * Xd**2 /2.0

            if dl < 0.5*Leff:
                lfact = Leff/(Leff-dl)
            else:
                lfact = 4.0*dl/Leff
        else:
            lfact = 1.0

        Ids = Ids * lfact
        if VGS<Von:
            Ids = Ids * exp(q0/kb/self.T * (VGS-Von)/Xn)
       
        # series resistance
        Rsd = self.RSDW / self.W
        if abs(VDS)>1e-6:
            Ids = Ids/(1. + Rsd*Ids/VDS)
        return Ids

    def calcFunJac(self, state):
        VD, VG, VS, VB = state.getVars(self.varIdx)
        VDS = VD - VS
        VGS = VG - VS
        VBS = VB - VS

        Ids = self._Ids(VGS, VDS, VBS)


if __name__=='__main__':
    mos = MOSLv3(
            L    =1.200e-6, 
            W    =10.00e-6,
            LD   =1.647e-7,
            WD   =1.000e-7,
            TOX  =2.120e-8,
            PHI  =0.6,
            NSUB =2.747e16,
            XJ   =2.000e-7,
            VT0  =0.7860,
            GAMMA=0.5863,
            DELTA=0.6967,
            ETA  =4.368e-2,
            U0   =591.7,
            KAPPA=1.396e-10,
            THETA=8.122e-2,
            VMAX =1.733e5,
            NFS  =1.98e12
           )

    for Vbs in [0, -1.0, -2.0] :
        for j in xrange(26):
            Vds = 0.2*j
            for k in xrange(26):
                Vgs = 0.2*k
                Ids = mos._Ids(Vgs,Vds,Vbs)
                print Vgs, Vds, Vbs, Ids



