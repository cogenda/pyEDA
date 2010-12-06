__all__=['MOSBSim3v3']

from pyEDA.PDE.AutoDeriv import *
from Elements import *
import math
import string

q0      = 1.6021918e-19    # C
kb      = 1.3806266e-23    # J/K
kboq    = 8.617087e-5      # kb/q0
eps0    = 8.85421487e-12   # F/m
epsSi   = 11.7 * eps0
epsOx   = 3.9 * eps0

def log1pexp(x):
    if x<37.0:
        return log(1.+exp(x))
    else:
        return x


class MOSBSim3v3(CircuitElem):
    SymbolPrefix = 'M'
    Terminals = ['D', 'G', 'S', 'B']

    # {{{
    Param_Defs = {
        "L":  1.0e-6,           # Length (m)
        "W":  1.0e-6,           # Width (m)
        "TOX":  1.5e-8,         # Gate oxide thickness (m)
        "TOXM":  None,          # Gate oxide thickness used in extraction
        "VTH0":  0.7,           # Threshold voltage @VBS=0 for large L
        "CDSC":  2.4e-4,        # Drain/Source and channel coupling capacitance (F/m^2)
        "CDSCB":  0.0,          # Body-bias dependence of cdsc (F/V/m^2)
        "CDSCD":  0.0,          # Drain-bias dependence of cdsc (F/V/m^2)
        "CIT":  0.0,            # Interface state capacitance
        "NFACTOR":  0.0,        # Subthreshold swing coefficient
        "XJ":  1.5e-7,          # Junction depth (m)
        "VSAT":  8.0e4,         # Saturationvelocity at tnom (m/sec)
        "AT":  3.3e4,           # Temperature coefficient of vsat (m/sec)
        "A0":  1.0,             # Non-uniform depletion width effect coefficient
        "AGS":  0.0,            # Gate bias coefficient of Abulk
        "A1":  0.0,             # Non-saturation effect coefficient
        "A2":  1.0,             # Non-saturation effect coefficient
        "KETA":  -0.047,        # Body-bias coefficient of non-uniform depletion width effect
        "NSUB":  6e16,          # Substrate doping concentration
        "NCH":  1.7e17,         # Channel doping concentration
        "GAMMA1":  None,        # Body effect coefficient near the surface
        "GAMMA2":  None,        # Body effect coefficient in the bulk
        "VBX":  None,           # Vbs at which the depletion width = XT
        "NGATE":  0.0,          # Poly-gate doping concentration (cm^-3)
        "VBM":  -3.0,           # Maximum body voltage (V)
        "XT":  1.55e-7,         # Doping depth
        "KT1":  -0.11,          # Temperature coefficient of Vth (V)
        "KT1L":  0.0,           # Temperature coefficient of Vth, channel length dependence (V.m)
        "KT2":  0.022,          # Body-coefficient of kt1 (-)
        "K3":  80.0,            # Narrow width effect coefficient
        "K3B":  0.0,            # Body effect coefficient of k3
        "W0":  2.5e-6,          # Narrow width effect parameter (m)
        "NLX":  1.74e-7,        # Lateral non-uniform doping effect (m)
        "DVT0":  2.2,           # Short channel effect coefficient 0 (-)
        "DVT1":  0.53,          # Short channel effect coefficient 1 (-)
        "DVT2":  -0.032,        # Short channel effect coefficient 2 (V^-1)
        "DVT0W":  0.0,          # Narrow width effect coefficient 0 (m^-1)
        "DVT1W":  5.3e6,        # Narrow width effect coefficient 1 (m^-1)
        "DVT2W":  -0.032,       # Narrow width effect body-bias coefficient (V^-1)
        "DROUT":  0.0,          # DIBL coefficient of output resistance
        "DSUB":  None,          # DIBL coefficient in the subthreshold region
        "UA":  2.25e-9,         # Linear gate dependence of mobility (m/V)
        "UB":  5.87e-19,        # Quadratic gate dependence of mobility (m/V)^2
        "UC":  -4.65e-11,       # Body-bias dependence of mobility (m/V^2)
        "U0":  670.0,           # Low-field mobility at Tnom (cm^2/V/s)
        "VOFF":  -0.08,         # Threshold voltage offset (V)
        "TNOM": 25.0,           # Parameter measurement temperature (deg C)
        "ELM":  0.0,            # Non-quasi-static Elmore Constant Parameter
        "DELTA":  0.01,         # Effective Vds parameter (V)
        "RDSW":  0.0,           # Sorce-drain resistance per width
        "PRWG":  0.0,           # Gate-bias effect on parasitic resistance
        "PRWB":  0.0,           # Body-effect on parasitic resistance
        "PRT":  0.0,            # Temperature coefficient of parasitic resistance
        "ETA0":  0.08,          # Subthreshold region DIBL coefficeint
        "ETAB":  -0.07,         # Subthreshold region DIBL body-bias coefficeint (1/V)
        "PCLM":  1.3,           # Channel length modulation coefficient (-)
        "PDIBLC1":  0.39,       # Drain-induced barrier lowering oefficient
        "PDIBLC2":  0.0086,     # Drain-induced barrier lowering oefficient
        "PDIBLCB":  0.0,        # Body-effect on drain induced barrier lowering
        "PSCBE1":  4.24e8,      # Substrate current body-effect coeffiecient (V/m)
        "PSCBE2":  1.0e-5,      # Substrate current body-effect coeffiecient (m/V)
        "PVAG":  0.0,           # Gate dependence of output resistance parameter (-)
        "VFB":  -1.0,           # Flat band voltage
        "ACDE":  0.0,           # Exponential coefficient for finite charge thickness
        "MOIN":  0.0,           # Coefficient for gate-bias dependent surface potential
        "NOFF":  0.0,           # C-V turn-on/off parameter
        "VOFFCV":  0.0,         # C-V lateral shift parameter
        "LINT":  0.0,           # Length reduction parameter (m)
        "LL":  0.0,             # Length reduction parameter
        "LLC":  None,           # Length reduction parameter for CV
        "LLN":  1.0,            # Length reduction parameter
        "LW":  0.0,             # Length reduction parameter
        "LWC":  None,           # Length reduction parameter for CV
        "LWN":  1.0,            # Length reduction parameter
        "LWL":  0.0,            # Length reduction parameter
        "LWLC":  0.0,           # Length reduction parameter for CV
        "WR":  1.0,             # Width dependence of rds (-)
        "WINT":  0.0,           # Width reduction parameter (m)
        "DWG":  0.0,            # Width reduction gate bias dependence (m/V)
        "DWB":  0.0,            # Width reduction body bias dependence (m/V)
        "WL":  0.0,             # Width reduction parameter
        "WLC":  None,           # Width reduction parameter for CV
        "WLN":  1.0,            # Width reduction parameter
        "WW":  0.0,             # Width reduction parameter
        "WWC":  None,           # Width reduction parameter for CV
        "WWN":  1.0,            # Width reduction parameter
        "WWL":  0.0,            # Width reduction parameter
        "WWLC": None,           # Width reduction parameter for CV
        "B0":  0.0,             # Abulk narrow width parameter
        "B1":  0.0,             # Abulk narrow width parameter
        "CLC":  0.0,            # Vdsat paramater for C-V model
        "CLE":  0.0,            # Vdsat paramater for C-V model
        "ALPHA0":  0.0,         # Substrate current model parameter (m/V)
        "ALPHA1":  0.0,         # Substrate current model parameter (1/V)
        "BETA0":  30.0,         # Diode limiting current (V)
        "UTE":  -1.5,           # Temperature exponent of mobility (-)
        "K1":  0.5,             # First order body effect coefficient (V^0.5)
        "K2":  0.0,             # Second order body effect coefficient (-)
        "TEMP": 100,            # Circuit temperature
        "UA1":  4.31e-9,        # Temperature coefficient for ua in m/V
        "UB1":  -7.61e-18,      # Temperature coefficient for ub in (m/V)^2
        "UC1":  -5.6e-11,       # Temperature coefficient for uc in m/V^2
        "IJTH": 0.1,            # Diode limiting current (A)
        "JS": 1e-4,            # Source drain junction saturation current per unit area (A/m^2)
        "JSSW": 0.0,           # Side wall saturation current density (A/m)
        "RSH": 0.0,             # Source drain sheet resistance (ohm/square)
        "NJ": 1.0,              # Emission coefficient of junction
        "XTI": 3.0,             # Junction current temperature exponent coefficient
        "HDIF": 0.0,            # length of heavy diffusion region (m)
        "WMLT": 1.0,            # width multiplier
        "A2KETA": 0.0,
        "XL": 0.0,
        "XW": 0.0,
        "AS": 0.0,
        "AD": 0.0,
        "PS": 0.0,
        "PD": 0.0,
        "BINUNIT": 1,
        "CALCACM": 1,
        }
    # }}}

    LWP_Deps = [
         'CDSC',    'CDSCB',    'CDSCD',
         'CIT',     'NFACTOR',  'VSAT',
         'A0',      'AGS',      'A1',       'A2',
         'KETA',    'NGATE',    'K1',       'K2',
         'VTH0',
         'UA',      'UB',       'UC',       'U0',
         'VOFF',    'DELTA',
         'RDSW',    'PRWG',     'PRWB',
         'ETA0',    'ETAB',
         'PCLM',    'PDIBLC1',  'PDIBLC2',
         'PSCBE1',  'PSCBE2',
         'WR',      'AT',       'KT1',     'KT2',
         'UTE',     'UA1',      'UB1',     'UC1',
         'PRT']
    L_Deps = [
         'CGSL',    'CGDL',
         'CKAPPA',  'CF',       'CLC',     'CLE',
         'VOFFCV',  'NOFF',     'ACDE',    'MOIN']

    def checkParamFix(self, var, check, fix):
        '''
        @param var      string, name of the param
        @param check    function that checks var
        @param fix      new value for the var if check fails
        '''

        v = self.__dict__[var]
        if check(v):
            return True
        else:
            if isinstance(v, ADVar):
                v.val = fix
            else:
                self.__dict__[var] = fix
            return False

    def __init__(self, **kwargs):
        kwArgs={}
        for k,v in kwargs.iteritems():
            kwArgs[string.upper(k)] = v

        for k,v in self.Param_Defs.iteritems():
            self.__dict__[k] = kwArgs.get(k, v)

        if not self.TOXM: self.TOXM = self.TOX
        if not self.DSUB: self.DSUB = self.DROUT
        self.TNOM += 273.15
        self.TEMP += 273.15

        # cap dW/dL params copied from DC params
        for k in ['LL', 'LW', 'LWL', 'WL', 'WW', 'WWL']:
            kk = k+'C'
            if not kwArgs.has_key(kk):
                self.__dict__[kk] = self.__dict__[k]

        fac_scaling    = sqrt(epsSi/epsOx * self.TOX)
        Vtm0           = kboq * self.TNOM              # nominal thermal voltage
        Vtm            = kboq * self.TEMP              # thermal voltage
        Eg0            = 1.16 - (7.02e-4 * self.TNOM * self.TNOM) / (self.TNOM + 1108.0)
        Eg             = 1.16 - (7.02e-4 * self.TEMP * self.TEMP) / (self.TEMP + 1108.0)
        ni             = 1.45e10 * (self.TNOM / 300.15) * sqrt(self.TNOM / 300.15) * exp(21.5565981 - Eg0 / (2.0 * Vtm0))

        ###########  1  ############
        #HSPICE XL/WL
        self.Ldrn = self.L + self.XL
        self.Wdrn = self.W + self.XW

        # effective length and width
        # see BSIM3v3 manual app-B-1.9
        t0 = Pow(self.L, self.LLN)
        t1 = Pow(self.W, self.LWN)
        tmp1 = self.LL / t0 + self.LW / t1 + self.LWL / (t0 * t1)
        dL = self.LINT + tmp1
        tmp2 = self.LLC / t0 + self.LWC / t1 + self.LWLC / (t0 * t1)
        #dLc = self.DLC + tmp2

        t2 = Pow(self.L, self.WLN)
        t3 = Pow(self.W, self.WWN)
        dW = self.WINT + self.WL / t2 + self.WW / t3 + self.WWL / (t2 * t3)
        tmp4 = self.WLC / t2 + self.WWC / t3 + self.WWLC / (t2 * t3)
        #dWc = self.DWC + tmp4

        Leff = self.Ldrn - 2.0 * dL                   # effective channel length
        Weff0 = self.Wdrn - 2.0 * dW                   # effective channel width

        #LeffCV = self.L - 2.0 * dLc
        ############################

        # {{{ LWP dependent parameters
        if self.BINUNIT==1 :
            invLeff = 1.0/Leff
            invWeff0 = 1.0/Weff0
            invLWeff0 = 1.0/(Leff * Weff0)
        else:
            invLeff = 1e-6/Leff
            invWeff0 = 1e-6/Weff0
            invLWeff0 = 1e-12/(Leff * Weff0)

        for k in self.LWP_Deps:
            if kwArgs.has_key('L'+k):
                ldep = kwArgs['L'+k]
                if not ldep==0.0:
                    self.__dict__[k] += ldep*invLeff
            if kwArgs.has_key('W'+k):
                wdep = kwArgs['W'+k]
                if not wdep==0.0:
                    self.__dict__[k] += wdep*invWeff0
            if kwArgs.has_key('P'+k):
                pdep = kwArgs['P'+k]
                if not pdep==0.0:
                    self.__dict__[k] += pdep*invLWeff0

        # }}}

        if self.U0>1.0: self.U0 *= 1e-4

        # {{{ S/D Capacitance area calculation
        if self.CALCACM==1:
            Weff_diode = self.W * self.WMLT + self.XW
            HDIFeff = self.HDIF * self.WMLT
            if kwArgs.has_key('AD'):
                self.AD *= self.WMLT*self.WMLT
            else:
                self.AD = 2. * HDIFeff * Weff_diode
            if kwArgs.has_key('AS'):
                self.AS *= self.WMLT*self.WMLT
            else:
                self.AS = 2. * HDIFeff * Weff_diode

            if kwArgs.has_key('PS'):
                self.PS *= self.WMLT
            else:
                self.PS = 4. * HDIFeff + 2. * Weff_diode
            if kwArgs.has_key('PD'):
                self.PD *= self.WMLT
            else:
                self.PD = 4. * HDIFeff + 2. * Weff_diode
        # }}}

        ###########  2  ############
        # saturation velocity
        t4 = (self.TEMP/self.TNOM - 1.0)                # TempRatio
        VsatTemp = self.VSAT - self.AT * t4
        Rds0 = (self.RDSW + self.PRT * t4) / Pow(Weff0 * 1e6, self.WR)
        ############################

        # gate capacitance
        Cox = epsOx / self.TOX

        ###########  3  ############
        # electrostatics
        if not kwArgs.has_key('NCH') and kwArgs.has_key('GAMMA1'):
            # see BSIM3v3 manual app-A, notes NI-4
            self.NCH = Pow(self.GAMMA1 * Cox, 2.) / ( 2.*q0*epsSi )

        if self.NCH <= ni:
            self.NCH = float(ni)*10.

        # bulk potential and its powers
        phi = 2.0 * Vtm0 * log(self.NCH / ni)
        sqrtPhi = sqrt(phi)

        # built-in potential, depletion length
        Xdep0   = sqrt(2.0 * epsSi / (q0 * self.NCH * 1.0e6)) * sqrtPhi # depletion width
        Cdep0   = epsSi / Xdep0
        litl    = sqrt(3.0 * self.XJ * self.TOX)
        Lt0     = fac_scaling * sqrt(Xdep0)

        Vbi     = Vtm0 * log(1.0e20 * self.NCH / (ni * ni))             # built-in potential

        t0 = exp(-0.5*self.DROUT*Leff/Lt0)
        thetaRout = self.PDIBLC1 * ( t0 + 2.*t0*t0 ) + self.PDIBLC2

        # see BSIM3v3 manual app-A, notes NI-7
        if not kwArgs.has_key('VBX'):
            self.VBX = phi - (q0 * self.NCH * self.XT * self.XT) / (2.*epsSi)

        # see BSIM3v3 manual app-A, notes NI-5,6
        if not kwArgs.has_key('GAMMA1'):
            gamma1 = 5.753e-12 * sqrt(self.NCH) / Cox                   # sqrt(2*q*epsSi*NCH)/Cox
        if not kwArgs.has_key('GAMMA2'):
            gamma2 = 5.753e-12 * sqrt(self.NSUB) / Cox                  # sqrt(2*q*epsSi*NSUB)/Cox


        if not (kwArgs.has_key('K1') and kwArgs.has_key('K2')):
            # if values of K1 and K2 are not supplied, calculate them
            # see BSIM3v3 manual app-A, notes NI-2
            t1 = sqrt(phi-self.VBX) - sqrtPhi
            t2 = 2. * sqrtPhi * (sqrt(phi-self.VBM) - sqrtPhi) + self.VBM
            K2 = (gamma1-gamma2) * t1 / t2
            K1 = gamma2 - 2*K2* sqrt( phi - self.VBM )
        else:
            K1 = self.K1
            K2 = self.K2

        # see BSIM3v3 manual app-A, notes NI-1
        if kwArgs.has_key('VTH0'):
            if not kwArgs.has_key('VFB'):
                self.VFB = self.VTH0 - phi - K1 * sqrtPhi
        else:
            self.VTH0 = self.VFB + phi + K1 * sqrtPhi

        # diode current
        t0 = Eg0/Vtm0 - Eg/Vtm + self.XTI * log(self.TEMP/self.TNOM)
        t1 = exp(t0/self.NJ)
        Js = self.JS * t1
        Jssw = self.JSSW * t1
        Isbs = self.AS * Js + self.PS * Jssw
        Isbd = self.AD * Js + self.PD * Jssw

        #########
        self.fac_scaling    = fac_scaling
        self.Leff           = Leff
        self.Weff0          = Weff0
        self.Cox            = Cox
        self.phi            = phi
        self.Vtm0           = Vtm0
        self.Vtm            = Vtm
        self.Eg0            = Eg0
        self.ni             = ni
        self.Vbi            = Vbi
        self.K1             = K1
        self.K2             = K2
        self.Xdep0          = Xdep0
        self.Lt0            = Lt0
        self.thetaRout      = thetaRout
        self.VsatTemp       = VsatTemp
        self.Rds0           = Rds0
        self.litl           = litl
        self.Cdep0          = Cdep0
        self.Isbs           = Isbs
        self.Isbd           = Isbd


        self.checkParamFix('DROUT', lambda x: x>=0., 0.0)
        self.checkParamFix('DVT1',  lambda x: x>=0., 0.0)
        self.checkParamFix('DVT1W', lambda x: x>=0., 0.0)
        self.checkParamFix('PDIBLC1',   lambda x: x>=0., 0.0)
        self.checkParamFix('PDIBLC2',   lambda x: x>=0., 0.0)

    def _DC_Curr(self, VGS, VDS, VBS):
        Leff                = self.Leff
        Weff0               = self.Weff0
        fac_scaling         = self.fac_scaling
        Cox                 = self.Cox
        phi                 = self.phi
        sqrtPhi             = sqrt(phi)
        phis3               = sqrtPhi * phi;
        Vtm0                = self.Vtm0
        Vtm                 = self.Vtm
        Eg0                 = self.Eg0
        ni                  = self.ni
        Vbi                 = self.Vbi
        K1                  = self.K1
        K2                  = self.K2
        Xdep0               = self.Xdep0
        VsatTemp            = self.VsatTemp
        Lt0                 = self.Lt0
        thetaRout           = self.thetaRout
        litl                = self.litl
        Cdep0               = self.Cdep0
        Isbs                = self.Isbs
        Isbd                = self.Isbd

        if K2<0.0:
            Vbc = 0.9 * ( self.phi - Pow(0.5*K1/K2, 2.) )
            if Vbc > -3.0:
                Vbc = -3.0
            elif Vbc < -30.0:
                Vbc = -30.0
        else:
            Vbc = -30.0
        if Vbc > self.VBM: Vbc = self.VBM

        # effective Vbs
        t0 = VBS - Vbc - 0.001
        t1 = sqrt(t0 * t0 - 0.004 * Vbc)
        Vbseff = Vbc + 0.5 * (t0 + t1)

        # Calculate Phis, sqrtPhis and Xdep
        if Vbseff>0:
            Phis = phi*phi/(phi+Vbseff)
            sqrtPhis = phis3 / (phi+0.5*Vbseff)
        else:
            Phis = phi - Vbseff
            sqrtPhis = sqrt(Phis)
        Xdep = Xdep0 * sqrtPhis / sqrtPhi

        # {{{ Calculation of Threshold voltage
        ######################
        t3 = sqrt(Xdep)
        V0 = Vbi - phi
        if self.DVT2*Vbseff + 0.5 > 0:
            t1 = 1.0 + self.DVT2*Vbseff
        else:
            t1 = (1.0 + 3.0*self.DVT2*Vbseff) / (3.0 + 8.0*self.DVT2*Vbseff)
        Ltl = fac_scaling * t3 * t1   # char. length of potential along channel

        if self.DVT2W*Vbseff + 0.5 > 0:
            t1 = 1.0 + self.DVT2W*Vbseff
        else:
            t1 = (1.0 + 3.0*self.DVT2W*Vbseff) / (3.0 + 8.0*self.DVT2W*Vbseff)
        Ltw = fac_scaling * t3 * t1 # char. length of potential across channel

        # length dependence of Vth roll-off
        t2 = exp(-0.5 * self.DVT1 * Leff / Ltl)
        t2 = t2 * ( 1. + 2. * t2)
        dVth_L = self.DVT0 * t2 * V0

        # width dependence of Vth roll-off
        t2 = exp(-0.5 * self.DVT1W * Weff0 * Leff / Ltw)
        t2 = t2 * ( 1. + 2. * t2)
        dVth_W = self.DVT0W * t2 * V0

        # DIBL Vth shift
        t2 = exp(-0.5 * self.DSUB * Leff / Lt0)
        t2 = t2 * ( 1. + 2. * t2)
        dVth_DIBL = VDS * t2 * ( self.ETA0 + self.ETAB * Vbseff )

        K1ox = self.K1 * self.TOX / self.TOXM
        K2ox = self.K2 * self.TOX / self.TOXM

        # lateral non-uniform doping
        dVth_lat = K1ox * (sqrt(1.+self.NLX/Leff) - 1.) * sqrtPhi

        # temperature effect
        dVth_temp = (self.KT1 + self.KT1L/Leff + self.KT2*Vbseff) * (self.TEMP/self.TNOM - 1.0)

        # narrow width effect
        dVth_nrw = (self.K3 + self.K3B * Vbseff) * self.TOX / (Weff0 + self.W0) * phi

        Vth = self.VTH0 - self.K1*sqrtPhi + K1ox*sqrtPhis - K2ox * Vbseff \
              + dVth_lat + dVth_nrw \
              - dVth_L - dVth_W - dVth_DIBL + dVth_temp
        
        ######################
        # }}}


        # {{{ Poly Depletion Effect
        t0 = VGS - (self.VFB + phi)
        if t0>0 and self.NGATE>1e18 and self.NGATE<1e25:
            t1 = 1e6 * q0 * epsSi * self.NGATE / ( Cox*Cox )
            t2 = t1 * ( sqrt(1. + 2. * t0 / t1) - 1.0 )
            Vgseff = self.VFB+phi+t2
        else:
            Vgseff = VGS
            
        # }}}

        # {{{ Calculatation of Vgsteff, the effective Vgs-Vth

        # factor n
        t0 = self.NFACTOR * epsSi / Xdep
        t1 = self.CDSC + self.CDSCB * Vbseff + self.CDSCD * VDS
        t2 = exp(-0.5 * self.DVT1 * Leff / Ltl)
        t2 = t2 * ( 1. + 2. * t2)
        t3 = (t0 + self.CIT + t1*t2) / Cox
        if t3 > -0.5:
            fac_n = 1. + t3
        else:
            fac_n = (1. + 3.*t3) * ( 1. / ( 3. + 8.*t3 ) )

        t0 = 1.0 / fac_n /Vtm
        t1 = 0.5 * (Vgseff - Vth) * t0
        t2 = exp(self.VOFF*t0 - t1)
        Vgsteff = 2.*fac_n * Vtm * log1pexp(t1) / ( 1. + 2.*fac_n * Cox / Cdep0 *t2 )
        # }}}


        # voltage dependant effective width
        Weff = Weff0 - 2. * (self.DWG*Vgsteff + self.DWB * (sqrtPhis - sqrtPhi))

        # voltage dependant source/drain resistance
        t0 = self.PRWG*Vgsteff + self.PRWB*(sqrtPhis-sqrtPhi)
        if t0>-0.9:
            Rds = self.Rds0 * (1.+t0)
        else:
            Rds = self.Rds0 * (0.8+t0) * (1.0 / (17.0 + 20.0 * t0))


        # {{{ Calculation of Abulk
        t1 = 0.5 * K1ox / sqrtPhis
        t2 = Leff / (Leff + 2.0 * sqrt(self.XJ * Xdep))
        t3 = (self.A0 * t2) + (self.B0 / (Weff0 + self.B1))
        t4 = self.AGS * self.A0 * t2*t2*t2

        Abulk0 = 1.0 + t1 * t3
        Abulk = Abulk0 + (-t1 * t4) * Vgsteff

        if Abulk0<0.1:
            Abulk0 = (0.2-Abulk0) / (3.0-20.0*Abulk0)
        if Abulk <0.1:
            Abulk = (0.2-Abulk) / (3.0-20.0*Abulk)

        t2 = self.KETA * Vbseff
        if t2>-0.9:
            t0 = 1. / (1. + t2)
        else:
            t0 = (17. + 20. * t2) / (0.8 + t2)

        Abulk *= t0
        Abulk0 *= t0
        # }}}

        # {{{ Mobility calculation (mobMOD=1)
        # temperature correction
        t0 = self.TEMP/self.TNOM - 1.0
        Ua = self.UA + self.UA1 * t0
        Ub = self.UB + self.UB1 * t0
        Uc = self.UC + self.UC1 * t0
        U0 = self.U0 * Pow(self.TEMP/self.TNOM, self.UTE)

        t2 = Ua + Uc * Vbseff
        t3 = (Vgsteff + Vth + Vth) / self.TOX
        t5 = t3 * (t2 + Ub * t3)
        if t5>-0.8:
            denomi = 1.0 + t5
        else:
            denomi = (0.6 + t5) / (7.0 + 10.0 * t5)
        Ueff = U0 / denomi
        # }}}

        # Saturation E-field
        Esat = 2.0 * VsatTemp / Ueff

        # Saturation Drain Voltage Vdsat
        lmd = self.A1 * Vgsteff + self.A2

        Vgst2Vtm = Vgsteff + 2.0 * Vtm;
        WVCoxRds = Weff * VsatTemp * Cox * Rds

        if lmd==1.0 and Rds==0.0:
            Vdsat = Esat*Leff*Vgst2Vtm / ( Abulk*Esat*Leff + Vgst2Vtm )
        else:
            t1 = Abulk*Abulk * WVCoxRds + (1./lmd -1.)*Abulk # a
            t2 = - ( Vgst2Vtm*(2./lmd-1.) + Abulk*Esat*Leff + 3.*Abulk*Vgst2Vtm*WVCoxRds ) # b
            t3 = Vgst2Vtm*Esat*Leff + 2. *Vgst2Vtm*Vgst2Vtm * WVCoxRds # c
            Vdsat = ( -t2 - sqrt(t2*t2 - 4.*t1*t3) ) / (2.*t1)


        # effective Vds
        t0 = Vdsat - VDS - self.DELTA
        Vdseff = Vdsat - 0.5*(t0 + sqrt(t0*t0 + 4.*self.DELTA*Vdsat) )
        if VDS==0.0:
            Vdseff=0.0

        diffVds = VDS - Vdseff

        # VAsat
        t0 = 1.0 - 0.5 * Abulk*Vdsat / Vgst2Vtm
        t1 = Esat*Leff + Vdsat + 2. * WVCoxRds * Vgsteff * t0
        t2 = 2./lmd - 1. + WVCoxRds*Abulk
        Vasat = t1/t2

        # VA_CLM
        if self.PCLM>0.0 and diffVds>1.0e-10:
            VaCLM = (Abulk*Esat*Leff + Vgsteff) / ( self.PCLM*Abulk*Esat*litl ) * diffVds
        else:
            VaCLM = 1e8

        # VA_DIBLC
        if thetaRout>0.0:
            t1 = Abulk*Vdsat
            VaDIBLC = Vgst2Vtm/(thetaRout * (1.+self.PDIBLCB*Vbseff)) * (1.-t1/(t1+Vgst2Vtm))
        else:
            VaDIBLC = 1e8

        # VA
        t0 = self.PVAG * Vgsteff / (Esat*Leff)
        if t0>-0.9:
            t0 = 1. + t0
        else:
            t0 = (0.8+t0) * (1. / (17. + 20.*t0) )
        Va = Vasat + t0 * (VaCLM*VaDIBLC)/(VaCLM + VaDIBLC)


        # one over VaSCBE
        #if 100.*diffVds > self.PSCBE1*litl:
        if diffVds > 1e-10:
            rcpVaSCBE = self.PSCBE2/Leff * exp(-self.PSCBE1*litl/diffVds)
        else:
            rcpVaSCBE = 0.0

        # Calculation of Ids
        CoxWovL = Cox * Weff / Leff
        beta = Ueff * CoxWovL
        fgche1 = Vgsteff * (1.0 - 0.5 * Abulk * Vdseff / Vgst2Vtm)
        fgche2 = 1.0 + Vdseff / (Esat*Leff)
        gche = beta * fgche1 / fgche2

        t0 = Vdseff / (1.0 + gche * Rds)
        Idl = gche * t0

        Idsa = Idl * (1.0 + diffVds / Va)
        Ids = Idsa * (1.0 + diffVds * rcpVaSCBE)

        # Substrate current
        if diffVds>1e-5:
            t0 = exp(-self.BETA0/diffVds)
        else:
            t0 = 0.0
        Isub = (self.ALPHA0 + self.ALPHA1*Leff)/Leff * diffVds * t0 * Idsa

        # diode current
        NVtm = self.NJ * Vtm
        if Isbs>0.0:
            if self.IJTH==0.0:
                Ibs = Isbs*(exp(VBS/NVtm)-1.)
            else:
                Vjsm = NVtm * log(self.IJTH/Isbs+1.)
                if VBS<Vjsm:
                    Ibs = Isbs*(exp(VBS/NVtm)-1.)
                else:
                    Ibs = self.IJTH + (self.IJTH+Isbs)/NVtm * (VBS-Vjsm)
        else:
            Ibs = 0.0

        if Isbd>0.0:
            VBD = VBS-VDS
            if self.IJTH==0.0:
                Ibd = Isbd*(exp(VBD/NVtm)-1.)
            else:
                Vjdm = NVtm * log(self.IJTH/Isbd+1.)
                if VBD<Vjdm:
                    Ibd = Isbd*(exp(VBD/NVtm)-1.)
                else:
                    Ibd = self.IJTH + (self.IJTH+Isbd)/NVtm * (VBD-Vjdm)
        else:
            Ibd = 0.0

        #res = [Ids+Isub-Ibd, -Ids-Ibs, Isub+Ibs+Ibd]  # Id, Is, Isub
        res = []
        for v in output:
            if v=='Id':
                res.append(Ids+Isub-Ibd)
            elif v=='Is':
                res.append(-Ids-Ibs)
            elif v=='Isub':
                res.append(Isub+Ibs+Ibd)
            elif v=='Vth':
                res.append(Vth)
            elif v=='fac_n':
                res.append(fac_n)
            elif v=='Vgsteff':
                res.append(Vgsteff)
            elif v=='Cox':
                res.append(Cox)
        return tuple(res)

    def calcFunJac(self, state):
        VD, VG, VS, VB = state.getVars(self.varIdx)
        VDS = VD - VS
        VGS = VG - VS
        VBS = VB - VS

        Id, Is, Isub = self._Ids(VGS, VDS, VBS)

if __name__=='__main__':
    mos = MOSBSim3v3(
        L = 0.19e-6, W=10e-6, TEMP=25.,
        TOX      = 3.87E-09,            TOXM     = 3.87E-09,            XJ       = 1.6000000E-07,
        NCH      = 3.8694000E+17,       LLN      = 1.1205959,           LWN      = 0.9200000,
        WLN      = 1.0599999,           WWN      = 0.8768474,           LINT     = 1.5757085E-08,
        LL       = 2.6352781E-16,       LW       = -2.2625584E-16,      LWL      = -2.0576711E-22,
        WINT     = -1.4450482E-09,      WL       = -2.3664573E-16,      WW       = -3.6409690E-14,
        WWL      = -4.0000000E-21,      MOBMOD   = 1,
        XL       = 1.8E-8,              XW       = 0.00,                DWG      = -5.9600000E-09,
        DWB      = 4.5000000E-09,
        ##### DIODE PARAMETERS
        CALCACM  = 1,
        ACM      = 12,                  HDIF     = 2.00E-07,
        RSH      = 7.08,                RD       = 0,                   RS       = 0,
        RSC      = 1.7,                 RDC      = 1.7,
        JS       = 3.52E-07,            JSW      = 3.0E-13,             NJ       = 1.0392,
        XTI      = 3.25, 
        ##### Vth parameter
        VTH0     = 0.39,                WVTH0    = -2.9709472E-08,      PVTH0    = 5.0000000E-16,
        K1       = 0.6801043,           WK1      = -2.4896840E-08,      PK1      = 1.3000000E-15,
        K2       = -4.9977830E-02,      K3       = 10.0000000,          DVT0     = 1.3000000,
        DVT1     = 0.5771635,           DVT2     = -0.1717554,          DVT0W    = 0.00,
        DVT1W    = 0.00,                DVT2W    = 0.00,                NLX      = 7.5451030E-08,
        W0       = 5.5820150E-07,       K3B      = -3.0000000,
        ##### Mobility
        VSAT     = 8.2500000E+04,       PVSAT    = -8.3000000E-10,      UA       = -1.0300000E-9,
        LUA      = 7.7349790E-19,       PUA      = -1.0000000E-24,      UB       = 2.3666682E-18,
        UC       = 1.2000000E-10,       PUC      = 1.5000000E-24,       RDSW     = 55.5497200,
        PRWB     = -0.2400000,          PRWG     = 0.4000000,           WR       = 1.0000000,
        U0       = 3.4000000E-02,       LU0      = 2.3057663E-11,       WU0      = -3.1009695E-09,
        A0       = 0.8300000,           KETA     = -3.0000000E-03,      LKETA    = -1.7000000E-09,
        A1       = 0.00,                A2       = 0.9900000,           AGS      = 0.3200000,
        B0       = 6.0000000E-08,       B1       = 0.00,
        ##### Subthreshold
        VOFF     = -0.1030000,          LVOFF    = -3.3000000E-09,      NFACTOR  = 1.2500000,
        LNFACTOR = 4.5000000E-08,       CIT      = 0.00,                CDSC     = 0.00,
        CDSCB    = 0.00,                CDSCD    = 1.0000000E-04,       ETA0     = 2.8000001E-02,
        ETAB     = -2.7000001E-02,      DSUB     = 0.4000000,
        ###### ROUT PARAMETERS
        PCLM     = 1.2000000,           PPCLM    = 2.9999999E-15,       PDIBLC1  = 2.5000000E-02,
        PDIBLC2  = 3.8000000E-03,       PPDIBLC2 = 2.7000001E-16,       PDIBLCB  = 0.00,
        DROUT    = 0.5600000,           PSCBE1   = 3.4500000E+08,       PSCBE2   = 1.0000000E-06,
        PVAG     = 0.00,                DELTA    = 1.0000000E-02,       ALPHA0   = 1.7753978E-08,
        ALPHA1   = 0.1764000,           LALPHA1  = 7.6250000E-09,       BETA0    = 11.1683940,
        ##### TEMPERATURE EFFECTS PARAMETERS
        KT1      = -0.2572866,          KT2      = -4.0000000E-02,      AT       = 3.7000000E+04,
        PAT      = -7.5000000E-10,      UTE      = -1.5500000,          UA1      = 1.7600000E-09,
        LUA1     = 6.0000000E-18,       WUA1     = -1.1000000E-16,      PUA1     = -5.0000000E-25,
        UB1      = -2.4000000E-18,      UC1      = -1.0000000E-10,      LUC1     = 1.6999999E-17,
        PUC1     = -3.0000000E-24,      KT1L     = -1.0000000E-09,      PRT      = -55.0000000,

            )

    for Vbs in [0.]:
        for j in xrange(1):
            Vds = 0.05
            for k in xrange(41):
                Vgs = 0.05*k
                Id, Is, Isub = mos._DC_Curr(Vgs,Vds,Vbs)
                print "%15g %15g %15g %15g %15g" % (Vgs, Vds, Vbs, Id, Isub)
                Vth, fac_n = mos._DC_Curr(Vgs,Vds,Vbs,output=('Vth','fac_n'))
                print '******', Vth, fac_n

#    for j in xrange(37):
#        Vbs = 0.
#        Vds = j*0.05
#        Vgs = 1.8
#        Id, Is, Isub = mos._DC_Curr(Vgs,Vds,Vbs)
#        print Vgs, Vds, Vbs, Id

#-1.20150e-09
