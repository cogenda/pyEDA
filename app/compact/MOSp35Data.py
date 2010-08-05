from pyEDA.Compact.MOSParamFit import *
from pyEDA.Compact.AuroraData import *

class MOSp35Data(MOS_FitProject):
    def __init__(self, param0):
        super(MOSp35Data, self).__init__(param0)

        print 'loading data files...'
        self.loadAuroraFile('data/level49.dat')
        print self.datasets.keys()
        print 'done.'

        # {{{
        # Idvg, long/wide, linear region, zero Vb
        dsrc = MOS_IV_FitData('Idvg, long/wide, linear region, zero Vb')
        ds = self.datasets['n15x15.gtr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':0.,    'Vds':0.1}, 'Id'))
        self.IdVg_LW_lin_b0 = dsrc

        # Idvg, mid/wide, linear region, zero Vb
        dsrc = MOS_IV_FitData('Idvg, mid/wide, linear region, zero Vb')
        ds = self.datasets['n15x1p8.gtr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':0.,    'Vds':0.1}, 'Id'))
        self.IdVg_MW_lin_b0 = dsrc

        # Idvg, short/wide, linear region, zero Vb
        dsrc = MOS_IV_FitData('Idvg, short/wide, linear region, zero Vb')
        ds = self.datasets['n15x0p6.gtr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':0.,    'Vds':0.1}, 'Id'))
        self.IdVg_SW_lin_b0 = dsrc

        # Idvg, long/narrow, linear region, zero Vb
        dsrc = MOS_IV_FitData('Idvg, long/narrow, linear region, zero Vb')
        ds = self.datasets['n1p8x15.gtr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':0.,    'Vds':0.1}, 'Id'))
        self.IdVg_LN_lin_b0 = dsrc

        # Idvg, short/narrow, linear region, zero Vb
        dsrc = MOS_IV_FitData('Idvg, short/narrow, linear region, zero Vb')
        ds = self.datasets['n1p8x0p6.gtr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':0.,    'Vds':0.1}, 'Id'))
        self.IdVg_SN_lin_b0 = dsrc

        # ---
        # Idvg, long/wide, saturation region, zero Vb
        dsrc = MOS_IV_FitData('Idvg, long/wide, saturation region, zero Vb')
        ds = self.datasets['n15x15.gtr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':0.,    'Vds':3.3}, 'Id'))
        self.IdVg_LW_sat_b0 = dsrc

        # Idvg, mid/wide, saturation region, zero Vb
        dsrc = MOS_IV_FitData('Idvg, mid/wide, saturation region, zero Vb')
        ds = self.datasets['n15x1p8.gtr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':0.,    'Vds':3.3}, 'Id'))
        self.IdVg_MW_sat_b0 = dsrc

        # Idvg, short/wide, saturation region, zero Vb
        dsrc = MOS_IV_FitData('Idvg, short/wide, saturation region, zero Vb')
        ds = self.datasets['n15x0p6.gtr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':0.,    'Vds':3.3}, 'Id'))
        self.IdVg_SW_sat_b0 = dsrc

        # Idvg, long/narrow, saturation region, zero Vb
        dsrc = MOS_IV_FitData('Idvg, long/narrow, saturation region, zero Vb')
        ds = self.datasets['n1p8x15.gtr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':0.,    'Vds':3.3}, 'Id'))
        self.IdVg_LN_sat_b0 = dsrc

        # Idvg, short/narrow, saturation region, zero Vb
        dsrc = MOS_IV_FitData('Idvg, short/narrow, saturation region, zero Vb')
        ds = self.datasets['n1p8x0p6.gtr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':0.,    'Vds':3.3}, 'Id'))
        self.IdVg_SN_sat_b0 = dsrc

        # ---
        # Idvg, long/wide, Vd=1, zero Vb
        dsrc = MOS_IV_FitData('Idvg, long/wide, Vd=1, zero Vb')
        ds = self.datasets['n15x15.gtr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':0.,    'Vds':1.0}, 'Id'))
        self.IdVg_LW_d1_b0 = dsrc

        # Idvg, mid/wide, Vd=1, zero Vb
        dsrc = MOS_IV_FitData('Idvg, mid/wide, Vd=1, zero Vb')
        ds = self.datasets['n15x1p8.gtr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':0.,    'Vds':1.0}, 'Id'))
        self.IdVg_MW_d1_b0 = dsrc

        # Idvg, short/wide, Vd=1, zero Vb
        dsrc = MOS_IV_FitData('Idvg, short/wide, Vd=1, zero Vb')
        ds = self.datasets['n15x0p6.gtr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':0.,    'Vds':1.0}, 'Id'))
        self.IdVg_SW_d1_b0 = dsrc

        # Idvg, long/narrow, Vd=1, zero Vb
        dsrc = MOS_IV_FitData('Idvg, long/narrow, Vd=1, zero Vb')
        ds = self.datasets['n1p8x15.gtr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':0.,    'Vds':1.0}, 'Id'))
        self.IdVg_LN_d1_b0 = dsrc

        # Idvg, short/narrow, saturation region, zero Vb
        dsrc = MOS_IV_FitData('Idvg, short/narrow, saturation region, zero Vb')
        ds = self.datasets['n1p8x0p6.gtr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':0.,    'Vds':1.0}, 'Id'))
        self.IdVg_SN_d1_b0 = dsrc

        # ---
        # IdVg, long/wide, linear region, all Vb
        dsrc = MOS_IV_FitData('IdVg, long/wide, linear region, all Vb')
        ds = self.datasets['n15x15.gtr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':0.,     'Vds':0.1}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-1.1,   'Vds':0.1}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-2.2,   'Vds':0.1}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-3.3,   'Vds':0.1}, 'Id'))
        self.IdVg_LW_lin_ba = dsrc

        # IdVg, long/narrow, linear region, all Vb
        dsrc = MOS_IV_FitData('IdVg, long/narrow, linear region, all Vb')
        ds = self.datasets['n1p8x15.gtr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':0.,     'Vds':0.1}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-1.1,   'Vds':0.1}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-2.2,   'Vds':0.1}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-3.3,   'Vds':0.1}, 'Id'))
        self.IdVg_LN_lin_ba = dsrc

        # IdVg, mid/wide, linear region, all Vb
        dsrc = MOS_IV_FitData('IdVg, mid/wide, linear region, all Vb')
        ds = self.datasets['n15x1p8.gtr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':0.,     'Vds':0.1}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-1.1,   'Vds':0.1}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-2.2,   'Vds':0.1}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-3.3,   'Vds':0.1}, 'Id'))
        self.IdVg_MW_lin_ba = dsrc

        # IdVg, short/wide, linear region, all Vb
        dsrc = MOS_IV_FitData('IdVg, short/wide, linear region, all Vb')
        ds = self.datasets['n15x0p6.gtr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':0.,     'Vds':0.1}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-1.1,   'Vds':0.1}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-2.2,   'Vds':0.1}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-3.3,   'Vds':0.1}, 'Id'))
        self.IdVg_SW_lin_ba = dsrc

        # IdVg, short/narrow, linear region, all Vb
        dsrc = MOS_IV_FitData('IdVg, short/narrow, linear region, all Vb')
        ds = self.datasets['n1p8x0p6.gtr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':0.,     'Vds':0.1}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-1.1,   'Vds':0.1}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-2.2,   'Vds':0.1}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-3.3,   'Vds':0.1}, 'Id'))
        self.IdVg_SN_lin_ba = dsrc

        #---
        # IdVg, long/wide , Vd=1, all Vb
        dsrc = MOS_IV_FitData('IdVg, long/wide , Vd=1, all Vb')
        ds = self.datasets['n15x15.gtr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':0.,     'Vds':1.0}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-1.1,   'Vds':1.0}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-2.2,   'Vds':1.0}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-3.3,   'Vds':1.0}, 'Id'))
        self.IdVg_LW_d1_ba = dsrc

        # IdVg, mid/wide, Vd=1, all Vb
        dsrc = MOS_IV_FitData('IdVg, mid/wide, Vd=1, all Vb')
        ds = self.datasets['n15x1p8.gtr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':0.,     'Vds':1.0}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-1.1,   'Vds':1.0}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-2.2,   'Vds':1.0}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-3.3,   'Vds':1.0}, 'Id'))
        self.IdVg_MW_d1_ba = dsrc

        # IdVg, short/wide, Vd=1, all Vb
        dsrc = MOS_IV_FitData('IdVg, short/wide, Vd=1, all Vb')
        ds = self.datasets['n15x0p6.gtr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':0.,     'Vds':1.0}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-1.1,   'Vds':1.0}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-2.2,   'Vds':1.0}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-3.3,   'Vds':1.0}, 'Id'))
        self.IdVg_SW_d1_ba = dsrc

        #---
        # IdVg, long/wide , saturation region, all Vb
        dsrc = MOS_IV_FitData('IdVg, long/wide , saturation region, all Vb')
        ds = self.datasets['n15x15.gtr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':0.,     'Vds':3.3}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-1.1,   'Vds':3.3}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-2.2,   'Vds':3.3}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-3.3,   'Vds':3.3}, 'Id'))
        self.IdVg_LW_sat_ba = dsrc

        # IdVg, mid/wide, saturation region, all Vb
        dsrc = MOS_IV_FitData('IdVg, mid/wide, saturation region, all Vb')
        ds = self.datasets['n15x1p8.gtr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':0.,     'Vds':3.3}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-1.1,   'Vds':3.3}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-2.2,   'Vds':3.3}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-3.3,   'Vds':3.3}, 'Id'))
        self.IdVg_MW_sat_ba = dsrc

        # IdVg, short/wide, saturation region, all Vb
        dsrc = MOS_IV_FitData('IdVg, short/wide, saturation region, all Vb')
        ds = self.datasets['n15x0p6.gtr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':0.,     'Vds':3.3}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-1.1,   'Vds':3.3}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-2.2,   'Vds':3.3}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vgs', {'Vbs':-3.3,   'Vds':3.3}, 'Id'))
        self.IdVg_SW_sat_ba = dsrc

        #---
        # IdVd, long/wide , zero Vb
        dsrc = MOS_IV_FitData('IdVd, long/wide , zero Vb')
        ds = self.datasets['n15x15.drr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vds', {'Vbs':0.,     'Vgs':0.9}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vds', {'Vbs':0.,     'Vgs':1.5}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vds', {'Vbs':0.,     'Vgs':2.1}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vds', {'Vbs':0.,     'Vgs':2.7}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vds', {'Vbs':0.,     'Vgs':3.3}, 'Id'))
        self.IdVd_LW_b0 = dsrc

        # IdVd, mid/wide , zero Vb
        dsrc = MOS_IV_FitData('IdVd, mid/wide , zero Vb')
        ds = self.datasets['n15x1p8.drr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vds', {'Vbs':0.,     'Vgs':0.9}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vds', {'Vbs':0.,     'Vgs':1.5}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vds', {'Vbs':0.,     'Vgs':2.1}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vds', {'Vbs':0.,     'Vgs':2.7}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vds', {'Vbs':0.,     'Vgs':3.3}, 'Id'))
        self.IdVd_MW_b0 = dsrc

        # IdVd, short/wide , zero Vb
        dsrc = MOS_IV_FitData('IdVd, short/wide , zero Vb')
        ds = self.datasets['n15x0p6.drr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vds', {'Vbs':0.,     'Vgs':0.9}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vds', {'Vbs':0.,     'Vgs':1.5}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vds', {'Vbs':0.,     'Vgs':2.1}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vds', {'Vbs':0.,     'Vgs':2.7}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vds', {'Vbs':0.,     'Vgs':3.3}, 'Id'))
        self.IdVd_SW_b0 = dsrc

        # IdVd, long/narrow, zero Vb
        dsrc = MOS_IV_FitData('IdVd, long/narrow, zero Vb')
        ds = self.datasets['n1p8x15.drr']
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vds', {'Vbs':0.,     'Vgs':0.9}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vds', {'Vbs':0.,     'Vgs':1.5}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vds', {'Vbs':0.,     'Vgs':2.1}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vds', {'Vbs':0.,     'Vgs':2.7}, 'Id'))
        dsrc.addCurve(ds.mosID(), ds.getCurve('Vds', {'Vbs':0.,     'Vgs':3.3}, 'Id'))
        self.IdVd_LN_b0 = dsrc

        # }}}

