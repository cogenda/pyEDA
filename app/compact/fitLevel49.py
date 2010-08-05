from pyEDA.Compact.MOSParamFit import *
from pyEDA.Compact.AuroraData import *
import math
import numpy as np
from matplotlib import pyplot

from BSIM3v3Param import *
from MOSp35Data import *

class MOSp35Proj(MOSp35Data):
    def __init__(self, param0):
        super(MOSp35Proj, self).__init__(param0)

    def step10(self):
        targets=['VTH0', 'U0', 'UA', 'UB']
        fit = MOS_IV_Fit(self.param, targets, 'Step 1')
        fit.setDataSource(self.IdVg_LW_lin_b0.copy())

        result,err = fit.doFit()
        fit.visualize(result, timeout=0.0)

        self.acceptParam(result, targets)

        return result, err


    def step20(self):
        targets=['K1', 'K2', 'UC']
        fit = MOS_IV_Fit(self.param, targets, 'Step 2')
        fit.setDataSource(self.IdVg_LW_lin_ba.copy())

        result,err = fit.doFit()
        fit.visualize(result, timeout=0.0)

        self.acceptParam(result, targets)

        return result, err


    def step30(self):
        targets = ['NFACTOR', 'VOFF']
        fit = MOS_IV_Fit(self.param, targets, 'Step 3')
        fit.setDataSource(1.0*self.IdVg_LW_lin_ba)

        result,err = fit.doFit()
        self.acceptParam(result, targets)

        return result, err

    def step35(self):
        targets = ['NFACTOR', 'VOFF']
        fit = MOS_IV_Fit(self.param, targets, 'Step 3-a')
        fit.setDataSource(self.IdVg_LW_lin_b0.copy())

        fit.dataSrc.setSubVth(True)
        result,err = fit.doFit()
        fit.visualize(result, timeout=0.0)
        self.acceptParam(result, targets)

        return result, err


    def step40(self):
        targets = ['K3', 'W0', 'WINT', 'DWG']
        fit = MOS_IV_Fit(self.param, targets, 'Step 4')
        fit.setDataSource(self.IdVg_LW_lin_b0 + 10.0 * self.IdVg_LN_lin_b0)

        result,err = fit.doFit()
        fit.visualize(result, timeout=0.0)

        self.acceptParam(result, targets)

        return result, err

    def step50(self):
        targets = ['K3B', 'DWB']
        fit = MOS_IV_Fit(self.param, targets, 'Step 5')
        fit.setDataSource(self.IdVg_LW_lin_ba + 15.0 * self.IdVg_LN_lin_ba)

        result,err = fit.doFit()
        fit.visualize(result, timeout=0.0)
        self.acceptParam(result, targets)

        return result, err

    def step70(self):
        targets = ['LINT', 'RDSW', 'DVT0', 'DVT1', 'DVT2', 'NLX']
        fit = MOS_IV_Fit(self.param, targets, 'Step 7')
        fit.setDataSource(self.IdVg_LW_lin_ba + 0.1*self.IdVg_MW_lin_ba + 0.05*self.IdVg_SW_lin_ba)

        result,err = fit.doFit()
        fit.visualize(result, timeout=0.0)
        self.acceptParam(result, targets)

        return result, err

    def step80(self):
        targets = ['LINT', 'RDSW', 'PRWG']
        fit = MOS_IV_Fit(self.param, targets, 'Step 8')
        fit.setDataSource(self.IdVg_LW_lin_b0 + 0.1*self.IdVg_MW_lin_b0 + 0.05*self.IdVg_SW_lin_b0)

        result,err = fit.doFit()
        fit.visualize(result, timeout=0.0)
        self.acceptParam(result, targets)

        return result, err

    def step90(self):
        targets = ['PRWB']
        fit = MOS_IV_Fit(self.param, targets, 'Step 9')
        fit.setDataSource(self.IdVg_LW_lin_ba + 0.1*self.IdVg_MW_lin_ba + 0.05*self.IdVg_SW_lin_ba)

        result,err = fit.doFit()
        fit.visualize(result, timeout=0.0)
        self.acceptParam(result, targets)

        return result, err

    def disable_step100(self):
        targets = ['DVT0W', 'DVT1W', 'RDSW', 'WR', 'PRWG']
        fit = MOS_IV_Fit(self.param, targets, 'Step 10')
        fit.setDataSource(     self.IdVg_LW_lin_ba + 
                          0.1 *self.IdVg_MW_lin_ba +
                          0.05*self.IdVg_SW_lin_ba +
                          10. *self.IdVg_LN_lin_ba +
                          2.0*self.IdVg_SN_lin_ba)

        result,err = fit.doFit()
        fit.visualize(result, timeout=0.0)
        self.acceptParam(result, targets)
        return result, err

    def step110(self):
        targets = ['RDSW', 'WR', 'PRWG']
        fit = MOS_IV_Fit(self.param, targets, 'Step 11')
        fit.setDataSource(     self.IdVg_LW_lin_b0 + 
                          0.1 *self.IdVg_MW_lin_b0 +
                          0.05*self.IdVg_SW_lin_b0 +
                          10. *self.IdVg_LN_lin_b0 +
                          2.0*self.IdVg_SN_lin_b0)

        result,err = fit.doFit()
        fit.visualize(result, timeout=0.0)
        self.acceptParam(result, targets)
        return result, err

    def step120(self):
        targets = ['UC', 'DWB', 'PRWB']
        fit = MOS_IV_Fit(self.param, targets, 'Step 12')
        fit.setDataSource(     self.IdVg_LW_lin_ba + 
                          0.1 *self.IdVg_MW_lin_ba +
                          0.05*self.IdVg_SW_lin_ba +
                          10. *self.IdVg_LN_lin_ba +
                          2.0*self.IdVg_SN_lin_ba)

        result,err = fit.doFit()
        fit.visualize(result, timeout=0.0)
        self.acceptParam(result, targets)
        return result, err

    def step130(self):
        targets = ['VSAT', 'A0', 'AGS']
        fit = MOS_IV_Fit(self.param, targets, 'Step 13')
        fit.setDataSource(     self.IdVd_LW_b0+ 
                          0.1 *self.IdVd_MW_b0 +
                          0.05*self.IdVd_SW_b0)

        result,err = fit.doFit()
        fit.visualize(result, timeout=0.0)
        self.acceptParam(result, targets)
        return result, err

    def step132(self):
        targets = ['ETA0', 'ETAB']
        fit = MOS_IV_Fit(self.param, targets, 'Step 13-2')
        fit.setDataSource(     self.IdVg_LW_sat_ba+ 
                          0.1 *self.IdVg_MW_sat_ba+
                          0.05*self.IdVg_SW_sat_ba+
                               self.IdVg_LW_d1_ba+ 
                          0.1 *self.IdVg_MW_d1_ba+
                          0.05*self.IdVg_SW_d1_ba)


        result,err = fit.doFit()
        fit.visualize(result, timeout=0.0)
        self.acceptParam(result, targets)
        return result, err

    def step140(self):
        targets = ['B0', 'B1']
        fit = MOS_IV_Fit(self.param, targets, 'Step 14')
        fit.setDataSource(     self.IdVd_LW_b0+ 
                          10. *self.IdVd_LN_b0)

        result,err = fit.doFit()
        fit.visualize(result, timeout=0.0)
        self.acceptParam(result, targets)
        return result, err

    def step150(self):
        targets = ['ETA0', 'DSUB']
        fit = MOS_IV_Fit(self.param, targets, 'Step 15')
        fit.setDataSource(     self.IdVd_LW_b0+ 
                          0.1 *self.IdVd_MW_b0 +
                          0.05*self.IdVd_SW_b0)

        result,err = fit.doFit()
        fit.visualize(result, timeout=0.0)
        self.acceptParam(result, targets)
        return result, err

    def step152(self):
        targets = ['ETA0', 'DSUB', 'PCLM', 'PDIBLC1', 'PDIBLC2', 'DROUT']
        fit = MOS_IV_Fit(self.param, targets, 'Step 15-2')
        fit.setDataSource(     self.IdVd_LW_b0+ 
                          0.1 *self.IdVd_MW_b0 +
                          0.05*self.IdVd_SW_b0)

        result,err = fit.doFit()
        fit.visualize(result, timeout=0.0)
        self.acceptParam(result, targets)
        return result, err

    def step160(self):
        targets = ['ETA0', 'DSUB', 'PCLM', 'PDIBLC1', 'PDIBLC2', 'A1', 'A2']
        fit = MOS_IV_Fit(self.param, targets, 'Step 16')
        fit.setDataSource(     self.IdVd_LW_b0+ 
                          0.1 *self.IdVd_MW_b0 +
                          0.05*self.IdVd_SW_b0)

        result,err = fit.doFit()
        fit.visualize(result, timeout=0.0)
        self.acceptParam(result, targets)
        return result, err


    def step170(self):
        targets = ['A1']
        fit = MOS_IV_Fit(self.param, targets, 'Step 17')
        fit.setDataSource(     self.IdVg_LW_sat_b0+ 
                          0.1 *self.IdVg_MW_sat_b0+
                          0.05*self.IdVg_SW_sat_b0+
                               self.IdVg_LW_d1_b0+ 
                          0.1 *self.IdVg_MW_d1_b0+
                          0.05*self.IdVg_SW_d1_b0)

        result,err = fit.doFit()
        fit.visualize(result, timeout=0.0)
        self.acceptParam(result, targets)
        return result, err

    def step180(self):
        targets = ['KETA', 'ETAB']
        fit = MOS_IV_Fit(self.param, targets, 'Step 18')
        fit.setDataSource(     self.IdVg_LW_sat_ba+ 
                          0.1 *self.IdVg_MW_sat_ba+
                          0.05*self.IdVg_SW_sat_ba+
                               self.IdVg_LW_d1_ba+ 
                          0.1 *self.IdVg_MW_d1_ba+
                          0.05*self.IdVg_SW_d1_ba)

        result,err = fit.doFit()
        fit.visualize(result, timeout=0.0)
        self.acceptParam(result, targets)
        return result, err

    def step190(self):
        targets = ['KETA', 'LKETA']
        fit = MOS_IV_Fit(self.param, targets, 'Step 19')
        fit.setDataSource(     self.IdVg_LW_sat_b0+ 
                          0.1 *self.IdVg_MW_sat_b0+
                          0.05*self.IdVg_SW_sat_b0+
                               self.IdVg_LW_d1_b0+ 
                          0.1 *self.IdVg_MW_d1_b0+
                          0.05*self.IdVg_SW_d1_b0)

        result,err = fit.doFit()
        fit.visualize(result, timeout=0.0)
        self.acceptParam(result, targets)
        return result, err

    def step200(self):
        targets = ['CDSCD']
        fit = MOS_IV_Fit(self.param, targets, 'Step 20')
        fit.setDataSource(     self.IdVg_LW_sat_b0+ 
                          0.1 *self.IdVg_MW_sat_b0+
                          0.05*self.IdVg_SW_sat_b0+
                               self.IdVg_LW_d1_b0+ 
                          0.1 *self.IdVg_MW_d1_b0+
                          0.05*self.IdVg_SW_d1_b0)

        fit.dataSrc.setSubVth(True)
        result,err = fit.doFit()
        fit.visualize(result, timeout=0.0)
        self.acceptParam(result, targets)

        return result, err

    def step205(self):
        targets = ['CDSCB']
        fit = MOS_IV_Fit(self.param, targets, 'Step 20-5')
        fit.setDataSource(     self.IdVg_LW_sat_ba+ 
                          0.1 *self.IdVg_MW_sat_ba+
                          0.05*self.IdVg_SW_sat_ba+
                               self.IdVg_LW_d1_ba+ 
                          0.1 *self.IdVg_MW_d1_ba+
                          0.05*self.IdVg_SW_d1_ba)

        fit.dataSrc.setSubVth(True)
        result,err = fit.doFit()
        fit.visualize(result, timeout=0.0)
        self.acceptParam(result, targets)

        return result, err



param0['TEMP']=25.
param0['TNOM']=25.

proj = MOSp35Proj(param0)
proj.run()

