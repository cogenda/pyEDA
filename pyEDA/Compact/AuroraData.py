__all__=['ARR_MOSFET_Instance', 'AuroraFile']

from pyparsing import *
import numpy as np
from DevMeasData import *
import re
import string

class ARR_MOSFET_Instance(MOSFET_Instance):
    def __init__(self, W=None, L=None, T=None):
        super(ARR_MOSFET_Instance, self).__init__(W, L, T)
        self.columns = []
        self.data = None

    def setData(self, cols, data):
        self.columns = cols
        self.data = data

    def getCurve(self, VScan, VConsts, IOut):
        colScan = self.columns.index(VScan)
        colOut  = self.columns.index(IOut)
        VConstsC = {}
        for k,v in VConsts.iteritems():
            idx = self.columns.index(k)
            VConstsC[idx] = v

        cnt, arr_size=0, 100
        res = np.zeros((arr_size),dtype=[('v', np.float),('i', np.float)])
        rows,cols = np.shape(self.data)
        for i in xrange(rows):
            agree=True
            for j in xrange(cols):
                if VConstsC.has_key(j) and abs(VConstsC[j]-self.data[i,j])>1e-8:
                    agree=False
            if agree:
                if cnt>=arr_size:
                    arr_size = arr_size*2
                    res = np.resize(res, arr_size)
                res[cnt] = tuple(self.data[i, [colScan, colOut]])
                cnt=cnt+1
        res = np.resize(res, cnt)
        res1 = np.zeros((cnt,2))
        for i in xrange(cnt):
            res1[i,:] = list(res[i])

        return MOSFET_IV_Curve(VScan, VConsts, IOut, res1)

class AuroraFile(object):
    INS_HEADER=0
    INS_VAR=1
    DATA_HEADER=2
    DATA_TABLE=3
    re_ins_h1 = re.compile('\$ Aurora File:.*')
    re_ins_h2 = re.compile('\$ Atem File:\s*(?P<name>.*)')
    re_real_num = '[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?'
    re_var = re.compile('VARIABLE\s(?P<name>\w+)\s*=\s*(?P<val>'+ re_real_num +')')

    col_map = {'VGS':'Vgs', 'VBS':'Vbs', 'VDS':'Vds', 'ID':'Id'}

    def __init__(self, fname):
        self.file = open(fname)
        self.instances = {}

    def make_instances(self):
        self.state = self.INS_HEADER
        line=self.file.readline()

        ins_vars={}
        ins_name=None
        cols = []
        data = None
        cnt, arr_size=0, 100

        while(len(line)>0):
            if self.state==self.INS_HEADER:
                m1 = self.re_ins_h1.match(line)
                m2 = self.re_ins_h2.match(line)
                if m1:
                    pass
                elif m2:
                    ins_name = m2.group('name')
                    self.state = self.INS_VAR
                else:
                    raise Exception
            elif self.state==self.INS_VAR:
                m = self.re_var.match(line)
                if m:
                    ins_vars[m.group('name')] = m.group('val')
                elif line=='\n':  # end of header section
                    self.state = self.DATA_HEADER
                else:
                    raise Exception
            elif self.state==self.DATA_HEADER:
                m = line.split(None)
                for col in m[1:]:
                    cols.append(self.col_map[string.upper(col)])
                data = np.zeros((arr_size,len(cols)))
                self.state=self.DATA_TABLE
            elif self.state==self.DATA_TABLE:
                m = line.split(None)
                if len(m)==len(cols):
                    if cnt>=arr_size:
                        arr_size = arr_size*2
                        data = np.resize(data, (arr_size,len(cols)))
                    data[cnt,:] = m
                    cnt+=1
                elif line=='\n': # end of data
                    data = np.resize(data, (cnt,len(cols)))

                    W = ins_vars.get('W', 1.0e-6)
                    L = ins_vars.get('L', 1.0e-6)
                    T = ins_vars.get('T', 25)
                    ins = ARR_MOSFET_Instance(W, L, T)
                    ins.setData(cols, data)
                    self.instances[ins_name]=ins
                    
                    # reset data
                    ins_name=None
                    ins_vars={}
                    cols = []
                    data = None
                    cnt, arr_size=0, 100
                    self.state = self.INS_HEADER
                else:
                    raise Exception

            line=self.file.readline()

        

                    


