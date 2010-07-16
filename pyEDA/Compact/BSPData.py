__all__=['CurveSeries', 'BSP_MOSFET_Instance', 'makeInstance']

from pyparsing import *
import numpy as np
from DevMeasData import *

class CurveSeries(object):
    def __init__(self, label):
        self.label = label
        self.curves = []
        self.VScan = None
        self.IOut  = None

    def __len__(self):
        return len(self.curves)

    def __getitem__(self, key):
        '''
        get a curve from the series,

        @param key      index of the curve or a filter function.
                        In the second case, each curve c is passed to the filter fn(c),
                        if fn() yields True, c is returned.
        '''
        if isinstance(key, int):
            idx = key
            if idx>=0 and idx<len(self):
                return self.curves[idx]
            else:
                raise IndexError
        elif isinstance(key, dict):
            cond = key
            for c in self.curves:
                flg = True
                for n,v in cond.iteritems():
                    if abs(c.VConsts[n]-v)>1e-8:
                        flg = False
                if flg:
                    return c
        elif callable(key):
            fn = key
            for c in self.curves:
                if fn(c):
                    return c
        return None

    def __setitem__(self, key, curve):
        '''
        '''
        if not isinstance(key, int):
            raise TypeError
        if not isinstance(curve, IV_Curve):
            raise TypeError

        if self.VScan and self.IOut:
            if not (curve.VScan==self.VScan and curve.IOut==self.IOut):
                raise ValueError
        else:
            self.VScan = curve.VScan
            self.IOut  = curve.IOut

        if key==len(self):
            self.curves.append(curve)
        elif key>=0 and key<len(self):
            self.curves[key] = curve
        else:
            raise IndexError

    def append(self, curve):
        self[len(self)] = curve


class BSP_MOSFET_Instance(MOSFET_Instance):
    def __init__(self, W=None, L=None, T=None):
        super(BSP_MOSFET_Instance, self).__init__(W, L, T)
        self.IVData = {}
    
    def addIVSeries(self, label, series):
        self.IVData[label] = series

    def getIVSeries(self, label):
        return self.IVData[label]

    def listIVSeries(self):
        return self.IVData.keys()

    def loadBSimProFile(self, fname):
        makeInstance(self, fname)


real_number = Combine(Optional(oneOf("+ -")) + Word(nums) + 
               Optional("." + Optional(Word(nums))) +
               Optional(oneOf("e E")+Optional(oneOf("+ -")) +Word(nums)))\
    .setName("real_number")\
    .setParseAction( lambda toks: float(toks[0]) )

identifier = Word(alphas, alphanums+'_')
delimitor  = oneOf(',\t ')
str_label  = Optional(Word(alphanums+'!"#$%&()*+,-./:;<=>?@[\\]^_`|~ '))

param_num  = Group ( identifier('name') + '=' + real_number('value') )
var_input  = Group ( identifier('name') + Suppress('=') + real_number('value') )
var_output = Group ( identifier('name') + Suppress('(') + param_num('param') + Suppress(')') )

file_header = Group(
    'ObjInfo'       + '{' + str_label('ObjInfo')       + '}' +
    'Version'       + '{' + str_label('Version')       + '}' +
    'ModelType'     + '{' + str_label('ModelType')     + '}' +
    'DataType'      + '{' + str_label('DataType')      + '}' +
    'Delimitor'     + '{' + delimitor('DataType')      + '}' +
    'Workingmode'   + '{' + str_label('Workingmode')   + '}' +
    'Instance'      + '{' + delimitedList(var_input)('Instance')  + '}' +
    'Input'         + '{' + delimitedList(identifier)('Input')    + '}' +
    'Output'        + '{' + delimitedList(identifier)('Output')   + '}'
    )

data_row = Group( delimitedList(real_number) + Suppress(lineEnd) )
data_series_header = Group( 
                Suppress('[') + 
                delimitedList(identifier^
                              var_input('var_input')^
                              var_output('var_output')) +
                Suppress(']') )
data_series = Group(data_series_header('header') + OneOrMore(data_row)('data'))

data_set_header = Combine(Suppress('{') + identifier + Suppress('}'))
data_set = Group(data_set_header('name') + OneOrMore(data_series)('data_series'))
        
BSimProFile = file_header('file_header') + OneOrMore(data_set)('data_sets')

def makeInstance(mos, fname):
    res = BSimProFile.parseFile(fname)

    res_header = res.file_header
    for var in res_header.Instance:
        if var.name in ['W', 'L', 'T']:
            mos.__setattr__(var.name, var.value)

    inputVars = res_header.Input.asList()
    outputVars = res_header.Output.asList()

    for dset in res.data_sets:
        curveSeries = CurveSeries(dset.name)
        for dseries in dset.data_series:

            VScan = None
            colScan = 0
            IOut  = None
            VConsts = {}

            data = np.array(dseries.data.asList())

            # first loop: scan voltage and consts
            col = 0
            for var in dseries.header:
                if isinstance(var, str):
                    if var in inputVars:
                        VScan = var
                        colScan = col
                    col = col+1
                elif var.getName()=='var_input':
                    VConsts[var.name] = var.value
                else:
                    col = col+1

            # second loop: output with parameter
            col = 0
            for var in dseries.header:
                if isinstance(var, str):
                    col = col+1
                elif var.getName()=='var_output' and var.name in outputVars:
                    IOut = var.name
                    VConsts_t = VConsts.copy()
                    VConsts_t[var.param.name] = var.param.value
                    
                    curve = MOSFET_IV_Curve(VScan, VConsts_t, IOut, data[:,[colScan, col]])
                    curveSeries.append(curve)

                    col = col+1

        mos.addIVSeries(dset.name, curveSeries)

#######
def parseBSimProFile(fname):
    res = BSimProFile.parseFile(fname)

    res_header = res.file_header
    print res_header

    print '----------'
    for param in res_header.Instance:
        print param.name, param.value
    print res_header.Input
    print res_header.Output

    print '----------'
    for i,dset in enumerate(res.data_sets):
        print '********* Dataset #%d : %s' % (i, dset.name)
        for dseries in dset.data_series:
            for var in dseries.header:
                if isinstance(var, str):
                  print 'simple var', var #simple var name
                elif var.getName()=='var_input':
                  print 'input var ', var.name, var.value
                elif var.getName()=='var_output':
                  print 'output var', var.name, var.param.name, var.param.value

            print len(dseries.data)


if __name__=='__main__':
    mos = BSP_MOSFET_Instance()
    mos.loadBSimProFile('data/NEA1X10.DAT')
    
    for sLabel in mos.listIVSeries():
        print 'Series:', sLabel
        series = mos.getIVSeries(sLabel)
        for i in xrange(len(series)):
            curve = series[i]

            print 'Curve --------------',len(curve)
            for vbias, curr in curve.iterData():
                print vbias, curr


