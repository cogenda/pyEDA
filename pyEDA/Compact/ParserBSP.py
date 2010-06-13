__all__=[]

from pyparsing import *
import numpy as np
from DevMeasData import MOSFET_Instance
from DevMeasData import MOSFET_IV_Curve
from DevMeasData import CurveSeries


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
  parseBSimProFile('data/NEA1X10.DAT')

