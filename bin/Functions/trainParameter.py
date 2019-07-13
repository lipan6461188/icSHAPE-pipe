#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import numpy
import getopt
import os
import sklearn, sklearn.metrics
import pandas as pd
import version

dirname = os.path.dirname(os.path.abspath(__file__))

Usage = """
trainParameter - Train different substract factor (default 0.25 in calcSHAPE) and window size
============================================================================================
\x1b[1mUSAGE:\x1b[0m 
  %s -d structure_dot -s sizeFile -D DMSO_tab_files -N NAI_tab_files -o report.pdf

\x1b[1mHELP:\x1b[0m
  -d                    <String>
                                Input a dot file record the standard structures (Usually 18S rRNA)
  -s                    <String>
                                Input a chromosome size file
  -D                    <String>
                                Comma-seperated DMSO tab files
  -N                    <String>
                                Comma-seperated NAI tab files
  -o                    <String>
                                Output a report file (default: report.pdf)
  --debug               <None>
                                Print details

  More options:
  --subFac              <Float,Float,Int>
                                Start, Step, Nums (default: 0.0,0.1,11)
  --Window              <Int,Int,Int>
                                Start, Step, Nums (default: 100,100,5)


\x1b[1mVERSION:\x1b[0m
    %s

\x1b[1mAUTHOR:\x1b[0m
    Li Pan

""" % (sys.argv[0], version.Version)


def init():
    params = { 'dotFile': None, 'DMSOFiles': None, 'NAIFiles': None, 'sizeFile': None, 'outReport': 'report.pdf', 'debug': False, 'subFac': (0.0,0.1,11), 'Window': (100,100,5) }
    opts, args = getopt.getopt(sys.argv[1:], 'hd:D:N:s:o:', ['debug', 'subFac=', 'Window='])
    for op, value in opts:
        if op == '-h':
            print(Usage)
            sys.exit(-1)
        # Basic Parameters
        elif op == '-d':
            params['dotFile'] = os.path.abspath(value)
        elif op == '-D':
            params['DMSOFiles'] = value
        elif op == '-N':
            params['NAIFiles'] = value
        elif op == '-s':
            params['sizeFile'] = os.path.abspath(value)
        elif op == '-o':
            params['outReport'] = os.path.abspath(value)
        
        elif op == '--debug':
            params['debug'] = True
        
        elif op == '--subFac':
            values = value.split(',')
            assert len(values) == 3
            params['subFac'] = ( float(values[0]), float(values[1]), int(values[2]) )
        elif op == '--Window':
            values = value.split(',')
            assert len(values) == 3
            params['Window'] = ( int(values[0]), int(values[1]), int(values[2]) )
        
        else:
            print >>sys.stderr, "Error: unrecognized parameter: "+op
            print(Usage)
            sys.exit(-1)
    
    # check
    if (not params['dotFile']):
        sys.stderr.writelines("Error: Please specify -d"+"\n")
        print(Usage)
        sys.exit(-1)
    if (not params['DMSOFiles']) or (not params['NAIFiles']):
        sys.stderr.writelines("Error: Please specify -D and -N"+"\n")
        print(Usage)
        sys.exit(-1)
    if (not params['sizeFile']):
        sys.stderr.writelines("Error: Please specify -s"+"\n")
        print(Usage)
        sys.exit(-1)
    
    return params

def calculateSHAPE(DFiles, NFiles, sizeFile, windowSize, Subfac, outSHAPE, debug):
    CMD = "%s/calc_sliding_shape TrtCont -D %s -N %s -size %s -wsize %s -sf %s -out %s -noparam"
    CMD = CMD % (dirname, DFiles, NFiles, sizeFile, windowSize, Subfac, outSHAPE)
    if not debug:
        CMD += " 1>/dev/null 2>/dev/null"
    else:
        print(CMD)
    os.system(CMD)

def genomeSHAPEToTransSHAPE(genSHAPE, sizeFile, transSHAPE, debug):
    CMD = "python %s/genSHAPEToTransSHAPE.py -i %s -s %s -o %s"
    CMD = CMD % (dirname, genSHAPE, sizeFile, transSHAPE)
    if not debug:
        CMD += " 1>/dev/null 2>/dev/null"
    else:
        print(CMD)
    os.system(CMD)

def loadicSHAPE(file_name):
    SHAPE = {}
    for line in open(file_name):
        arr = line.strip().split()
        trans_id = arr[0]
        shape = arr[3:]
        SHAPE[ trans_id ] = shape
    return SHAPE

def readDot(inFile):
    structure = {}
    IN = open(inFile)
    line = IN.readline()
    cur_trans_id = ""
    while line:
        if line[0] == '>':
            cur_trans_id = line[1:].strip().split()[0]
            structure[cur_trans_id] = [ IN.readline().strip(), IN.readline().strip() ]
        line = IN.readline()
    IN.close()
    return structure

def Shape_positive_rate(ss_code, shape_list, cutoff):
    Pos_Num = 0
    True_Pos = 0
    False_Pos = 0
    Neg_Num = 0
    for idx, code in enumerate(list(ss_code)):
        if shape_list[idx] != 'NULL':
            if code != ".":
                Pos_Num += 1
                if float(shape_list[idx]) <= cutoff:
                    True_Pos += 1
                else:
                    pass
            else:
                Neg_Num += 1
                if float(shape_list[idx]) <= cutoff:
                    False_Pos += 1
                else:
                    pass
    return 1.0*True_Pos/Pos_Num, 1.0*False_Pos/Neg_Num


def calc_shape_ROC(ss_code, shape_list, step=0.01):
    assert(len(ss_code)==len(shape_list))
    ROC = []
    cutoff = -step
    while cutoff < 1.0 + step:
        TPR, FPR = Shape_positive_rate(ss_code, shape_list, cutoff)
        ROC.append( (FPR, TPR) )
        cutoff += step
    return ROC

def calc_AUC(ROC):
    x = [it[0] for it in ROC]
    y = [it[1] for it in ROC]
    return sklearn.metrics.auc(x, y)

def init_rect(rowNum, colNum, rowNames=[], colNames=[], init_value=None):
    import pandas as pd
    import numpy as np
    
    if colNames:
        assert(len(colNames)==colNum)
    else:
        colNames = np.arange(colNum)
    
    if rowNames:
        assert(len(rowNames)==rowNum)
    else:
        rowNames = np.arange(rowNum)
    
    df = pd.DataFrame(np.zeros((rowNum, colNum)), index=rowNames, columns=colNames)
    if init_value == None:
        return df
    else:
        df.iloc[:,:] = init_value
        return df

def construct_array(start, step, nums):
    array = [ start ]
    i = 1
    while i < nums:
        if type(start) == float:
            array.append( round(array[-1]+step,3) )
        else:
            array.append( array[-1]+step )
        i += 1
    return array

def main():
    params = init()
    
    Dot = readDot(params['dotFile'])
    out_tmp_prefix = params['outReport'].rstrip()
    
    wSizes = construct_array(params['Window'][0], params['Window'][1], params['Window'][2])
    subFacs = construct_array(params['subFac'][0], params['subFac'][1], params['subFac'][2])
    assert subFacs[-1] <= 1.0

    record = init_rect(len(wSizes), len(subFacs), rowNames=wSizes, colNames=subFacs, init_value=0.0)
    
    print("Will train with window size: " + str(wSizes))
    print("Will train with sunstract factor: " + str(subFacs))
    
    for wSize in wSizes:
        for subFac in subFacs:
            sys.stdout.writelines( "Start to train (%s, %s)..." % (wSize, subFac) ); sys.stdout.flush()
            
            outGenomeSHAPE = out_tmp_prefix+".tmp_%s_%s.gTab" % (wSize, subFac)
            outTransSHAPE = out_tmp_prefix+".tmp_%s_%s.shape" % (wSize, subFac)
            calculateSHAPE(params['DMSOFiles'], params['NAIFiles'], params['sizeFile'], wSize, subFac, outGenomeSHAPE, params['debug'])
            genomeSHAPEToTransSHAPE(outGenomeSHAPE, params['sizeFile'], outTransSHAPE, params['debug'])
            
            SHAPE = loadicSHAPE(outTransSHAPE)
            commonTrans = list(set(Dot) & set(SHAPE))
            #print "Common transcript: ", commonTrans
            
            if len(commonTrans) == 0:
                sys.stderr.writelines("Error: no common transcript in SHAPE and dot\n")
                exit(-1)
            
            combineSHAPE = []
            combineDot = ""
            for tid in commonTrans: 
                combineSHAPE += SHAPE[tid]
                combineDot += Dot[tid][1]
            
            assert len(combineSHAPE) == len(combineDot)
            
            ROC = calc_shape_ROC(combineDot, combineSHAPE, step=0.01)
            AUC = calc_AUC(ROC)
            
            record.loc[wSize, subFac] = AUC
            sys.stdout.writelines( " AUC=%.4f\n" % (AUC, ) ); sys.stdout.flush()
            
            os.remove(outGenomeSHAPE)
            os.remove(outTransSHAPE)
    
    print(record)
    record.to_csv(params['outReport']+'.csv')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    plt.figure(figsize=(15, 6))
    sns.heatmap(data=record, annot=True, fmt=".4f", cmap=sns.color_palette("coolwarm", 50))
    plt.xlabel("Substract factor")
    plt.ylabel("Window size")
    plt.title("AUC with different parameters")
    plt.tight_layout()
    plt.savefig(params['outReport'])
    plt.close()

if __name__ == '__main__':
    main()



