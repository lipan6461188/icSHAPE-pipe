#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import numpy
import getopt
import os
import random
import version
import sklearn, sklearn.metrics
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import General, Structure

Usage = """
evaluateSHAPE - Calculate SHAPE AUC and plot ROC corve with known structure
===========================================================================
\x1b[1mUSAGE:\x1b[0m 
  %s -i trans_shape.out -s structure.dot -o report.pdf
\x1b[1mHELP:\x1b[0m
  -i                    <String>
                                Input transcript-based SHAPE
  -s                    <String>
                                Input known structure with dot-bracked format
  -o                    <String>
                                Output a PDF report (default: no output)

 More options:
  --accessiblity        <String>
                                A file provide the acceessibility for each base
  --min_area            <Float>
                                Provide minimun area to consider each base (default: 5.0)
                                Only useful when --accessFn specified
  
  --ignore_double_strand        The double-stranded bases are not filtered with acceessibility values

 Accessibility file format
    18S 1       T       0.0
    18S 2       A       0.0
    18S 3       C       0.0
    18S 4       C       0.0
    18S 5       T       0.0
    ....

\x1b[1mWARNING:\x1b[0m
    typical dot-bracked format:
        >seq1
        ATCGAGTAGCATCGTACGAT
        ....((((....))))....
        >seq2
        GCTGAGTCAGCTAGCTAGCTAAGA
        .((((....).)).((...))...

\x1b[1mVERSION:\x1b[0m
    %s

\x1b[1mAUTHOR:\x1b[0m
    Li Pan

""" % (sys.argv[0], version.Version)


def init():
    params = { 'inSHAPE': None, 'inDotbracket': None, 'outPDF': None, 
    'step': 0.01, 'accessFn': None, 'minArea': 5.0, 'ignore_double_strand': False }
    opts, args = getopt.getopt(sys.argv[1:], 'hi:s:o:', ['step=','accessiblity=','min_area=', 'ignore_double_strand'])
    for op, value in opts:
        if op == '-h':
            sys.stdout.writelines(Usage+"\n");
            sys.exit(-1)
        # Basic Parameters
        elif op == '-i':
            params['inSHAPE'] = os.path.abspath(value)
        elif op == '-s':
            params['inDotbracket'] = os.path.abspath(value)
        elif op == '-o':
            params['outPDF'] = os.path.abspath(value)
        elif op == '--step':
            params['step'] = float(value)
            assert 0 < params['step'] < 1
        elif op == '--accessiblity':
            params['accessFn'] = os.path.abspath(value)
        elif op == '--min_area':
            params['minArea'] = float(value)
        elif op == '--ignore_double_strand':
            params['ignore_double_strand'] = True
        else:
            sys.stderr.writelines("Error: unrecognized parameter: "+op+"\n")
            sys.stdout.writelines(Usage+"\n");
            sys.exit(-1)
    
    # check
    if (not params['inSHAPE']):
          sys.stderr.writelines("Error: Please specify -i\n")
          sys.stdout.writelines(Usage+"\n")
          sys.exit(-1)
    if (not params['inDotbracket']):
          sys.stderr.writelines("Error: Please specify -s\n")
          sys.stdout.writelines(Usage+"\n")
          sys.exit(-1)
    
    return params

def readAccessibility(inFn):
    """
    The accessibility file format:
    18S 1       T       0.0
    18S 2       A       0.0
    18S 3       C       0.0
    18S 4       C       0.0
    18S 5       T       0.0
    """
    Access = {}
    Fasta = {}
    for line in open(inFn):
        data = line.strip().split()
        if len(data)!=4:
            print("Error: accessibility file format error")
            exit(-1)
        tid = data[0]
        pos = int(data[1])
        base = data[2]
        value = float(data[3]) if data[3]!='NULL' else 'NULL'
        if tid not in Access:
            Access[tid] = []
            Fasta[tid] = []
        while pos > len(Access[tid]):
            Access[tid].append('NULL')
            Fasta[tid].append('N')
        Access[tid][pos-1] = value
        Fasta[tid][pos-1] = base
    for tid in Fasta:
        Fasta[tid] = "".join(Fasta[tid]).upper().replace('U', 'T')
    return Access, Fasta

def filter_seq_shape(shape, dot, params, accessibility=None, aligned_shape_seq=None, aligned_access_seq=None):
    assert len(shape)==len(dot)
    if accessibility:
        assert len(shape)==len(accessibility)
    if aligned_shape_seq:
        assert len(shape)==len(aligned_shape_seq)
    if aligned_access_seq:
        assert len(shape)==len(aligned_access_seq)
    new_shape = []
    new_dot = ""
    for i in range(len(shape)):
        if aligned_shape_seq and aligned_shape_seq[i]=='-':
            continue
        if aligned_access_seq and aligned_access_seq[i]=='-':
            continue
        if accessibility:
            if accessibility[i] == 'NULL':
                new_shape.append(shape[i])
                new_dot += dot[i]
            elif dot[i]!='.' and params['ignore_double_strand']:
                new_shape.append(shape[i])
                new_dot += dot[i]
            elif accessibility[i]>=params['minArea']:
                new_shape.append(shape[i])
                new_dot += dot[i]
            else:
                pass
        else:
            new_shape.append(shape[i])
            new_dot += dot[i]
    return new_shape, new_dot


def main():
    params = init()
    
    dotbracket = General.load_dot(params['inDotbracket'])
    transSHAPE = General.load_shape(params['inSHAPE'])
    
    common_tid = list(set(dotbracket) & set(transSHAPE))
    print("Common transcript in structure file and SHAPE file: "+str(common_tid))
    
    Access = {}
    Access_Fasta = {}
    if params['accessFn']:
        Access, Access_Fasta = readAccessibility(params['accessFn'])
    
    total_shape = []
    total_dot = ""
    for tid in common_tid:
        dot_seq, dot = dotbracket[tid]
        shape = transSHAPE[tid]
        assert len(shape)==len(dot), f"structure length {len(dot)} != shape length {len(shape)} in {tid}"
        
        if tid in Access:
            access = Access[tid]
            #assert len(access)==len(dot), f"structure length {len(dot)} != accessiblity length {len(access)} in {tid}"
            shape_seq = dot_seq.upper().replace('U', 'T')
            aligned_shape_seq, aligned_access_seq = Structure.multi_alignment([shape_seq, Access_Fasta[tid]])
            aligned_shape = Structure.shape_to_alignSHAPE(shape, aligned_shape_seq)
            aligned_access = Structure.shape_to_alignSHAPE(access, aligned_access_seq)
            aligned_dot = Structure.dot_to_alignDot(dot, aligned_shape_seq)
            new_shape, new_dot = filter_seq_shape(aligned_shape, aligned_dot, params, accessibility=aligned_access, aligned_shape_seq=aligned_shape_seq, aligned_access_seq=aligned_access_seq)
        else:
            new_shape, new_dot = filter_seq_shape(shape, dot, params)
        
        unpair_num = new_dot.count('.')
        pair_num = len(new_dot)-unpair_num
        AUC = General.calc_AUC_v2(new_dot, new_shape)
        print(f"{tid} AUC: {AUC:.3}; {pair_num} paired bases {unpair_num} unpaired bases")
        total_shape += new_shape
        total_dot += new_dot
    
    if params['outPDF']:
        ROC = General.calc_shape_structure_ROC(total_dot, total_shape)
        AUC = General.calc_AUC_v2(total_dot, total_shape)
        
        fig = plt.figure(1, figsize=(7, 6))
        ax = fig.add_subplot(111)
        
        x = [ i[0] for i in ROC ]
        y = [ i[1] for i in ROC ]
        ax.plot(x, y, '-', color='black')
        ax.plot([0,1], [0,1], '-', color='gray')
        ax.set_xlim(0,1)
        ax.set_ylim(0,1)
        ax.set_xlabel("False positive rate")
        ax.set_ylabel("True positive rate")
        ax.set_title("ROC of " + ",".join(common_tid)+" AUC="+str(round(AUC, 2)))
        fig.tight_layout()
        fig.savefig(params['outPDF'])

if __name__ == '__main__':
    main()

