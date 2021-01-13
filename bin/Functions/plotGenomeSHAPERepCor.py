#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import numpy
import getopt
import os
import random
import version
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

Usage = """
plotGenomeSHAPERepCor - Calculate replicate correlation for genome RT
=====================================================================
\x1b[1mUSAGE:\x1b[0m 
  %s -i inputFile -o report.pdf
\x1b[1mHELP:\x1b[0m
  -i                    <String>
                                Input countRT file record the RT (prduced by combine_gTab_SHAPE.py)
  -o                    <String>
                                Output a PDF report (default: report.pdf)

 More options:
  --winSize             <Int>
                                Window size for each replicate calculation (default: 100)

\x1b[1mVERSION:\x1b[0m
    %s

\x1b[1mAUTHOR:\x1b[0m
    Li Pan

""" % (sys.argv[0], version.Version)


def init():
    params = { 'inFile': None, 'outPDF': 'report.pdf', 'winSize': 100 }
    opts, args = getopt.getopt(sys.argv[1:], 'hi:o:', ['winSize='])
    for op, value in opts:
        if op == '-h':
            print(Usage)
            sys.exit(-1)
        # Basic Parameters
        elif op == '-i':
            params['inFile'] = os.path.abspath(value)
        elif op == '-o':
            params['outPDF'] = os.path.abspath(value)
        elif op == '--winSize':
            params['winSize'] = int(value)
        else:
            sys.stderr.writelines("Error: unrecognized parameter: "+op+"\n")
            print(Usage)
            sys.exit(-1)
    
    # check
    if (not params['inFile']):
          sys.stderr.writelines("Error: Please specify -i"+"\n")
          print(Usage)
          sys.exit(-1)
    
    return params


def calcSHAPEReplicateCorrelation(inFile, windowsize=100):
    import scipy
    import scipy.stats
    
    shape1 = []
    shape2 = []
    true_cor = []
    flip1_cor = []
    rand_cor = []
    
    for line in open(inFile):
        data = line.strip().split()
        s1, s2 = float(data[3]), float(data[4])
        shape1.append(s1)
        shape2.append(s2)
        if len(shape1) == windowsize:
            p, v = scipy.stats.pearsonr(shape1, shape2)
            if v < 0.05:
                true_cor.append( p )
            
            p, v = scipy.stats.pearsonr(shape1, [shape2[-1]]+shape2[:-1])
            if v < 0.05:
                flip1_cor.append( p )
            
            shuffle_shape = shape2[:]
            random.shuffle(shuffle_shape)
            p, v = scipy.stats.pearsonr(shape1, shuffle_shape)
            if v < 0.05:
                rand_cor.append( p )
            
            shape1 = []
            shape2 = []
    
    return true_cor, flip1_cor, rand_cor


def boxplot(data_list, ax, width=0.4, labels=None, title=None):
    """
    Example:
        fig = plt.figure(1, figsize=(5, 6))
        ax = fig.add_subplot(111)
        boxplot(data_list, ax)
        fig.show()
    """
    obj = ax.boxplot(x=data_list, showfliers=False, patch_artist=True, widths = width)
    for box in obj['boxes']:
        box.set( color='#7570b3', linewidth=0.5)
        box.set( facecolor = '#1b9e77' )
    
    for whisker in obj['whiskers']:
        whisker.set(color='#7570b3', linewidth=2)
    
    for cap in obj['caps']:
        cap.set(color='#7570b3', linewidth=2)
    
    for median in obj['medians']:
        median.set(color='#b2df8a', linewidth=2)
    
    for flier in obj['fliers']:
        flier.set(marker='o', color='#e7298a', alpha=0.5)
    
    if labels:
        ax.set_xticklabels(labels)
    
    if title:
        ax.set_title(title).set_weight("bold")
    
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    
    return obj


def main():
    params = init()
    true_cor, flip1_cor, rand_cor = calcSHAPEReplicateCorrelation(params['inFile'], windowsize=100)
    
    fig = plt.figure(1, figsize=(5, 6))
    ax = fig.add_subplot(111)
    boxplot([true_cor, flip1_cor, rand_cor], ax, width=0.4, labels=['Pearson\ncorrelation', 'Flip-1bp\ncorrelation', 'Shuffled\ncorrelation'])
    ax.set_ylim(-0.05, 1.05)
    ax.set_title("SHAPE pearson correlation for replicates").set_weight("bold")
    ax.set_ylabel("Pearson correlation for each window")
    fig.savefig(params['outPDF'])

if __name__ == '__main__':
    main()

