#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import numpy
import getopt
import os
import re
from matplotlib.gridspec import GridSpec
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import version

def gene_type(raw_type):
    valid_gene_type = ('pseudogene', 'snoRNA', 'snRNA', 'miRNA', 'misc_RNA', 'rRNA', 'mRNA')
    lncRNA_class = ('3prime_overlapping_ncrna','antisense','lincRNA','non_coding','sense_intronic','sense_overlapping','processed_transcript')
    if raw_type in valid_gene_type: return raw_type;
    if re.match('.*pseudogene',raw_type): return 'pseudogene';
    if raw_type == 'protein_coding': return 'mRNA';
    if raw_type in lncRNA_class: return 'lncRNA';
    return 'other'

def getGeneTypeColor(gType):
    #colors = [(0.12156862745098039, 0.4666666666666667, 0.7058823529411765), (1.0, 0.4980392156862745, 0.054901960784313725), (0.17254901960784313, 0.6274509803921569, 0.17254901960784313), (0.8392156862745098, 0.15294117647058825, 0.1568627450980392), (0.5803921568627451, 0.403921568627451, 0.7411764705882353), (0.5490196078431373, 0.33725490196078434, 0.29411764705882354), (0.8901960784313725, 0.4666666666666667, 0.7607843137254902), (0.4980392156862745, 0.4980392156862745, 0.4980392156862745), (0.7372549019607844, 0.7411764705882353, 0.13333333333333333), (0.09019607843137255, 0.7450980392156863, 0.8117647058823529)]
    colors = ['#f44336', '#2196f3', '#3f51b5', '#9c27b0', '#009688', '#673ab7', '#ff5722', '#e91e63', '#ff9800'] * 2
    gTypes = ['mRNA', 'lncRNA', 'pseudogene', 'snoRNA', 'snRNA', 'miRNA', 'misc_RNA', 'rRNA']
    if gType == 'other':
        return 'gray'
    i = gTypes.index(gType)
    return colors[i]

def count_valid_ratio(SHAPE):
    ratio_list = []
    for tid in SHAPE:
        total_len = len(SHAPE[tid])
        valid_len = total_len - SHAPE[tid].count('NULL')
        ratio_list.append( 1.0*valid_len/total_len )
    ratio_list.sort()
    return ratio_list

def pie(num_list, textColors=None, textColorFollows=False, explodes=None, colors=None, labels=None, format=None, labeldistance=1.05, fontweight=12):
    import matplotlib.pyplot as plt
    obj = plt.pie(num_list, explode=explodes, colors=colors, labels=labels, shadow=False, autopct=format, labeldistance=labeldistance)
    if textColorFollows:
        for b, t in zip(obj[0], obj[1]):
            t.set_color( b.get_facecolor() )
            t.set_fontsize(12)
    if textColors:
        i = 0
        for b, t in zip(obj[0], obj[1]):
            t.set_color( textColors[i] )
            i += 1
    
    return obj

def rename(rawDict):
    newDict = {}
    for k in rawDict:
        GType = gene_type(k)
        newDict[GType] = newDict.get(GType, 0) + rawDict[k]
    return newDict

def prepare_gene_pie_elements(inDict):
    newDict = rename(inDict)
    
    geneCount = []
    for k in newDict:
        geneCount.append( (k, newDict[k]) )
    geneCount.sort(key=lambda x: x[1], reverse=True)
    gTypes = [ it[0] for it in geneCount ]
    counts = [ it[1] for it in geneCount ]
    ratios = [ 1.0*counts[k]/sum(counts) for k in range(len(counts)) ]
    colors = [ getGeneTypeColor(it) for it in gTypes ]
    
    labels = gTypes[:]
    for i in range(len(ratios)):
        if ratios[i] < 0.01:
            labels[i] = ""
    
    df = []
    for gtype,label,count,color in zip(gTypes, labels, counts, colors):
        df.append( (gtype,label,count,color) )
    df = pd.DataFrame(df, columns=['gtype','label','count','color'])
    df['ratio'] = df['count']/df['count'].sum()
    df['ratioText'] = [ "%.2f%%" % (100*it, ) for it in list(df['ratio']) ]
    
    return df

def classify_trans(SHAPE, Parser):
    TransDict = {}
    BaseDict = {}
    for tid in SHAPE:
        try:
            ft = Parser.getTransFeature(tid)
            gt = ft['gene_type']
            TransDict[ gt ] = TransDict.get(gt, 0) + 1
            valid_num = len(SHAPE[tid]) - SHAPE[tid].count("NULL")
            BaseDict[ gt ] = BaseDict.get(gt, 0) + valid_num
        except KeyError:
            continue
    return TransDict, BaseDict

def cdf(data_list, color='red', topdown=False, label=None, plotMedian=True):
    import re
    import numpy as np
    import matplotlib.pyplot as plt
    
    data_list = np.sort(data_list)
    if topdown:
        p = 1 - 1.0 * np.arange(len(data_list))/(len(data_list) - 1)
    else:
        p = 1.0 * np.arange(len(data_list))/(len(data_list) - 1)
    plt.plot(data_list, p, color=color, label=label)
    if plotMedian:
        median_x = data_list[ int(len(data_list)/2) ]
        median_y = p[ int(len(p)/2) ]
        plt.plot([median_x], [median_y], 'bo')
        plt.axvline(x=median_x, ymin=0, ymax=1, linewidth=2, color='r')

def PlotTransSHAPEStatistics(SHAPE, Parser, outPDF):
    import matplotlib.pyplot as plt
    
    TransDict, BaseDict = classify_trans(SHAPE, Parser)
    df_trans_element = prepare_gene_pie_elements(TransDict)
    df_base_element = prepare_gene_pie_elements(BaseDict)
    
    ratio_list = count_valid_ratio(SHAPE)
    
    print("Plot figures...")
    fig = plt.figure(figsize=(14, 14))
    grids = GridSpec(4, 3)
    
    #### Block 1
    plt.subplot(grids[0, 0], aspect=1)
    objs1 = pie(list(df_trans_element['count']), colors=list(df_trans_element['color']), labels=list(df_trans_element['label']), textColorFollows=True, labeldistance=1.01)
    plt.title("Transcripts statistics", fontdict={'fontweight': 'bold'})
    
    #### Block 2
    plt.subplot(grids[0, 1:], aspect=1)
    plt.axis('off')
    plt.axis('tight')
    cellColors = [ [it, it] for it in list(df_trans_element['color']) ]
    table1 = plt.table(cellText=df_trans_element.loc[:, ('gtype','count')].values, loc='center', cellLoc='left', colLabels=['geneType', 'count'], cellColours=cellColors)
    cell_dict = table1.get_celld()
    for k in cell_dict:
        cell_dict[k].set_width(0.3)
        cell_dict[k].set_height(0.1)
        cell_dict[k].set_linewidth(0.2)
    
    #### Block 3
    plt.subplot(grids[1, 0], aspect=1)
    objs2 = pie(list(df_base_element['count']), colors=list(df_base_element['color']), labels=list(df_base_element['label']), textColorFollows=True, labeldistance=1.01)
    plt.title("Bases statistics", fontdict={'fontweight': 'bold'})
    
    #### Block 4
    plt.subplot(grids[1, 1:], aspect=1)
    plt.axis('off')
    plt.axis('tight')
    cellColors = [ [it, it] for it in list(df_base_element['color']) ]
    table2 = plt.table(cellText=df_base_element.loc[:, ('gtype','count')].values, loc='center', cellLoc='left', colLabels=['geneType', 'count'], cellColours=cellColors)
    cell_dict = table2.get_celld()
    for k in cell_dict:
        cell_dict[k].set_width(0.3)
        cell_dict[k].set_height(0.1)
        cell_dict[k].set_linewidth(0.2)
    
    #### Block 5
    plt.subplot(grids[2, 0], aspect=1)
    cdf(ratio_list, color='black', topdown=True, label=None, plotMedian=True)
    plt.xlabel("Cover ratio")
    plt.ylabel("Sorted transcript")
    plt.title("Covered ratio of sorted transcripts")
    
    #### Block 7
    ax = fig.add_subplot(grids[2, 1:])
    
    GINI_list, GINI_type = genetype_gini(SHAPE, Parser)
    violin(ax, GINI_list, GINI_type, colors=["#C44E52"]*len(GINI_type), rem_ext=0)
    ax.set_ylabel("Gini (Reactivity score)")
    
    #### Block 6
    ax = fig.add_subplot(grids[3, :])
    
    start_codon, stop_codon = calc_period(SHAPE, Parser)
    start = [ int(sum(it)/len(it)) for it in start_codon ]
    stop = [ int(sum(it)/len(it)) for it in stop_codon ]
    
    ax.plot( range(1, 52), start, '-' )
    ax.plot( range(61, 112), stop, '-' )
    
    ax.axvline(x=26, ymin=0, ymax=1, linewidth=1, color='gray')
    ax.axvline(x=86, ymin=0, ymax=1, linewidth=1, color='gray')
    ax.set_xticks([1, 26, 51, 61, 86, 111])
    ax.set_xticklabels(['-25', 'start', '25', '-25', 'stop', '25'])
    ax.set_xlabel("Nucleotide position")
    ax.set_ylabel("Reactivity score")
    
    fig.savefig(outPDF)
    grids.update(wspace=0.05, hspace=0.2)
    plt.close()

#############################
### Code for Gini calculation
#############################

def adjacent_values(vals, q1, q3):
    import numpy
    vals = sorted(vals)
    
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = numpy.clip(upper_adjacent_value, q3, vals[-1])
    
    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = numpy.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


def set_axis_style(ax, labels):
    import numpy
    
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(numpy.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)

def violin(ax, data_list, labels, colors=None, rem_ext=0):
    """
    fig, axs = tools.plt.subplots(nrows=5, ncols=1, figsize=(6, 12), sharey=True)
    
    axs[0].set_title('smartSHAPE 1 ng')
    axs[0].set_ylabel('smart SHAPE score')
    data = [ [],[],[],[] ]
    violin(axs[0], data, ['A', 'T', 'C', 'G'])
    
    fig.tight_layout()
    fig.show()
    """
    import numpy
    
    if colors == None:
        colors = ['#D43F3A'] * len(data_list)
    else:
        assert len(colors) == len(data_list)
        colors = colors[::-1]
    
    if rem_ext:
        assert 0.0 <= rem_ext <= 0.5
        import copy
        data_list = copy.deepcopy(data_list)
        for idx in range(len(data_list)):
            data_list[idx].sort()
            remNum = int(len(data_list[idx]) * rem_ext)
            start = remNum; end = len(data_list[idx]) - remNum
            data_list[idx] = data_list[idx][start:end]
    
    parts = ax.violinplot(data_list, showmeans=False, showmedians=False, showextrema=False)
    
    for pc in parts['bodies']:
        pc.set_facecolor(colors.pop())
        pc.set_edgecolor('black')
        pc.set_alpha(1)
    
    quartile1 = []; medians = []; quartile3 = []
    for data in data_list:
        quartile1.append( numpy.percentile(data, 25) )
        medians.append( numpy.percentile(data, 50) )
        quartile3.append( numpy.percentile(data, 75) )
    
    whiskers = numpy.array([
        adjacent_values(sorted_array, q1, q3)
        for sorted_array, q1, q3 in zip(data_list, quartile1, quartile3)])
    
    whiskersMin, whiskersMax = whiskers[:, 0], whiskers[:, 1]
    
    inds = numpy.arange(1, len(medians) + 1)
    ax.scatter(inds, medians, marker='o', color='white', s=20, zorder=3)
    ax.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=8)
    ax.vlines(inds, whiskersMin, whiskersMax, color='k', linestyle='-', lw=1)
    
    set_axis_style(ax, labels)


def calcGINI(my_list, valid_cutoff=20):
    def GINI(list_of_values):
        length = len(list_of_values)
        total = sum(list_of_values)
        if total == 0: 
            return None
        Sorted_Array = sorted(list_of_values)
        accum, giniB = 0, 0
        for i in Sorted_Array:
            accum += i
            giniB += accum - i / 2.0
        fair_area = accum * length / 2.0
        return (fair_area - giniB) / fair_area
    
    floatArr = []
    for item in my_list:
        if item != 'NULL':
            floatArr.append(float(item))
    
    if len(floatArr) > valid_cutoff:
        return GINI(floatArr)
    return None


def genetype_gini(SHAPE, Parser):
    GINI = { 'pseudogene':[], 'snoRNA':[], 'snRNA':[], 'miRNA':[], 'misc_RNA':[], 'mRNA':[], 'lncRNA': [], 'UTR5':[], 'CDS':[], 'UTR3': [] }
    
    for tid in SHAPE:
        try:
            ft = Parser.getTransFeature(tid)
        except KeyError:
            continue
        gt = gene_type(ft['gene_type'])
        if gt == 'other':
            continue
        
        if gt == 'rRNA':
            continue
        
        shape = SHAPE[tid]
        if gt == 'mRNA':
            cds_s, cds_e = ft['cds_start'], ft['cds_end']
            gini_5 = calcGINI(shape[:cds_s])
            gini_cds = calcGINI(shape[cds_s:cds_e])
            gini_3 = calcGINI(shape[cds_e:])
            
            if gini_5: GINI['UTR5'].append( gini_5 )
            if gini_cds: GINI['CDS'].append( gini_cds )
            if gini_3: GINI['UTR3'].append( gini_3 )
        
        gini = calcGINI(shape)
        if gini: GINI[gt].append( gini )
    
    
    GINI_type = ['mRNA', 'UTR5', 'CDS', 'UTR3', 'pseudogene', 'lncRNA', 'miRNA', 'snoRNA', 'misc_RNA']
    new_GINI_type = ["mRNA", "5'UTR", 'CDS', "3'UTR", 'Pseudogene', 'lncRNA', 'miRNA', 'snoRNA', 'misc_RNA']
    GINI_list = [ GINI[gt] for gt in GINI_type ]
    for i in range(len(GINI_list)):
        GINI_type[i] = new_GINI_type[i] + "\n("+str(len(GINI_list[i]))+")"
        if len(GINI_list[i]) < 20:
            GINI_list[i] = [0.5] * 20
    
    return GINI_list, GINI_type


#############################
### Code CDS period
#############################

def calc_period(SHAPE, Parser):
    start_codon = []
    stop_codon = []
    
    for i in range(51):
        start_codon.append([])
        stop_codon.append([])
    
    for tid in SHAPE:
        try:
            ft = Parser.getTransFeature(tid)
        except KeyError:
            continue
        if gene_type(ft['gene_type']) == 'mRNA':
            shape = SHAPE[tid]
            le = ft['trans_len']
            if le < 100: continue
            cds_s, cds_e = ft['cds_start'], ft['cds_end']
            for i in range( max(cds_s-25,0), min(cds_s+26,le) ):
                if shape[i] != 'NULL':
                    start_codon[ i-max(cds_s-25,0) ].append( float(shape[i]) )
            for i in range( max(cds_e-25,0), min(cds_e+26,le) ):
                if shape[i] != 'NULL':
                    stop_codon[ i-max(cds_e-25,0) ].append( float(shape[i]) )
    
    for i in range(51):
        if len(start_codon[i]) < 20:
            start_codon[i] = [0.2] * 20
        if len(stop_codon[i]) < 20:
            stop_codon[i] = [0.2] * 20
    
    return start_codon, stop_codon


