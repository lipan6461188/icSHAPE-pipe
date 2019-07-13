#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import numpy
import getopt
import os, sys
import version

Usage = """
genSHAPEToBedGraph - Convert tab-separated genome SHAPE (.gTab) to bedGraph files
================================================================================
\x1b[1mUSAGE:\x1b[0m 
  %s -i input.gTab -t TrtCont -o out_file

\x1b[1mHELP:\x1b[0m
  -i                    <String>
                                Input a gTab file (produced by calc_sliding_shape)
  -o                    <String>
                                Specify a path of directory to save files
  -c                    <Int>
                                Minimun coverage for SHAPE (default: 200 for TrtCont and 100 for Trt)

\x1b[1mVERSION:\x1b[0m
    %s

\x1b[1mAUTHOR:\x1b[0m
    Li Pan

""" % (sys.argv[0], version.Version)



def init():
    params = { }
    opts, args = getopt.getopt(sys.argv[1:], 'hi:o:c:')
    for op, value in opts:
        if op == '-h':
            print(Usage)
            sys.exit(-1)
        # Basic Parameters
        elif op == '-i':
            params['input'] = os.path.abspath(value)
        elif op == '-o':
            params['outdir'] = value
        elif op == '-c':
            params['min_cov'] = int(value)
        else:
            sys.stderr.writelines("Error: unrecognized parameter: "+op+"\n")
            print(Usage)
            sys.exit(-1)
    
    # check
    if ('input' not in params) or ('outdir' not in params):
          sys.stderr.writelines("Error: Please specify -i -o"+"\n")
          print(Usage)
          sys.exit(-1)
    
    return params


def read_gTab_head(IN):
    last_pos = IN.tell()
    
    gTab_head = {}
    
    line = IN.readline()
    while line:
        if not line.startswith('@'):
            ### Check column nums
            data = line.strip().split()
            if len(data) != gTab_head['ColNum']:
                sys.stderr.writelines("Error: actual column number (%s) != labeled number (%s)" % (len(data), gTab_head['ColNum'])+"\n")
                exit(-1)
            IN.seek(last_pos)
            break
        tag, num = line.strip()[1:].split()
        if tag == 'ColNum':
            gTab_head['ColNum'] = int(num)
        elif tag == 'ChrID':
            gTab_head['ChrID'] = int(num)
        elif tag == 'Strand':
            gTab_head['Strand'] = int(num)
        elif tag == 'ChrPos':
            gTab_head['ChrPos'] = int(num)
        elif tag == 'N_RT':
            gTab_head['N_RT'] = int(num)
        elif tag == 'N_BD':
            gTab_head['N_BD'] = int(num)
        elif tag == 'D_RT':
            gTab_head['D_RT'] = int(num)
        elif tag == 'D_BD':
            gTab_head['D_BD'] = int(num)
        elif tag == 'Shape':
            gTab_head['Shape'] = int(num)
        elif tag == 'ShapeNum':
            gTab_head['ShapeNum'] = int(num)
        elif tag == 'WindowShape':
            gTab_head['WindowShape'] = int(num)
        else:
            sys.stderr.writelines("Warning: Unknown head tag: "+line.strip()+"\n")
        
        last_pos = IN.tell()
        line = IN.readline()
    
    return gTab_head

def variance(inputlist):
    import numpy
    if len(inputlist) < 5:
        return max(inputlist) - min(inputlist)
    inputlist = sorted(inputlist)
    return numpy.mean(inputlist[-2:]) - numpy.mean(inputlist[:2])

def prepare_bedGraph(OUT, param, strand):
    name = param['name']
    
    description = param['desc']
    if description == False:
        description = name
    
    if strand == '+':
        name += '_plus'
        description += ' plus strand'
    else:
        name += '_minus'
        description += ' minus strand'
    
    color = param['color']
    
    head = "track type=bedGraph name=\"%s\" description=\"%s\" color=\"%s\" autoScale=off smoothingWindow=off graphType=bar" % (name, description, color)
    OUT.writelines(head+"\n")

def sortBedGraph(inBedGraph):
    bedGraph = []
    headline = ""
    for line in open(inBedGraph):
        if line.startswith('track'):
            headline = line.strip()
            continue
        data = line.strip().split()
        bedGraph.append( (data[0], int(data[1]), data[2], data[3]) )
    bedGraph.sort(key=lambda x: (x[0], x[1]) )
    
    OUT = open(inBedGraph, 'w')
    if headline:
        OUT.writelines(headline+"\n")
    for data in bedGraph:
        OUT.writelines("%s\t%s\t%s\t%s\n" % tuple(data))
    
    OUT.close()

def TrtCont_mode_To_bedGraph(IN, gTab_head, outdir, min_cov=0):
    import numpy
    
    Handlers = {}
    outdir = outdir.rstrip('/') + '/'
    
    SHAPE_PLUS = open(outdir+'shape.plus.bedGraph', 'w')
    SHAPE_MINUS = open(outdir+'shape.minus.bedGraph', 'w')
    SHAPE_VAR_PLUS = open(outdir+'shape_var.plus.bedGraph', 'w')
    SHAPE_VAR_MINUS = open(outdir+'shape_var.minus.bedGraph', 'w')
    N_RT_PLUS = open(outdir+'n_rt.plus.bedGraph', 'w')
    N_RT_MINUS = open(outdir+'n_rt.minus.bedGraph', 'w')
    D_BD_PLUS = open(outdir+'d_BD.plus.bedGraph', 'w')
    D_BD_MINUS = open(outdir+'d_BD.minus.bedGraph', 'w')
    D_RT_PLUS = open(outdir+'d_rt.plus.bedGraph', 'w')
    D_RT_MINUS = open(outdir+'d_rt.minus.bedGraph', 'w')
    
    prepare_bedGraph(SHAPE_PLUS, {'name': 'TrtCont', 'desc': 'TrtCont', 'color': '83,169,102'}, '+')
    prepare_bedGraph(SHAPE_MINUS, {'name': 'TrtCont', 'desc': 'TrtCont', 'color': '83,169,102'}, '-')
    prepare_bedGraph(SHAPE_VAR_PLUS, {'name': 'TrtCont_var', 'desc': 'TrtCont_var', 'color': '83,169,102'}, '+')
    prepare_bedGraph(SHAPE_VAR_MINUS, {'name': 'TrtCont_var', 'desc': 'TrtCont_var', 'color': '83,169,102'}, '-')
    
    chr_i = gTab_head['ChrID']-1
    strand_i = gTab_head['Strand']-1
    pos_i = gTab_head['ChrPos']-1
    N_RT_i = gTab_head['N_RT']-1
    N_BD_i = gTab_head['N_BD']-1
    D_RT_i = gTab_head['D_RT']-1
    D_BD_i = gTab_head['D_BD']-1
    shape_i = gTab_head['Shape']-1
    window_i = gTab_head['WindowShape']-1
    
    for line in IN:
        #if not line.startswith('chr'):
        #    continue
        if line[0] == '@': continue
        data = line.strip().split()
        pos = int(data[pos_i])
        N_RT, N_BD = int(data[N_RT_i]), int(data[N_BD_i])
        D_RT, D_BD = int(data[D_RT_i]), int(data[D_BD_i])
        shape_Score = data[shape_i]
         
        locus = (data[chr_i], data[strand_i])
        
        if D_BD >= min_cov and shape_Score != '-1':
            window_shape = [ float(it) for it in data[window_i].strip(',').split(',') if it != '-1' ]
            if len(window_shape) == 0:
                var = 0
            else:
                var = variance(window_shape)
            
            if locus[1] == '+':
                SHAPE_PLUS.writelines("%s\t%s\t%s\t%.3f\n" % (locus[0], pos-1, pos, float(shape_Score)))
                SHAPE_VAR_PLUS.writelines("%s\t%s\t%s\t%.3f\n" % (locus[0], pos-1, pos, var))
            elif locus[1] == '-':
                SHAPE_MINUS.writelines("%s\t%s\t%s\t%.3f\n" % (locus[0], pos-1, pos, float(shape_Score)))
                SHAPE_VAR_MINUS.writelines("%s\t%s\t%s\t%.3f\n" % (locus[0], pos-1, pos, var))
        
        if D_BD > 20:
            if locus[1] == '+':
                D_BD_PLUS.writelines("%s\t%s\t%s\t%s\n" % (locus[0], pos-1, pos, D_BD))
                N_RT_PLUS.writelines("%s\t%s\t%s\t%s\n" % (locus[0], pos-1, pos, N_RT))
                D_RT_PLUS.writelines("%s\t%s\t%s\t%s\n" % (locus[0], pos-1, pos, D_RT))
            elif locus[1] == '-':
                D_BD_MINUS.writelines("%s\t%s\t%s\t%s\n" % (locus[0], pos-1, pos, D_BD))
                N_RT_MINUS.writelines("%s\t%s\t%s\t%s\n" % (locus[0], pos-1, pos, N_RT))
                D_RT_MINUS.writelines("%s\t%s\t%s\t%s\n" % (locus[0], pos-1, pos, D_RT))
    
    SHAPE_PLUS.close(); SHAPE_VAR_PLUS.close();
    SHAPE_MINUS.close(); SHAPE_VAR_MINUS.close();
    D_BD_PLUS.close(); N_RT_PLUS.close(); D_RT_PLUS.close();
    D_BD_MINUS.close(); N_RT_MINUS.close(); D_RT_MINUS.close();
    
    sortBedGraph(outdir+'shape.plus.bedGraph')
    sortBedGraph(outdir+'shape.minus.bedGraph')
    sortBedGraph(outdir+'shape_var.plus.bedGraph')
    sortBedGraph(outdir+'shape_var.minus.bedGraph')

def Trt_mode_To_bedGraph(IN, gTab_head, outdir, min_cov=0):
    import numpy
    
    Handlers = {}
    outdir = outdir.rstrip('/') + '/'
    
    SHAPE_PLUS = open(outdir+'shape.plus.bedGraph', 'w')
    SHAPE_MINUS = open(outdir+'shape.minus.bedGraph', 'w')
    SHAPE_VAR_PLUS = open(outdir+'shape_var.plus.bedGraph', 'w')
    SHAPE_VAR_MINUS = open(outdir+'shape_var.minus.bedGraph', 'w')
    N_RT_PLUS = open(outdir+'n_rt.plus.bedGraph', 'w')
    N_RT_MINUS = open(outdir+'n_rt.minus.bedGraph', 'w')
    N_BD_PLUS = open(outdir+'n_BD.plus.bedGraph', 'w')
    N_BD_MINUS = open(outdir+'n_BD.minus.bedGraph', 'w')
    
    prepare_bedGraph(SHAPE_PLUS, {'name': 'Trt', 'desc': 'Trt', 'color': '202,75,78'}, '+')
    prepare_bedGraph(SHAPE_MINUS, {'name': 'Trt', 'desc': 'Trt', 'color': '202,75,78'}, '-')
    prepare_bedGraph(SHAPE_VAR_PLUS, {'name': 'Trt_var', 'desc': 'Trt_var', 'color': '202,75,78'}, '+')
    prepare_bedGraph(SHAPE_VAR_MINUS, {'name': 'Trt_var', 'desc': 'Trt_var', 'color': '202,75,78'}, '-')
    
    chr_i = gTab_head['ChrID']-1
    strand_i = gTab_head['Strand']-1
    pos_i = gTab_head['ChrPos']-1
    N_RT_i = gTab_head['N_RT']-1
    N_BD_i = gTab_head['N_BD']-1
    shape_i = gTab_head['Shape']-1
    window_i = gTab_head['WindowShape']-1

    for line in IN:
        #if not line.startswith('chr'):
        #    continue
        if line[0] == '@': continue
        data = line.strip().split()
        pos = int(data[pos_i])
        N_RT, N_BD = int(data[N_RT_i]), int(data[N_BD_i])
        shape_Score = data[shape_i]
         
        locus = (data[chr_i], data[strand_i])
        
        if N_BD >= min_cov and shape_Score != '-1':
            window_shape = [ float(it) for it in data[window_i].strip(',').split(',') if it != '-1' ]
            if len(window_shape) == 0:
                var = 0
            else:
                var = variance(window_shape)
            
            if locus[1] == '+':
                SHAPE_PLUS.writelines("%s\t%s\t%s\t%.3f\n" % (locus[0], pos-1, pos, float(shape_Score)))
                SHAPE_VAR_PLUS.writelines("%s\t%s\t%s\t%.3f\n" % (locus[0], pos-1, pos, var))
            elif locus[1] == '-':
                SHAPE_MINUS.writelines("%s\t%s\t%s\t%.3f\n" % (locus[0], pos-1, pos, float(shape_Score)))
                SHAPE_VAR_MINUS.writelines("%s\t%s\t%s\t%.3f\n" % (locus[0], pos-1, pos, var))
        
        if N_BD > 20:
            if locus[1] == '+':
                N_BD_PLUS.writelines("%s\t%s\t%s\t%s\n" % (locus[0], pos-1, pos, N_BD))
                N_RT_PLUS.writelines("%s\t%s\t%s\t%s\n" % (locus[0], pos-1, pos, N_RT))
            elif locus[1] == '-':
                N_BD_MINUS.writelines("%s\t%s\t%s\t%s\n" % (locus[0], pos-1, pos, N_BD))
                N_RT_MINUS.writelines("%s\t%s\t%s\t%s\n" % (locus[0], pos-1, pos, N_RT))
    
    SHAPE_PLUS.close(); SHAPE_VAR_PLUS.close();
    SHAPE_MINUS.close(); SHAPE_VAR_MINUS.close();
    N_RT_PLUS.close(); N_BD_PLUS.close();
    N_RT_MINUS.close(); N_BD_MINUS.close();
    
    sortBedGraph(outdir+'shape.plus.bedGraph')
    sortBedGraph(outdir+'shape.minus.bedGraph')
    sortBedGraph(outdir+'shape_var.plus.bedGraph')
    sortBedGraph(outdir+'shape_var.minus.bedGraph')


def main():
    parameters = init()
    
    IN = open(parameters['input'], 'r')
    gTab_head = read_gTab_head(IN)
    mode = 'TrtCont' if 'D_RT' in gTab_head else 'Trt'
    
    print("Start to transform genome SHAPE to bedGraph files...")
    if mode == 'Trt':
        min_cov = parameters['min_cov'] if 'min_cov' in parameters else 100
        Trt_mode_To_bedGraph(IN, gTab_head, parameters['outdir'], min_cov)
    
    elif mode == 'TrtCont':
        min_cov = parameters['min_cov'] if 'min_cov' in parameters else 200
        TrtCont_mode_To_bedGraph(IN, gTab_head, parameters['outdir'], min_cov)
    
    print("Success!")

main();



