#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import numpy
import getopt
import os
import version

Usage = """
combine_gTab_SHAPE - Combine SHAPE score from two input gTab files
==================================================================
 \x1b[1mUSAGE:\x1b[0m 
  %s input_1.gTab input_2.gTab output.txt

\x1b[1mVERSION:\x1b[0m
    %s

\x1b[1mAUTHOR:\x1b[0m
    Li Pan

""" % (sys.argv[0], version.Version)


def init():
    
    if len(sys.argv) != 4:
        print(Usage)
        sys.exit(-1)
    
    params = { 'input1': os.path.abspath(sys.argv[1]), 'input2': os.path.abspath(sys.argv[2]), 'output': os.path.abspath(sys.argv[3]) }

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


def combine_gTab(gTab1, gTab2, outFile):
    IN1 = open(gTab1)
    IN2 = open(gTab2)
    OUT = open(outFile, 'w')
    
    gTab_head1 = read_gTab_head(IN1)
    gTab_head2 = read_gTab_head(IN2)
    
    chrID1 = gTab_head1['ChrID'] - 1
    chrpos1 = gTab_head1['ChrPos'] - 1
    chrstrand1 = gTab_head1['Strand'] - 1
    shapepos1 = gTab_head1['Shape'] - 1
    
    chrID2 = gTab_head2['ChrID'] - 1
    chrpos2 = gTab_head2['ChrPos'] - 1
    chrstrand2 = gTab_head2['Strand'] - 1
    shapepos2 = gTab_head2['Shape'] - 1
    
    line1 = IN1.readline()
    line2 = IN2.readline()
    while line1 and line2:
        data1 = line1.strip().split()
        data2 = line2.strip().split()
        
        chr1 = data1[chrID1]
        chr2 = data2[chrID2]
        if chr1 < chr2:
            line1 = IN1.readline()
            continue
        elif chr1 > chr2:
            line2 = IN2.readline()
            continue
        else:
            strand1 = data1[chrstrand1]
            strand2 = data2[chrstrand2]
            if strand1 < strand2:
                line1 = IN1.readline()
                continue
            elif strand1 > strand2:
                line2 = IN2.readline()
                continue
            else:
                pos1 = int(data1[chrpos1])
                pos2 = int(data2[chrpos2])
                if pos1 < pos2:
                    line1 = IN1.readline()
                    continue
                elif pos1 > pos2:
                    line2 = IN2.readline()
                    continue
                else:
                    shape1 = data1[shapepos1]
                    shape2 = data2[shapepos2]
                    if shape1 != '-1' and shape2 != '-1':
                        shape1 = "%.3f" % (float(shape1), )
                        shape2 = "%.3f" % (float(shape2), )
                        OUT.writelines( "%s\t%s\t%s\t%s\t%s\n" % (chr1, strand1, pos1, shape1, shape2) )
                    line1 = IN1.readline()
                    line2 = IN2.readline()
    
    IN1.close()
    IN2.close()
    OUT.close()

def main():
    params = init()
    combine_gTab(params['input1'], params['input2'], params['output'])

if __name__ == '__main__':
    main()

