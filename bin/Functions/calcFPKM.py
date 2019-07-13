#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import os
import version

Usage = """
calcFPKM - Calculate FPKM with cufflinks
========================================================
\x1b[1mUSAGE:\x1b[0m 
  %s [-p 1] -i mapGenome.bam -o outDir -G annotation.gtf
\x1b[1mHELP:\x1b[0m
  -i                    <String>
                            Input a sorted bam file
  -o                    <String>
                            Output directory
  -G                    <String>
                            Input a GTF annotation file

  More options:
  -p                    <Int>
                            How many threads to use (default: 1)

 \x1b[1mWARNING\x1b[0m
    Input bam file should be sorted by coordination;

\x1b[1mVERSION:\x1b[0m
    %s

\x1b[1mAUTHOR:\x1b[0m
    Li Pan

""" % (sys.argv[0], version.Version)

def init():
    import getopt
    
    Params = { 'inBam': None, 'outDir': None, 'gtfFile': None, 'threads': 1 }
    
    opts, args = getopt.getopt(sys.argv[1:], 'hi:o:G:p:')
    
    for op, value in opts:
        if op == '-h':
            print(Usage)
            exit(-1)
        # Basic Parameters
        elif op == '-i':
            Params['inBam'] = os.path.abspath(value)
        elif op == '-o':
            Params['outDir'] = os.path.abspath(value)
        elif op == '-G':
            Params['gtfFile'] = os.path.abspath(value)
        elif op == '-p':
            Params['threads'] = int(value)

        else:
            sys.stderr.writelines("parameter Error: unrecognized parameter: "+op+"\n")
            print(Usage)
            sys.exit(-1)
    
    if not Params['inBam']:
        sys.stderr.writelines("Error: please specify -i"+"\n")
        print(Usage)
        exit(-1)
    if not Params['outDir']:
        sys.stderr.writelines("Error: please specify -o"+"\n")
        print(Usage)
        exit(-1)
    if not Params['gtfFile']:
        sys.stderr.writelines("Error: please specify -G"+"\n")
        print(Usage)
        exit(-1)
    
    return Params




def main():
    params = init()
    CMD = "cufflinks -p %s %s -o %s -G %s"
    CMD = CMD % (params['threads'], params['inBam'], params['outDir'], params['gtfFile'])
    print("Start to calculate FPKM:\n\t%s" % (CMD, ))
    os.system(CMD)

if __name__ == "__main__":
    main()





