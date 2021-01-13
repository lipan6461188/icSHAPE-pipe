#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import os
import version

Usage = """
cleanFq - Removed reads mapped to a given genome with bowtie2
=============================================================
\x1b[1mUSAGE:\x1b[0m 
  %s [-p 1 --mode Local --sam map.sam] -i inFastq -o outFastq -x index
\x1b[1mHELP:\x1b[0m
  -i                    <String>
                            Input a fastq file
  -o                    <String>
                            Output a fastq file can not map to given genome
  -x                    <String>
                            Input a index to exclude

  More options:
  -p                    <Int>
                            How many threads to use (default: 1)
  --mode                <Local/EndToEnd>
                            Mapping the reads to reference with end-to-end mode or local mode (default: Local)
  --sam                 <String>
                            Input a file to save the mapped reads
  --bowparam            <String>
                            More parameters for bowtie2


\x1b[1mVERSION:\x1b[0m
    %s

\x1b[1mAUTHOR:\x1b[0m
    Li Pan

""" % (sys.argv[0], version.Version)

def init():
    import getopt
    
    Params = { 'inFastq': None, 'outFastq': None, 'mode': "Local", 'index': None, 'threads': 1, 'samFile': '/dev/null', 'bow_params':"" }
    
    opts, args = getopt.getopt(sys.argv[1:], 'hi:o:x:p:', ['mode=', 'sam=', 'bowparam='])
    
    for op, value in opts:
        if op == '-h':
            print(Usage)
            exit(-1)
        # Basic Parameters
        elif op == '-i':
            Params['inFastq'] = os.path.abspath(value)
        elif op == '-o':
            Params['outFastq'] = os.path.abspath(value)
        elif op == '-x':
            Params['index'] = os.path.abspath(value)
        elif op == '-p':
            Params['threads'] = int(value)
        elif op == '--mode':
            assert value in ('Local', 'EndToEnd')
            Params['mode'] = value
        elif op == '--sam':
            Params['samFile'] = os.path.abspath(value)
        elif op == '--bowparam':
            Params['bow_params'] = value
        
        else:
            sys.stderr.writelines("parameter Error: unrecognized parameter: "+op+"\n")
            print(Usage)
            sys.exit(-1)
    
    if not Params['inFastq']:
        sys.stderr.writelines("Error: please specify -i"+"\n")
        print(Usage)
        exit(-1)
    if not Params['outFastq']:
        sys.stderr.writelines("Error: please specify -o"+"\n")
        print(Usage)
        exit(-1)
    if not Params['index']:
        sys.stderr.writelines("Error: please specify -x"+"\n")
        print(Usage)
        exit(-1)
    
    return Params

def main():
    params = init()
    if params['mode'] == 'Local':
        Bowtie_More = "--local"
    else:
        Bowtie_More = "--end-to-end"
    
    CMD = """bowtie2 %s -U %s -x %s -p %s --reorder %s | awk '{ if(substr($0,1,1)=="@"||$2==0){print $0 > "%s"}else if($2==4){print "@"$1; print $10; print "+"; print $11;}  }' > %s"""
    CMD = CMD % (Bowtie_More, params['inFastq'], params['index'], params['threads'], params['bow_params'], params['samFile'], params['outFastq'])
    
    print("Start to clean fastq:\n\t%s" % (CMD, ))
    os.system(CMD)

if __name__ == "__main__":
    main()

