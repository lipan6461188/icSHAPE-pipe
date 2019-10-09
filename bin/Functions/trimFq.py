#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import os
import version
import gzip
import random
import IOFile

Usage = """
trimFq - Trim 5' barcode and 3' adaptor
=========================================
\x1b[1mUSAGE:\x1b[0m 
  %s [-p 1 -m 25] -i inFastq -o outFastq -l trim_leading_len -a adaptor
\x1b[1mHELP:\x1b[0m
  -i                    <String>
                            Input a fastq file
  -o                    <String>
                            Output a processed fastq file
  -l                    <Int>
                            How many base to trim 5' of reads (typical icSHAPE: 13)
  -a                    <String>
                            Input a adaptor fasta file

  More options:
  -p                    <Int>
                            How many threads to use (default: 1)
  -m                    <Int>
                            Minimun length to preserve (default: 25)
  --simplify            <None>
                        Simplify fastq information and save space:
                        Example: 
                            line 1n:"@SRR1057939.1 other information" => "@SRR1057939.1"
                            line 3n:"+ other information" => "+"

\x1b[1mVERSION:\x1b[0m
    %s

\x1b[1mAUTHOR:\x1b[0m
    Li Pan

""" % (sys.argv[0], version.Version)

dirname = os.path.dirname(os.path.abspath(__file__))
trimmomatic = os.path.join(dirname, 'trimmomatic-0.38.jar')

def init():
    import getopt
    
    Params = { 'inFastq': None, 'outFastq': None, 'leading': 13, 'adator': None, 'threads': 1, 'minLen': 25, 'simplify':False }
    
    opts, args = getopt.getopt(sys.argv[1:], 'hi:o:l:a:p:m:', ['simplify'])
    
    for op, value in opts:
        if op == '-h':
            print(Usage)
            exit(-1)
        # Basic Parameters
        elif op == '-i':
            Params['inFastq'] = os.path.abspath(value)
        elif op == '-o':
            Params['outFastq'] = os.path.abspath(value)
        elif op == '-l':
            Params['leading'] = int(value)
        elif op == '-a':
            Params['adator'] = os.path.abspath(value)
        elif op == '-p':
            Params['threads'] = int(value)
        elif op == '-m':
            Params['minLen'] = int(value)
        elif op == '--simplify':
            Params['simplify'] = True
        
        else:
            sys.stderr.writelines( "parameter Error: unrecognized parameter: "+op+"\n" )
            print(Usage)
            sys.exit(-1)
    
    if not Params['inFastq']:
        sys.stderr.writelines( "Error: please specify -i"+"\n" )
        print(Usage)
        exit(-1)
    if not Params['outFastq']:
        sys.stderr.writelines( "Error: please specify -o"+"\n" )
        print(Usage)
        exit(-1)
    if not Params['adator']:
        sys.stderr.writelines( "Error: please specify -a"+"\n" )
        print(Usage)
        exit(-1)
    
    return Params


def need_to_simplify(inFQFn):
    """
    Test the head 100 reads to see if it is possible to simplify
    """
    IN = IOFile.IOFile(inFQFn, 'r')
    readCount = 0
    key = IN.readline()
    while key:
        seq = IN.readline()
        tag = IN.readline()
        qual = IN.readline()
        if key.strip()==key.split()[0] and tag.strip() == "+":
            readCount += 1
            if readCount == 100:
                return False
        else:
            return True
        key = IN.readline()
    
    return False
    IN.close()

def simplify_fq(inFQFn, outFQFn):
    IN = IOFile.IOFile(inFQFn, 'r')
    OUT = IOFile.IOFile(outFQFn, 'w')
    
    key = IN.readline()
    while key:
        seq = IN.readline()
        tag = IN.readline()
        qual = IN.readline()
        key = key.split()[0] + "\n"
        tag = "+\n"
        OUT.writelines(key+seq+tag+qual)
        key = IN.readline()
    
    IN.close()
    OUT.close()


def main():
    params = init()
    inFq = params['inFastq']
    simplified = False
    if params['simplify'] and need_to_simplify(inFq):
        simplified = True
        rid = random.randint(10000,99999)
        tmp_file = params['outFastq'] + '.' + str(rid) + '.fastq'
        print("Start to simplify fastq...")
        simplify_fq(params['inFastq'], tmp_file)
        inFq = tmp_file
    
    CMD = "java -mx256m -jar %s SE -threads %s -phred33 %s %s HEADCROP:%s ILLUMINACLIP:%s:2:30:4 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:%s" % (trimmomatic, params['threads'], inFq, params['outFastq'], params['leading'], params['adator'], params['minLen'])
    print("Start to run trimmomatic:\n\t%s" % (CMD, ))
    os.system(CMD)
    if simplified:
        os.system("rm "+tmp_file)

if __name__ == "__main__":
    main()

