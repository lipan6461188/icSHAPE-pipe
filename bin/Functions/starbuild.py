#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import os
import re
import numpy
import version

Usage = """
starbuild - Build a STAR index
=========================================================================================
\x1b[1mUSAGE:\x1b[0m 
  %s [--gtf GTF_File -p 1] -i genome.fa -o out_dir
\x1b[1mHELP:\x1b[0m
  -i                    <String>
                            Input a genome fasta file
  -o                    <String>
                            Output a path to save index

  More options:
  --gtf                 <String>
                            Provide a GTF file to build index
  -p                    <Int>
                            How many threads to use (default: 1)
  --noscaffold          <None>
                            Don't build index for scaffold, scaffolds are defined as 
                            those chromosomes with id length > 6 and not startswith chr and NC_

\x1b[1mVERSION:\x1b[0m
    %s

\x1b[1mAUTHOR:\x1b[0m
    Li Pan

""" % (sys.argv[0], version.Version)

def load_getoutput():
    import subprocess
    if 'getoutput' in dir(subprocess):
        return subprocess.getoutput
    else:
        import commands
        return commands.getoutput

getoutput = load_getoutput()

def init():
    import getopt
    
    Params = { 'inFile': None, 'outDir': None, 'gtfFile': None, 'threads': 1, 'noscaffold': False }
    
    opts, args = getopt.getopt(sys.argv[1:], 'hi:o:p:', ['gtf=', 'noscaffold'])
    
    for op, value in opts:
        if op == '-h':
            print(Usage)
            exit(-1)
        # Basic Parameters
        elif op == '-i':
            Params['inFile'] = os.path.abspath(value)
        elif op == '-o':
            Params['outDir'] = os.path.abspath(value.rstrip('/'))
        
        elif op == '--gtf':
            Params['gtfFile'] = os.path.abspath(value)
        elif op == '-p':
            Params['threads'] = int(value)
        elif op == '--noscaffold':
            Params['noscaffold'] = True
        
        else:
            sys.stderr.writelines("parameter Error: unrecognized parameter: "+op+"\n")
            print(Usage)
            sys.exit(-1)
    
    if Params['inFile'] == None:
        sys.stderr.writelines("Error: please specify -i"+"\n")
        print(Usage)
        exit(-1)
    if Params['outDir'] == None:
        sys.stderr.writelines("Error: please specify -o"+"\n")
        print(Usage)
        exit(-1)
    
    return Params

def build_noscaffold_fasta(input_fa, output_fa, verbose=True):
    OUT = open(output_fa, 'w')
    write_cur_trans = False
    removed_chr_list = []
    for line in open(input_fa):
        if line[0] == '>':
            chrID = line[1:].split()[0]
            if chrID.startswith('chr') or chrID.startswith('NC_') or len(chrID)<=6:
                write_cur_trans = True
                OUT.writelines(line)
                if verbose:
                    sys.stderr.writelines("Writing chromosome " + chrID + "..."+"\n")
            else:
                write_cur_trans = False
                removed_chr_list.append(chrID)
        elif write_cur_trans:
            OUT.writelines(line)
    
    if verbose:
        sys.stderr.writelines("Warning: Removed scaffolds: "+ "\t".join(removed_chr_list)+"\n")
    
    OUT.close()

def count_fasta(inFasta):
    ref_num = 0
    base_num = 0
    for line in open(inFasta):
        if line[0] == '>':
            ref_num += 1
        else:
            base_num += len(line) - 1
    return ref_num, base_num

def main():
    params = init()
    
    if params['noscaffold']:
        import random
        randID = random.randint(100000,999999)
        tmp_fa = "/tmp/tmp_genome_%s.fa" % (randID, )
        print("Start build temp genome: "+tmp_fa)
        build_noscaffold_fasta(params['inFile'], tmp_fa, verbose=True)
        params['inFile'] = tmp_fa
    
    CMD = "STAR --runMode genomeGenerate --genomeFastaFiles %s --genomeDir %s --runThreadN %s" % (params['inFile'], params['outDir'], params['threads'])
    if params['gtfFile']:
        CMD += " --sjdbGTFfile " + params['gtfFile']
    
    ref_num, base_num = count_fasta(params['inFile'])
    genomeSAindexNbases = min( 14, numpy.log2(base_num)/2-1 )
    if genomeSAindexNbases < 14:
        CMD += " --genomeSAindexNbases %.1f" % (genomeSAindexNbases, )
    
    if ref_num > 5000:
        genomeChrBinNbits = min( 18, numpy.log2(1.0*base_num/ref_num) )
        if genomeChrBinNbits < 18:
            CMD += " --genomeChrBinNbits %.1f" % (genomeChrBinNbits, )
    
    print("Start to build STAR index:\n\t%s" % (CMD, ))
    os.system(CMD)
    
    if params['noscaffold']:
        os.remove(tmp_fa)

if __name__ == "__main__":
    main()
