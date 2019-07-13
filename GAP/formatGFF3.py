#!/usr/bin/env python
#-*- coding:utf-8 -*- 

import re, sys, os, getopt, time, datetime
import GTFParserFunc

Usage = """
formatGFF3 -  Transform GFF3 NC_* code to chr* code
=============================================================
\x1b[1mUSAGE:\x1b[0m 
  %s inputGFF3 outputGFF3

\x1b[1mVERSION:\x1b[0m
    %s

\x1b[1mAUTHOR:\x1b[0m
    Li Pan

""" % (sys.argv[0], "Test Version")

def main():
    if len(sys.argv) != 3:
        print(Usage)
        exit(-1)
    
    reformat_GFF3(sys.argv[1], sys.argv[2])

def reformat_GFF3(inGFF3, outGFF3):
    
    pairser = GTFParserFunc.read_ncbi_gff3(inGFF3)
    NC_To_chr_dict = GTFParserFunc.build_NC_To_chr_dict(pairser)
    
    OUT = open(outGFF3, 'w')
    for line in open(inGFF3):
        if line[0] == '#':
            OUT.writelines(line)
        else:
            data = line.strip().split('\t')
            data[0] = NC_To_chr_dict.get(data[0], data[0])
            OUT.writelines("\t".join(data)+"\n")
    
    OUT.close()

if __name__ == '__main__':
    main()




