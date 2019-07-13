#!/usr/bin/env python
#-*- coding:utf-8 -*-

import os,sys

if len(sys.argv) < 3:
    print("Usage: {0} inputfasta.fa outputfasta.fa [maxlength]".format(sys.argv[0]))
    exit(0)

IN = open(sys.argv[1], 'r')
OUT = open(sys.argv[2], 'w')

maxLen = 200
if len(sys.argv) >= 4:
    maxLen = int(sys.argv[3])

headLine = IN.readline()
seq = ""
seqLen = 0
for line in IN:
    if line[0] == '>':
        if seqLen <= maxLen:
            OUT.writelines(headLine+seq)
        seq = ""
        headLine = line
        seqLen = 0
    else:
        seq += line
        seqLen += len(seq.strip())

if seqLen <= maxLen:
    OUT.writelines(headLine+seq)

OUT.close()
IN.close()


