#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys, re

Usage = """
combineGenome - Combine chromosomes from NCBI and convert the
                NC_ code to chrX code
=============================================================
\x1b[1mUSAGE:\x1b[0m 
  %s chr1.fasta chr2.fasta chr3.fasta... > outFile.fa

\x1b[1mDATE:\x1b[0m
    2018-12-23

\x1b[1mAUTHOR:\x1b[0m
    Li Pan

""" % (sys.argv[0], )


def main():
    if len(sys.argv)<=1 or sys.argv[0] in ("--help", "-h"):
        print(Usage)
        exit(-1)
    
    inFiles = sys.argv[1:]
    for infile in inFiles:
        sys.stderr.writelines("Process "+infile+"..."+"\n")
        for line in open(infile):
            if line[0] == '>':
                data = line[1:].rstrip().split('|')
                
                annotation = data[-1]
                chrMatch = re.findall("chromosome (\\w+)", annotation)
                
                NC_code = ""
                for i in range(len(data)):
                    if data[i] == "ref":
                        NC_code = data[i+1]
                
                if 'unplaced' in line or 'unlocalized' in line:
                    head = ">%s\n" % (NC_code, )
                    
                elif 'mitocho' in line:
                    data = line.split("|")
                    head = ">chrM %s\n" % (NC_code, )
                
                elif chrMatch:
                    chr_id = chrMatch[0]
                    head = ">chr%s %s\n" % (chr_id, NC_code)
                
                elif '|' in line:
                    head = ">%s\n" % (NC_code, )
                
                else:
                    pass
                
                sys.stdout.writelines(head)
            
            else:
                sys.stdout.writelines(line)


if __name__ == "__main__":
    main()


