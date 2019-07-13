#!/usr/bin/env python
#-*- coding:utf-8 -*- 

import re, sys, os, getopt, time, datetime
import GTFParserFunc

file_information = """
GTF: 
    1. Emsembl: https://asia.ensembl.org/info/data/ftp/index.html
    2. Gencode: https://www.gencodegenes.org/

GFF3:
    1. ftp://ftp.ncbi.nlm.nih.gov: ftp://ftp.ncbi.nlm.nih.gov/genomes/Equus_caballus/GFF/

"""

Usage = """
parseGTF -  Parse genome coordinate information and transcript
            coordinate information from Gencode GTF and NCBI GFF3 files
=============================================================
\x1b[1mUSAGE:\x1b[0m 
  %s -g GFF3/GTF -o output_prefix -s [ensembl|gencode|ncbi]
\x1b[1mHELP:\x1b[0m
     -g             Genome annotation
     -o             Output file prefix
     -s             <ensembl/gencode/ncbi> data source

     --genome       Fetch transcriptome from genome file, produce a prefix_transcriptome.fa file
     --noscaffold   Remove scaffolds. Scaffolds are defined as those chromosomes with id length > 6 and not startswith chr and NC_
     --rawchr       Use raw chromosome ID, don't convert NC_* to chrXX. Only useful for NCBI source GFF3

\x1b[1mVERSION:\x1b[0m
    %s

\x1b[1mAUTHOR:\x1b[0m
    Li Pan

""" % (sys.argv[0], "Test Version")



def init():
    params = { 'genomeFile': None, 'noscaffold': False, 'rawchr': False }
    opts, args = getopt.getopt(sys.argv[1:], 'hg:o:s:', ['genome=', 'noscaffold', 'rawchr'])
    for op, value in opts:
        if op == '-h':
            print(Usage)
            sys.exit(-1)
        elif op == '--genome':
            params['genomeFile'] = value
        elif op == '--noscaffold':
            params['noscaffold'] = True
        elif op == '--rawchr':
            params['rawchr'] = True
        elif op == '-g':
            params['input'] = os.path.abspath(value)
        elif op == '-o':
            params['output'] = value
        elif op == '-s':
            params['source'] = value
            assert value in ('ensembl', 'ncbi', 'gencode'), 'Error: -t <ensembl/gencode/ncbi>'
    if 'input' not in params:
        print('Error: specify -g')
        print(Usage)
        sys.exit(-1)
    if 'output' not in params:
        print('Error: specify -o')
        print(Usage)
        sys.exit(-1)
    if 'source' not in params:
        print('Error: specify -s')
        print(Usage)
        sys.exit(-1)
    return params

def main():
    params = init()
    genomeCoorFn = params['output']+'.genomeCoor.bed'
    transCoorFn = params['output']+'.transCoor.bed'
    transcriptomeFn = params['output']+'_transcriptome.fa'
    if params['source'] in ('gencode', 'ensembl'):
        gtf_container = GTFParserFunc.read_ensembl_gtf(params['input'], rem_scaffold=params['noscaffold'])
        GTFParserFunc.write_gtf_genomeCoor_bed(gtf_container, genomeCoorFn)
    elif params['source'] == 'ncbi':
        gff3_container = GTFParserFunc.read_ncbi_gff3(params['input'], rem_scaffold=params['noscaffold'], raw_chrID=params['rawchr'])
        GTFParserFunc.write_gff3_genomeCoor_bed(gff3_container, genomeCoorFn, raw_chrID=params['rawchr'])
    
    GTFParserFunc.genomeCoorBed_To_transCoorBed(genomeCoorFn, transCoorFn)
    if params['genomeFile']:
        GTFParserFunc.writeTranscriptome(genomeCoorFn, params['genomeFile'], transcriptomeFn, verbose=True, showAttr=False)

if __name__ == '__main__':
    main()

