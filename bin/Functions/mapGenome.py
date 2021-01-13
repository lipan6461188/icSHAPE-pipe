#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import os
import version

Usage = """
mapGenome - Map reads to genome with STAR
=============================================================
\x1b[1mUSAGE:\x1b[0m 
  %s [-p 1 --maxMMap 1 --maxMisMatch 2 --noMut5] -i inFastq -o outprefix -x index
\x1b[1mHELP:\x1b[0m
  -i                    <String>
                            Input a fastq file
  -o                    <String>
                            Prefix (with full path) of all outputfiles (like D1, D2...)
                            Final bam file will be Prefix.sorted.bam
  -x                    <String>
                            Input a index to map

  More options:

  --tool                <STAR/hisat2>
                            Which tool to be used for mapping (default: STAR)

  [option for STAR]
  --maxMMap             <Int>
                            Maximun multiple map to be allowed (default: 1, unique map)
  --maxMisMatch         <Int>
                            Maximun mismatch to be allowed (default: 2)
  --alignMode           <Local/EndToEnd>
                            Mapping the reads to genome with end-to-end mode or local mode (default: EndToEnd)
  --noWithin            <None>
                            Unmapped reads not within the sam file (default: within)
  
  [option for Both]
  --maxReport           <Int>
                            Maximum number of alignment to report (default: 1)
                            For hisat2: May not be unique map when multi-mapped
  -p                    <Int>
                            How many threads to use (default: 1)
  --noMut5              <None>
                            Remove reads with mutation at the first base in 5' (such as MD:Z:0A)
                            The removed reads will not appear in the sorted.bam file
  --moreparams          <String>
                            More parameters. For example: "--outReadsUnmapped\\ Fastx"


\x1b[1mVERSION:\x1b[0m
    %s

\x1b[1mAUTHOR:\x1b[0m
    Li Pan

""" % (sys.argv[0], version.Version)

def init():
    import getopt
    
    Params = { 'inFastq': None, 'outPrefix': None, 'index': None, 'threads': 1, 
                'maxMMap': 1, 'maxMisMatch': 2, 'noMut5':False, 'alignMode':'EndToEnd', 
                'noWithin': False, 'tool': 'STAR', 'maxReport':1, 'moreparams':'' }
    
    opts, args = getopt.getopt(sys.argv[1:], 'hi:o:x:p:', 
        ['maxMMap=', 'maxMisMatch=', 'noMut5', 'alignMode=', 'noWithin', 'maxReport=', 'tool=', 'moreparams='])
    
    for op, value in opts:
        if op == '-h':
            print(Usage)
            exit(-1)
        # Basic Parameters
        elif op == '-i':
            Params['inFastq'] = os.path.abspath(value)
        elif op == '-o':
            Params['outPrefix'] = os.path.abspath(value)
        elif op == '-x':
            Params['index'] = os.path.abspath(value)
        elif op == '-p':
            Params['threads'] = int(value)
        
        elif op == '--maxMMap':
            Params['maxMMap'] = int(value)
        elif op == '--maxMisMatch':
            Params['maxMisMatch'] = int(value)
        elif op == '--noMut5':
            Params['noMut5'] = True
        elif op == '--alignMode':
            assert value in ("Local", "EndToEnd")
            Params['alignMode'] = value
        elif op == '--noWithin':
            Params['noWithin'] = True
        
        elif op == '--tool':
            if value not in ("STAR", "hisat2"):
                sys.stderr.writelines("Error: --tool must be one of STAR/hisat2. You provide {}\n".format(value))
                exit(-1)
            Params['tool'] = value
        
        elif op == '--maxReport':
            Params['maxReport'] = int(value)

        elif op == '--moreparams':
            Params['moreparams'] = value.strip()
        
        else:
            sys.stderr.writelines("parameter Error: unrecognized parameter: "+op+"\n")
            print(Usage)
            sys.exit(-1)
    
    if not Params['inFastq']:
        sys.stderr.writelines("Error: please specify -i"+"\n")
        print(Usage)
        exit(-1)
    if not Params['outPrefix']:
        sys.stderr.writelines("Error: please specify -o"+"\n")
        print(Usage)
        exit(-1)
    if not Params['index']:
        sys.stderr.writelines("Error: please specify -x"+"\n")
        print(Usage)
        exit(-1)
    
    return Params


def build_STAR_cmd(params):
    unsorted_bam = params['outPrefix'] + ".unsorted.bam"
    
    CMD_1 = "STAR --readFilesIn %s \
        --outFileNamePrefix %s. \
        --genomeDir %s \
        --runThreadN %s \
        --genomeLoad NoSharedMemory \
        --runMode alignReads \
        --outSAMtype BAM Unsorted \
        --outSAMmultNmax %s \
        --outFilterMultimapNmax %s \
        --outFilterMismatchNmax %s \
        --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
        --outSAMstrandField intronMotif \
        --outSJfilterOverhangMin 30 12 12 12 \
        --alignEndsType %s \
        --outSAMattributes All \
        --outSAMunmapped %s \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJDBoverhangMin 1 \
        --outStd BAM_Unsorted"
    
    if params['inFastq'].endswith(".gz"): 
        CMD_1  += " --readFilesCommand zcat"
    
    if params['moreparams']:
        CMD_1  += " " + params['moreparams']
    
    if params['noWithin']: within = "None"
    else: within = "Within"
    
    final_cmd = CMD_1 % (params['inFastq'], params['outPrefix'], params['index'], params['threads'], params['maxReport'], params['maxMMap'], params['maxMisMatch'], params['alignMode'], within) + " > " + unsorted_bam
    return final_cmd

def build_hisat2_cmd(params):
    unsorted_bam = params['outPrefix'] + ".unsorted.bam"
    summary_file = params['outPrefix'] + ".summary"
    
    CMD_1 = "hisat2 \
        -U %s \
        --summary-file %s \
        --quiet \
        -k %s \
        -x %s \
        --reorder \
        --rna-strandness F \
        -p %s " + params['moreparams'] + " | samtools view --threads %s -bh -o %s -" % (params['inFastq'], summary_file, params['maxReport'], params['index'], params['threads'], params['threads'], unsorted_bam)
    
    return CMD_1

def main():
        
    CMD_sort_1 = "samtools sort -m 2G --threads %s %s -o %s"
    CMD_sort_2 = r"""samtools view -h %s | awk '$0~/^@/{print $0}$0!~/^@/{for(i=12;i<NF;i++){if(substr($i,1,4)=="MD:Z"){if(and(16,$2)==0){ if( $i!~/^MD:Z:0/ ) print $0; }else{ if($i!~/^MD:Z:.*0$/) print $0; }}}}' | samtools view --threads %s -bh - | samtools sort -m 2G --threads %s -o %s -"""
    
    params = init()
    unsorted_bam = params['outPrefix'] + ".unsorted.bam"
    sorted_bam = params['outPrefix'] + ".sorted.bam"
    
    if params['tool'] == 'STAR':
        CMD_1 = build_STAR_cmd(params)
    else:
        CMD_1 = build_hisat2_cmd(params)

    CMD_sort_1 = CMD_sort_1 % (params['threads'], unsorted_bam, sorted_bam)
    CMD_sort_2 = CMD_sort_2 % (unsorted_bam, params['threads'], params['threads'], sorted_bam)
    
    print("Start to map to genome:\n\t%s" % (CMD_1, ))
    os.system(CMD_1)
    if params['noMut5']:
        print("Start to sort bam:\n\t%s" % (CMD_sort_2, ))
        os.system(CMD_sort_2)
    else:
        print("Start to sort bam:\n\t%s" % (CMD_sort_1, ))
        os.system(CMD_sort_1)

if __name__ == "__main__":
    main()

