#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import os
import gzip
import version

Usage = """
readCollapse - Remove duplicate reads from fastq
================================================
\x1b[1mUSAGE:\x1b[0m
  %s -1 fastq_PE_reads_1 -2 fastq_PE_reads_2 -U fastq_SE_reads
\x1b[1mHELP:\x1b[0m
  -U                Single ends read
  -1                Paired ends read 1
  -2                Paired ends read 2

 More options
  -o                Single ends read output 
  -p                PE read output 1
  -q                PE read output 2

  -f                Unique fasta after collapse
  -l                <1-10> barcode length (default: guess)
                    Longer barcode, less memory, more time

  --simplify        Simplify fastq information and save space:
                    Example: 
                        line 1n:"@SRR1057939.1 other information" => "@SRR1057939.1"
                        line 3n:"+ other information" => "+"

\x1b[1mVERSION:\x1b[0m
    %s

\x1b[1mAUTHOR:\x1b[0m
    Li Pan

""" % (sys.argv[0], version.Version)

dirname = os.path.dirname(os.path.abspath(__file__))
SplitBarcode = "python " + os.path.join(dirname, 'splitByBarcode.py')
collapseFq = "python " + os.path.join(dirname, 'collapseSingleFq.py')

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
    
    Params = { 'outFasta': None, 'bcLen': None, 'simplify': False }
    
    opts, args = getopt.getopt(sys.argv[1:], 'hU:1:2:o:p:q:f:l:', ['simplify'])
    
    for op, value in opts:
        if op == '-h':
            print(Usage)
            exit(-1)
        # Basic Parameters
        elif op == '-U':
            Params['inFile'] = os.path.abspath(value)
        elif op == '-o':
            Params['outFile'] = os.path.abspath(value)
        
        elif op == '-1':
            Params['inFile_1'] = os.path.abspath(value)
        elif op == '-2':
            Params['inFile_2'] = os.path.abspath(value)
        elif op == '-p':
            Params['outFile_1'] = os.path.abspath(value)
        elif op == '-q':
            Params['outFile_2'] = os.path.abspath(value)
        
        elif op == '-l':
            Params['bcLen'] = int(value)
            assert 1 <= Params['bcLen'] <= 10
        
        elif op == '-f':
            Params['outFasta'] = os.path.abspath(value)
        elif op == '--simplify':
            Params['simplify'] = True
        
        else:
            sys.stderr.writelines("parameter Error: unrecognized parameter: "+op+"\n")
            print(Usage)
            sys.exit(-1)
    
    if len(sys.argv) == 1:
        print(Usage)
        sys.exit(-1)
    
    if 'inFile' in Params:
        Params['mode'] = 'single'
        sys.stderr.writelines("Single-end mode..."+"\n")
        if 'outFile' not in Params:
            sys.stderr.writelines("Error: \n\tSigle-end mode: specify -U, -o\n\tPair-end mode: specify -1, -2, -p and -q"+"\n")
            exit(-1)
        outFile, outDir, fileSuffix = fileparse(Params['outFile'])
        Params['outDir'] = outDir
    else:
        Params['mode'] = 'pair'
        sys.stderr.writelines("Pair-end mode..."+"\n")
        if ('inFile_1' not in Params) or ('inFile_2' not in Params) or ('outFile_1' not in Params) or ('outFile_2' not in Params):
            sys.stderr.writelines("Error: \n\tSigle-end mode: specify -U, -o\n\tPair-end mode: specify -1, -2, -p and -q"+"\n")
            exit(-1)
        outFile, outDir, fileSuffix = fileparse(Params['outFile_1'])
        Params['outDir'] = outDir
    
    return Params

def main():
    Params = init()
    
    tmpDir = Params['outDir'] + "/tmp_"+str(os.getpid())
    if not os.path.exists(tmpDir):
        sys.stderr.writelines(getoutput("mkdir "+tmpDir)+"\n")
    
    if Params['mode'] == 'single':
        sys.stderr.writelines("Collapsing file %s...\n\t" % (Params['inFile'], )+"\n")
        if os.path.exists(Params['outFile']):
            sys.stderr.writelines("Warning! %s exisits, will be overwritten. \n\t" % (Params['outFile'], )+"\n")
            sys.stderr.writelines(getoutput("rm "+Params['outFile'])+"\n")
        inFile = Params['inFile']
        outFile = Params['outFile']
    else:
        if os.path.exists(Params['outFile_1']):
            sys.stderr.writelines("Warning! %s exisits, will be overwritten. \n\t" % (Params['outFile_1'], )+"\n")
            sys.stderr.writelines(getoutput("rm "+Params['outFile_1'])+"\n")
        if os.path.exists(Params['outFile_2']):
            sys.stderr.writelines("Warning! %s exisits, will be overwritten. \n\t" % (Params['outFile_2'], )+"\n")
            sys.stderr.writelines(getoutput("rm "+Params['outFile_2'])+"\n")
        inFile = tmpDir + "/input.fastq"
        outFile = tmpDir + "/tmpOut.fastq"
        mergePairEndReads(Params['inFile_1'], Params['inFile_2'], inFile)
    
    if Params['outFasta'] and os.path.exists(Params['outFasta']):
        sys.stderr.writelines("Warning! %s exisits, will be overwritten. \n\t" % (Params['outFasta'], )+"\n")
        sys.stderr.writelines(getoutput("rm "+Params['outFasta'])+"\n")
    
    file2Collapse = [ inFile ]
    if Params['bcLen'] != None:
        cPos = 5
        cLen = Params['bcLen']
    else:
        cPos, cLen = estimateSplit(inFile)
    
    #print cPos, cLen
    if cLen>0:
        sys.stderr.writelines("File %s too large, will be splitted..." % ( inFile,  )+"\n")
        file2Collapse.pop()
        tmpOutDir = tmpDir
        CMD = SplitBarcode+" -i %s -o %s -p %s -l %s --mode new" % (inFile, tmpOutDir, cPos, cLen)
        code = os.system(CMD)
        if code != 0:
            sys.stderr.writelines("Error! splitting file failed. quiting..."+"\n")
            exit(-1)
        else:
            file2Collapse = getFqFiles( tmpOutDir )
            #print file2Collapse
    
    totalReads = 0
    uniqReads = 0
    outputFastas = []
    for file in file2Collapse:
        collapseResults = ""
        if Params['outFasta'] != None:
            CMD = collapseFq+" -i %s -o %s --mode append --fasta %s" % (file, outFile, Params['outFasta'])
        else:
            CMD = collapseFq+" -i %s -o %s --mode append" % (file, outFile)
        
        if Params['simplify']:
            CMD += " --simplify"
        #print CMD
        code = os.system(CMD)
        #volReads, volUniqReads, volUniqRatio, volFasta = _parseCollapseOutput ( "......." );
        #totalReads += volReads
        #uniqReads += volUniqReads
        #outputFastas.append( volFasta )
    
    #uniqRatio = ".2f" % (1.0*uniqReads/totalReads)
    if Params['mode'] == 'pair':
        splitPairEndReads ( outFile, Params['outFile_1'], Params['outFile_2'] )
        sys.stderr.writelines(getoutput("rm -f "+inFile)+"\n")
        sys.stderr.writelines(getoutput("rm -f "+outFile)+"\n")
        sys.stderr.writelines("Collapsing file %s and %s finished." % (Params['inFile_1'], Params['inFile_2'])+"\n")
    else:
        sys.stderr.writelines("Collapsing file %s finished." % (Params['inFile'], )+"\n")
    
    #print >>sys.stderr, "Read collapse successful! Total count: %s, unique count: %s, unique ratio: %s." % (totalReads, uniqReads, uniqRatio)
    os.system("rm -r "+tmpDir)

def mergePairEndReads(readFile1, readFile2, peFile):
    ## should test whether they are of the same length
    sys.stderr.writelines("merge two PE fastq files..."+"\n")
    inFile1 = readFile1
    inFile2 = readFile2
    if readFile1.endswith('.gz'):
        inFile1 = readFile1 + ".tmp.fq"
        sys.stderr.writelines(getoutput("gzip -d -c %s > %s" % (readFile1, inFile1))+"\n")
    elif readFile2.endswith('.gz'):
        inFile2 = readFile2 + ".tmp.fq"
        sys.stderr.writelines(getoutput("gzip -d -c %s > %s" % (readFile2, inFile2))+"\n")
    
    sys.stderr.writelines(getoutput("paste %s %s > %s" % (inFile1, inFile2, peFile))+"\n")

def estimateSplit(inputFile):
    import math
    
    fileSize = os.path.getsize(inputFile)
    if inputFile.endswith(".gz"):
        fileSize *= 5
    
    pos = 1
    length = int( math.log(fileSize/1000000000.0, 4.0) )
    if inputFile.endswith('.gz'):
        IN = gzip.open(inputFile, 'rt')
    else:
        IN = open(inputFile, 'r')
    line = IN.readline; line = IN.readline()
    IN.close()
    readLen = len(line)
    length = min(length, 4)
    if readLen <= 10:
        length = 0
    else:
        pos = int( readLen/4.0 + 10 )
        while pos + length > readLen:
            pos -= length
            if pos < 1:
                pos = 1
                break
    
    return pos, length

def splitPairEndReads(peFile, readFile1, readFile2):
    sys.stderr.writelines("split into two PE fastq files..."+"\n")
    PE = gzip.open(peFile, 'rt') if peFile.endswith('.gz') else open(peFile)
    R1 = gzip.open(readFile1, 'wt') if readFile1.endswith('.gz') else open(readFile1, 'w')
    R2 = gzip.open(readFile2, 'wt') if readFile2.endswith('.gz') else open(readFile2, 'w')
    for line in PE:
        r1, r2 = line.split('\t')
        R1.writelines(r1)
        R2.writelines(r2)
    PE.close()
    R1.close()
    R2.close()

def getFqFiles(Dir):
    files = os.listdir(Dir)
    fastqFiles = [ Dir.rstrip('/')+"/"+file for file in files if file.endswith(".fastq") or file.endswith(".fq") ]
    return fastqFiles

def fileparse(fillFilePath):
    fileDir = '/'.join(fillFilePath.split('/')[:-1])
    fileName = fillFilePath.split('/')[-1]
    pureFileName = '.'.join(fileName.split('.')[:-1])
    fileSuffix = fileName.split('.')[-1]
    return pureFileName, fileDir, fileSuffix

if __name__ == '__main__':
    main()

