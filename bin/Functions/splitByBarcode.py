#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import os
import gzip
import version

Usage = """
splitByBarcode - Split raw sequencing fastq file to multiple files with different barcode
=========================================================================================
\x1b[1mUSAGE:\x1b[0m 
  %s -i input_fastq -o output_dir -p barcode_pos -l barcode_len --mode mode --lib library --gzip
\x1b[1mHELP:\x1b[0m
    -p              <Int>
                        The start site of the barcode
    -l              <Int>
                        The length of the barcode
    --mode          <append/new>. 
                        [append] will lump reads into files in the outputdirectory with names like LIB_DMSO1.fastq, LIB_DMSO2.fastq...
                        [new] will create separate files for each input file, i.e., names like LIB_DMSO1_dataset1.fastq, LIB_DMSO2_dataset1.fastq...
    --lib           <String>
                        barcode1:libName1::barcode2:libName2::barcode3:libName3...
    --gzip          <None>
                        Output files will be compressed

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
    
    Params = { 'mode': 'new', 'library': None, 'gzip': False }
    
    opts, args = getopt.getopt(sys.argv[1:], 'hi:o:p:l:', ['mode=', 'lib=', 'gzip'])
    
    for op, value in opts:
        if op == '-h':
            print(Usage)
            exit(-1)
        # Basic Parameters
        elif op == '-i':
            Params['inFile'] = value
        elif op == '-o':
            Params['outDir'] = value
        elif op == '-p':
            Params['barPos'] = int(value)-1
        elif op == '-l':
            Params['barLen'] = int(value)
        elif op == '--mode':
            Params['mode'] = value
        elif op == '--lib':
            Params['library'] = value
        elif op == '--gzip':
            Params['gzip'] = True
        else:
            sys.stderr.writelines("parameter Error: unrecognized parameter: "+op+"\n")
            print(Usage)
            sys.exit(-1)
    
    if 'inFile' not in Params:
        sys.stderr.writelines("Error: specify -i"+"\n")
        print(Usage)
        exit(-1)
    if 'outDir' not in Params:
        sys.stderr.writelines("Error: specify -o"+"\n")
        print(Usage)
        exit(-1)
    if 'barPos' not in Params:
        sys.stderr.writelines("Error: specify -p"+"\n")
        print(Usage)
        exit(-1)
    if 'barLen' not in Params:
        sys.stderr.writelines("Error: specify -l"+"\n")
        print(Usage)
        exit(-1)

    return Params


def main():
    Params = init()
    
    sys.stderr.writelines("split fastq files %s ...\n\t" % (Params['inFile'], )+"\n")
    
    lib_info = {}
    if Params['library'] != None:
        lib_barcode = {}
        libs = Params['library'].split('::')
        for lib in libs:
            barcode, libName = lib.split(':')
            lib_barcode[barcode] = libName
        lib_info = lib_barcode
    else:
        lib_info = genCodeLib( Params['barLen'] );
    
    inFileName, inFileDir, inFileSuffix = fileparse(Params['inFile'])
    lib_outFile, bc_count = splitByLibrary(Params['inFile'], Params['outDir'], inFileName, \
        inFileSuffix, lib_info, Params['barPos'], Params['barLen'], Params['library'], Params['mode'], Params['gzip'] );
    
    statFile = Params['outDir'] + "/splitFastq.stat";
    outputStat ( statFile, bc_count, lib_info, Params['inFile'], Params['library'], Params['mode'] )    
    
    sys.stderr.writelines("File %s successfully splited into:" % (Params['inFile'], )+"\n")
    for bc in lib_outFile:
        if os.path.exists(lib_outFile[bc]) and lib_outFile[bc] != "":
            sys.stderr.writelines(" "+lib_outFile[bc]+"\n")
    
    sys.stderr.writelines("Splitting finished.\n")
    
    return 0;

def splitByLibrary(inFastqFile, outDir, fileName, fileSuffix, lib_info, barPos, barLen, library=None, mode='new', compress=False):
    lib_outFile = {}
    prepareLibraryOutputFiles( outDir, fileName, fileSuffix, lib_outFile, lib_info, library, mode, compress )
    
    count = 0
    bc_count = {}
    lib_content = {}
    key1 = ""; seq1 = ""; tag1 = ""; quality1 = ""
    if inFastqFile.endswith('gz'):
        IN = gzip.open(inFastqFile, 'rt')
    else:
        IN = open(inFastqFile)
    key1 = IN.readline()
    while key1:
        seq1 = IN.readline()
        tag1 = IN.readline()
        quality1 = IN.readline()
        
        count += 1
        #print count
        if count % 100000 == 0:
            sys.stderr.writelines("  %s" % (count, )+"\n")
            writeCachedOutput ( lib_outFile, lib_content )
            cleanLibContent ( lib_content )
        
        bc = seq1[barPos:barPos+barLen]
        lib = lib_info[bc] if (bc in lib_info) else "unmatched"
        assignLib ( key1, seq1, tag1, quality1, bc, lib, lib_content, bc_count )
        
        key1 = IN.readline()
    IN.close()

    writeCachedOutput( lib_outFile, lib_content )
    return lib_outFile, bc_count

def writeCachedOutput(lib_outFile, lib_content):
    for lib in lib_content:
        outputFile = lib_outFile[lib]
        if outputFile.endswith('.gz'):
            OUT = gzip.open(outputFile, 'ab')
        else:
            OUT = open(outputFile, 'a')
        if lib in lib_content:
            for read in lib_content[lib]:
                OUT.writelines(read)
        OUT.close()

def assignLib(key, seq, tag, quality, bc, lib, lib_content, bc_count):
    bc_count[bc] = bc_count.get(bc, 0) + 1
    
    if lib not in lib_content:
        lib_content[lib] = [ "".join([key, seq, tag, quality]) ]
    else:
        lib_content[lib].append( "".join([key, seq, tag, quality]) )

def cleanLibContent(lib_content):
    lib_content['unmatched'] = []
    for lib in lib_content:
        lib_content[lib] = []

def prepareLibraryOutputFiles ( outDir, fileName, fileSuffix, lib_outFile, lib_info, library, mode, compress ):
    name = "unmatched"
    if mode == "new":
        name = ("unmatched_"+fileName) if (library != None) else (fileName+"_unmatched")
    
    lib_outFile['unmatched'] = outDir + "/" + name + "." + fileSuffix
    
    if mode == "new" and os.path.exists(lib_outFile['unmatched']):
        sys.stderr.writelines("Warning! %s exists...will be erased\n" % (lib_outFile['unmatched'], )+"\n")
        sys.stderr.writelines(getoutput("rm "+lib_outFile['unmatched'])+"\n")
    
    for barcode in lib_info:
        lib = lib_info[barcode]
        #name = lib
        if mode == 'new':
            name = (lib+'_'+fileName) if (library != None) else (fileName+'_'+lib)
            outFile = outDir + '/' + name + ".fastq"
            if mode == 'new' and os.path.exists(outFile):
                sys.stderr.writelines("Warning! %s exists...will be erased\n" % (outFile, )+"\n")
                sys.stderr.writelines(getoutput("rm "+outFile)+"\n")
            lib_outFile[lib] = outFile
    
    if compress:
        for lib in lib_outFile:
            lib_outFile[lib] += '.gz'


def outputStat(statFile, bc_count, lib_info, inFile, library=None, mode='new'):
    if mode == 'new':
        if os.path.exists(statFile):
            sys.stderr.writelines("Warning! %s exists and will be overwitten.\n" % (statFile, )+"\n")
            sys.stderr.writelines(getoutput("rm "+statFile)+"\n")
    
    BCT = open(statFile, 'a')
    BCT.writelines( '# '+inFile+"\n-----------\n" )
    totalCount = 0
    for bc in bc_count.keys():
        totalCount += bc_count[bc]
    BCT.writelines( "total\t%d\n" % (totalCount, ) )
    
    sorted_bc = sorted(bc_count.keys(), key=lambda x: bc_count[x], reverse=True)
    for bc in sorted_bc:
        BCT.writelines( bc+"\t"+str(bc_count[bc]) )
        if (library != None) and (bc in lib_info):
            BCT.writelines("\t"+lib_info[bc])
        BCT.writelines("\n")
    BCT.close()


def fileparse(fillFilePath):
    fileDir = '/'.join(fillFilePath.split('/')[:-1])
    fileName = fillFilePath.split('/')[-1]
    pureFileName = '.'.join(fileName.split('.')[:-1])
    fileSuffix = fileName.split('.')[-1]
    return pureFileName, fileDir, fileSuffix


def genCodeLib(codeLength):
    alphabets = ['A', 'T', 'G', 'C']
    
    lib_code = { 'unmatched': 'unmatched' }
    
    # k-mer generator from michael eisen
    words = alphabets
    newwords = []
    for i in range(1, codeLength):
        newwords = []
        for w in words:
            for b in alphabets:
                newwords.append( w+b )
        words = newwords
    
    for word in words:
        lib_code[word] = word;
    
    return lib_code

if __name__ == "__main__":
    main()



