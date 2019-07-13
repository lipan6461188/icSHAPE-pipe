#-*- coding:utf-8 -*-

""" This is a function library to covert coordination system


                                                                                 Genome
                                   geneCoor2genomeCoor(Parser, GeneID, Pos)  -           -  genomeCoor2transCoor(TransBin, Parser, Chr, Start, End, Strand)
         genomeCoor2geneCoor(GeneBin, Parser, Chr, Start, End, Strand)     -               -  transCoor2genomeCoor(Parser, TransID, Pos)
                                                                         -                   - 
                                                                       -                       -
                                                                     -                           -
                                                                   -                               -
                                                                 Gene ------------------------- Transcript    
                                                        geneCoor2transCoor(TransBin, Parser, GeneID, Start, End)
                                                        transCoor2geneCoor(GeneBin, Parser, TransID, Start, End)


* generate a GeneBin: 
    binize_Gene(Parser, bw=100000)

* generate a TransBin: 
    binize_Trans(Parser, bw=100000)

* helper functions:
    - genomeRange2TransCoor(TransBin, Parser, Chr, Start, End, Strand, bw=100000)
    - overlapTrans(transID, exonArray, strand, absStart, absEnd)
    - genomeRange2geneCoor(GeneBin, Parser, Chr, Start, End, Strand, bw=100000)
    - overlapGene(geneParser, geneID, strand, absStart, absEnd)

"""

import sys

class out_of_range(Exception):
    def __init__(self, message=""):
        super(out_of_range, self).__init__(message)


#### Common Functions


def exonstr2exonlist(exon_str):
    exon_items = exon_str.split(',')
    exon_pair = [ [int(it.split('-')[0]), int(it.split('-')[1])] for it in exon_items]
    return exon_pair


#### genome2Trans

def binize_Trans(Parser, bw=100000):
    "把基因组上的区域分段，加速查找过程。bw越小，占内存越高，越快"
    def statisticPsuedoChrSize(Parser):
        ChrLen = {}
        for trans_id in Parser.GAPer:
            RNA = Parser.GAPer[trans_id]
            Chr = RNA['chr']
            End = int(RNA['end'])
            try:
                ChrLen[Chr] = End if End > ChrLen[Chr] else ChrLen[Chr]
            except KeyError:
                ChrLen[Chr] = End
        return ChrLen
    chr_size = statisticPsuedoChrSize(Parser)
    Bin = {}
    for Chr in chr_size.keys():
        Bin[Chr] = {}; Bin[Chr]['+'] = {}; Bin[Chr]['-'] = {}
        count = int( chr_size[Chr]/bw ) + 1
        for idx in range(count):
            Bin[Chr]['+'][idx] = []; Bin[Chr]['-'][idx] = []
    for trans_id in Parser.GAPer:
        RNA = Parser.GAPer[trans_id]
        Chr = RNA['chr']
        Strand = RNA['strand']
        Start = int( int(RNA['start']) / bw )
        End = int( int(RNA['end']) / bw )
        for idx in range(Start, End+1):
            Bin[Chr][Strand][idx].append( trans_id )
    return Bin

def genomeRange2TransCoor(TransBin, Parser, Chr, Start, End, Strand, bw=100000):
    "左闭右开，0-based"
    idxBin = int(Start/bw)
    overlapRegions = []
    checkedTrans = []
    while idxBin <= int(End/bw):
        if idxBin > list(TransBin[Chr][Strand].keys())[-1]:
            break
        for trans_id in TransBin[Chr][Strand][idxBin]:
            if trans_id not in checkedTrans:
                checkedTrans.append(trans_id)
            else:
                continue
            RNA = Parser.GAPer[trans_id]
            if End < int(RNA['start']) or Start > int(RNA['end']):
                continue
            exonArray = exonstr2exonlist(RNA['exon_str'])
            overlapRegion = overlapTrans(trans_id, exonArray, Strand, Start, End)
            for item in overlapRegion: item.insert(0, Chr)
            if len(overlapRegion) > 0:
                overlapRegions += overlapRegion
        idxBin += 1
    return overlapRegions

def overlapTrans(transID, exonArray, strand, absStart, absEnd):
    "查找一个基因组区段与某个转录本的交集"
    absStart += 1
    numExon = len( exonArray )
    relExonStart = 0; startPosInExon = 0; endPosInExon = 0; absStartInExon = 0; absEndInExon = 0;
    overlapRegions = []
    for idxExon in range(numExon):
        exon = exonArray[idxExon]
        exonStart = int(exon[0])
        exonEnd = int(exon[1])
        exonLength = exonEnd - exonStart + 1
        if exonStart <= absEnd and exonEnd >= absStart:
            "有交集"
            if strand == '+':
                startPosInExon = absStart - exonStart + 1
                if startPosInExon < 1: startPosInExon = 1
                endPosInExon = absEnd - exonStart + 1
                if endPosInExon > exonLength: endPosInExon = exonLength
            elif strand == '-':
                startPosInExon = exonEnd - absEnd + 1
                if startPosInExon < 1: startPosInExon = 1
                endPosInExon = exonEnd - absStart + 1
                if endPosInExon > exonLength: endPosInExon = exonLength
            relStart = relExonStart + startPosInExon
            relEnd = relExonStart + endPosInExon
            absStartInExon = absStart if absStart >= exonStart else exonStart
            absEndInExon = exonEnd if absEnd >= exonEnd else absEnd
            "add to overlapRegions"
            overlapRegions.append( [absStartInExon, absEndInExon, transID, relStart, relEnd] )
        relExonStart += exonLength
    return overlapRegions

def genomeCoor2transCoor(TransBin, Parser, Chr, Start, End, Strand, pureTransID=False, bw=100000):
    # covert 1-based to 0-based
    assert(Start<=End)
    if Start <= 0: raise out_of_range()
    raw_transCoorList = genomeRange2TransCoor(TransBin, Parser, Chr=Chr, Start=Start-1, End=End, Strand=Strand, bw=bw)
    clean_transCoorList = []
    for (Chr, chr_start, chr_end, rnaID, trans_start, trans_end) in raw_transCoorList:
        if pureTransID: rnaID = rnaID.strip().split('.')[0]
        if chr_start <= chr_end:
            clean_transCoorList.append( [Chr, chr_start, chr_end, rnaID, trans_start, trans_end] )
    return clean_transCoorList

#### trans2Genome

def transCoor2genomeCoor(Parser, TransID, Pos):
    if Pos <= 0: raise out_of_range()
    RNA = Parser.GAPer[TransID]
    if Pos > Parser.getTransFeature(TransID)['trans_len']: raise out_of_range()
    exon_list = exonstr2exonlist(RNA['exon_str'])
    Strand = RNA['strand']
    accu = 0
    for i in range(len(exon_list)):
        range_of_tuple = abs( exon_list[i][1] - exon_list[i][0] ) + 1
        if accu + range_of_tuple + 1 > Pos:
            if Strand == '+':
                return [ RNA['chr'], exon_list[i][0]+Pos-accu-1, Strand ]
            else:
                return [ RNA['chr'], exon_list[i][1]-(Pos-accu-1), Strand ]
        accu += range_of_tuple
    sys.stderr.writelines("Warning: %s %s cannot be found genome Coordinate\n" % (TransID, Pos))
    return [False, False, False]


#### gene2Genome

def geneCoor2genomeCoor(Parser, GeneID, Pos):
    if Pos <= 0: raise out_of_range()
    geneInfo = Parser.getGeneParser(showAttr=False)
    Gene = geneInfo[GeneID]
    (Chr, Chr_Start, Chr_End, Strand, GeneName) = (Gene['chr'], Gene['start'], Gene['end'], Gene['strand'], Gene['gene_name'])
    geneLen = Chr_End - Chr_Start + 1
    if Pos > geneLen: raise out_of_range()
    if Strand == '+':
        return [Chr, Chr_Start+Pos-1, Strand]
    if Strand == '-':
        return [Chr, Chr_End-Pos+1, Strand]


#### genome2Gene

def binize_Gene(Parser, bw=100000):
    "把基因组上的区域分段，加速查找过程。bw越小，占内存越高，越快"
    def statisticPsuedoChrSize(geneParser):
        ChrLen = {}
        for geneID in geneParser:
            Gene = geneParser[geneID]
            Chr = Gene['chr']
            End = int(Gene['end'])
            try:
                ChrLen[Chr] = End if End > ChrLen[Chr] else ChrLen[Chr]
            except KeyError:
                ChrLen[Chr] = End
        return ChrLen
    geneParser = Parser.getGeneParser(showAttr=False)
    chr_size = statisticPsuedoChrSize(geneParser)
    Bin = {}
    for Chr in chr_size.keys():
        Bin[Chr] = {}; Bin[Chr]['+'] = {}; Bin[Chr]['-'] = {}
        count = int( chr_size[Chr]/bw ) + 1
        for idx in range(count):
            Bin[Chr]['+'][idx] = []; Bin[Chr]['-'][idx] = []
    for geneID in geneParser:
        Chr = geneParser[geneID]['chr']
        Strand = geneParser[geneID]['strand']
        Start = int( int(geneParser[geneID]['start']) / bw )
        End = int( int(geneParser[geneID]['end']) / bw )
        for idx in range(Start, End+1):
            Bin[Chr][Strand][idx].append( geneID )
    return Bin

def genomeRange2geneCoor(GeneBin, Parser, Chr, Start, End, Strand, bw=100000):
    "左闭右开，0-based"
    idxBin = int(Start/bw)
    overlapRegions = []
    checkedGenes = []
    geneParser = Parser.getGeneParser(showAttr=False)
    while idxBin <= int(End/bw):
        if idxBin > list(GeneBin[Chr][Strand].keys())[-1]:
            break
        for geneID in GeneBin[Chr][Strand][idxBin]:
            if geneID not in checkedGenes:
                checkedGenes.append(geneID)
            else:
                continue
            Gene = geneParser[geneID]
            if End < int(Gene['start']) or Start > int(Gene['end']):
                continue
            overlapRegion = overlapGene(geneParser, geneID, Strand, Start, End)
            if overlapRegion != None:
                overlapRegion.insert(0, Chr)
                overlapRegions += [overlapRegion]
        idxBin += 1
    return overlapRegions


def overlapGene(geneParser, geneID, strand, absStart, absEnd):
    "查找一个基因组区段与基因的交集"
    absStart += 1
    Gene = geneParser[geneID]
    geneStart = int(Gene['start'])
    geneEnd = int(Gene['end'])
    geneLength = geneEnd - geneStart + 1
    if absStart <= geneEnd and geneStart <= absEnd:
        if strand == '+':
            startPosInGene = absStart - geneStart + 1
            if startPosInGene < 1: startPosInGene = 1
            endPosInGene = absEnd - geneStart + 1
            if endPosInGene > geneLength: endPosInGene = geneLength
        elif strand == '-':
            startPosInGene = geneEnd - absEnd + 1
            if startPosInGene < 1: startPosInGene = 1
            endPosInGene = geneEnd - absStart + 1
            if endPosInGene > geneLength: endPosInGene = geneLength
        absStartInExon = absStart if absStart >= geneStart else geneStart
        absEndInExon = geneEnd if absEnd >= geneEnd else absEnd
        return [absStartInExon, absEndInExon, geneID, startPosInGene, endPosInGene]
    return None


def genomeCoor2geneCoor(GeneBin, Parser, Chr, Start, End, Strand, pureGeneID=False, bw=100000):
    # covert 1-based to 0-based
    assert(Start<=End)
    if Start <= 0: raise out_of_range()
    raw_geneCoorList = genomeRange2geneCoor(GeneBin, Parser, Chr=Chr, Start=Start-1, End=End, Strand=Strand, bw=bw)
    clean_geneCoorList = []
    for (Chr, chr_start, chr_end, rnaID, trans_start, trans_end) in raw_geneCoorList:
        if pureGeneID: rnaID = rnaID.strip().split('.')[0]
        if chr_start <= chr_end:
            clean_geneCoorList.append( [Chr, chr_start, chr_end, rnaID, trans_start, trans_end] )
    return clean_geneCoorList



#### gene2Trans

def geneCoor2transCoor(TransBin, Parser, GeneID, Start, End):
    assert(Start<=End)
    if Start <= 0: raise out_of_range()
    (genome_chr, genome_pos_start, genome_strand) = geneCoor2genomeCoor(Parser, GeneID, Start)
    (genome_chr, genome_pos_end, genome_strand) = geneCoor2genomeCoor(Parser, GeneID, End)
    if genome_chr == False:
        return [ ]
    if genome_strand == '-':
        genome_pos_start, genome_pos_end = genome_pos_end, genome_pos_start
    return genomeCoor2transCoor(TransBin, Parser, genome_chr, genome_pos_start, genome_pos_end, genome_strand)

#### trans2Gene

def transCoor2geneCoor(GeneBin, Parser, TransID, Start, End):
    assert(Start<=End)
    if Start <= 0: raise out_of_range()
    (genome_chr, genome_pos_start, genome_strand) = transCoor2genomeCoor(Parser, TransID, Start)
    (genome_chr, genome_pos_end, genome_strand) = transCoor2genomeCoor(Parser, TransID, End)
    if genome_chr == False:
        return [ ]
    if genome_strand == '-':
        genome_pos_start, genome_pos_end = genome_pos_end, genome_pos_start
    return genomeCoor2geneCoor(GeneBin, Parser, genome_chr, genome_pos_start, genome_pos_end, genome_strand)







