#-*- coding:utf-8 -*-

import GTFParserFunc
import CoorFunc
import sys
import GAP_Colors as Colors

___ParseTrans = {
    'first_biuldGeneParser': 1
}

def showBeautifulExample(exmaple, title):
    sys.stderr.writelines('# =-=-=-=-=-=-=-=-=-=-=-=-=-%s=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n' % (title, ))
    Pattern = "\t%15s:\t%-25s"
    for attr in sorted(exmaple.keys()):
        sys.stderr.writelines(Pattern % (attr, exmaple[attr])+"\n")
    sys.stderr.writelines('# =-=-=-=-=-=-=-=-=-=-=-=-=-%s=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n' % ( int(len(title)/2) * "=-", ))

def format_Exon_UTR_str(raw_str, strand):
    regions_list = [ (int(region.split('-')[0]), int(region.split('-')[1])) for region in raw_str.split(',') ]
    
    if strand == '+':
        regions_list.sort(key=lambda x: x[0])
    else:
        regions_list.sort(key=lambda x: x[0], reverse=True)
    
    return ",".join([str(region[0])+"-"+str(region[1]) for region in regions_list])

def parseUTR(utr_str, trans_len):
    """
        Parse UTR string
        -> without 5'UTR
            utr_5_start=utr_5_end=0
        -> without 3'UTR
            utr_3_start=utr_3_end=trans_len+1
    """
    utr_5_start = utr_5_end = 0 
    utr_3_start = utr_3_end = trans_len + 1
    cds_start = cds_end = 0
    if utr_str:
        utr_bak = utr_str.split(',')
        utr = []
        for utr_item in utr_bak:
            utr.append([int(i) for i in utr_item.split('-')])
        if utr[0][0] == 1:
            # have 5'UTR
            utr_5_start = 1
            for i in range(len(utr)-1):
                if utr[i][1] != utr[i+1][0] - 1:
                    utr_5_end = utr[i][1]
                    utr_3_start = utr[i+1][0]
                    utr_3_end = utr[-1][1]
                    break
            if utr_5_end == 0: 
                utr_5_end = utr[-1][1]    
        else:
            # have 5'UTR, have no 3'UTR
            utr_3_start = utr[0][0]
            utr_3_end = utr[-1][1]
    
    cds_start = utr_5_end + 1
    cds_end = utr_3_start - 1
    elem_coor = {'utr_5_start':utr_5_start, 'utr_5_end':utr_5_end, 'utr_3_start':utr_3_start, 'utr_3_end':utr_3_end, 'cds_start':cds_start, 'cds_end':cds_end}
    
    return elem_coor

def ParseGenomeCoorBedFileLine(GAPer, Line, remove_tid_version=False, remove_gid_version=False):
    arr = Line.strip().split('\t')
    (Chr, Start, End, Strand, GeneName_GeneID, TransID, GeneType, Exon) = arr[:8]
    
    Exon = format_Exon_UTR_str(Exon, Strand)
    (GeneName, GeneID) = (GeneName_GeneID.split('=')[0], GeneName_GeneID.split('=')[1])
    
    if remove_tid_version and '.' in TransID: TransID = ".".join(TransID.split('.')[:-1])
    if remove_gid_version and '.' in GeneID: GeneID = ".".join(GeneID.split('.')[:-1])
    
    exonList = GTFParserFunc.norm_exons(Exon)[0]
    TransLen = int(exonList[-1][1])
    utrString = ''
    if len(arr) == 9:
        arr[8] = format_Exon_UTR_str(arr[8], Strand)
        utrList = GTFParserFunc.norm_utr(Exon, arr[8], Strand)[1]
        utrString = ','.join([str(item[0])+'-'+str(item[1]) for item in utrList])
    
    utr_coor = parseUTR(utrString, TransLen)
    GAPer[TransID] = {}
    GAPer[TransID]['gene_id'] = GeneID
    GAPer[TransID]['gene_name'] = GeneName
    GAPer[TransID]['gene_type'] = GeneType
    GAPer[TransID]['trans_len'] = TransLen
    GAPer[TransID]['chr'] = Chr
    GAPer[TransID]['start'] = int(Start)
    GAPer[TransID]['end'] = int(End)
    GAPer[TransID]['strand'] = Strand
    GAPer[TransID]['exon_str'] = Exon
    for it in utr_coor:
        GAPer[TransID][it] = utr_coor[it]

def ParseGenomeCoorBedFile(genomeCoorBedFile, showAttr=True, remove_tid_version=False, remove_gid_version=False):
    """ 
        Get genome-annotation parser
    """
    GAPer = dict();
    
    for line in open(genomeCoorBedFile):
        ParseGenomeCoorBedFileLine(GAPer, line, remove_tid_version=remove_tid_version, remove_gid_version=remove_gid_version)
    
    # show attribuets
    if showAttr and len(GAPer) != 0:
        example_tid = list(GAPer.keys())[0]
        exmaple = GAPer[example_tid]
        showBeautifulExample(exmaple, title='Object Example (%s)' % (example_tid, ))
    
    return GAPer

def readSeq(seqFileName, remove_tid_version=False):
    Seq = {}
    cur_id = ""
    for line in open(seqFileName):
        if line[0] == '>':
            cur_id = line[1:].split()[0]
            if remove_tid_version and '.' in cur_id:
                cur_id = '.'.join(cur_id.split('.')[:-1])
            Seq[ cur_id ] = ''
        else:
            Seq[ cur_id ] += line.strip()
    
    return Seq

def biuldGeneParser(GAPer, showAttr=True):
    if ___ParseTrans['first_biuldGeneParser']:
        sys.stderr.writelines('Warning: Gene Coordinate System: [1_based_start, 1_based_end]\n')
        ___ParseTrans['first_biuldGeneParser'] = 0
    
    geneParser = {}
    for transID in GAPer:
        geneID = GAPer[transID]['gene_id']
        Chr = GAPer[transID]['chr']
        Start = GAPer[transID]['start']
        End = GAPer[transID]['end']
        Strand = GAPer[transID]['strand']
        geneName = GAPer[transID]['gene_name']
        geneType = GAPer[transID]['gene_type']
        try:
            geneParser[geneID]['transcript'] += [transID]
            if geneType not in geneParser[geneID]['gene_type']: geneParser[geneID]['gene_type'].append(geneType)
            if Start < geneParser[geneID]['start']: 
                geneParser[geneID]['start'] = Start
                geneParser[geneID]['length'] = geneParser[geneID]['end'] - geneParser[geneID]['start'] + 1
            if End > geneParser[geneID]['end']: 
                geneParser[geneID]['end'] = End
                geneParser[geneID]['length'] = geneParser[geneID]['end'] - geneParser[geneID]['start'] + 1
        except KeyError:
            geneParser[geneID] = {}
            geneParser[geneID]['transcript'] = [ transID ]
            geneParser[geneID]['chr'] = Chr
            geneParser[geneID]['start'] = Start
            geneParser[geneID]['end'] = End
            geneParser[geneID]['length'] = End - Start + 1
            geneParser[geneID]['strand'] = Strand
            geneParser[geneID]['gene_name'] = geneName
            geneParser[geneID]['gene_type'] = [ geneType ]
    
    # show attribuets
    if showAttr and len(geneParser) != 0:
        exmaple_gid = list(geneParser.keys())[0]
        exmaple = geneParser[exmaple_gid]
        showBeautifulExample(exmaple, title='Object Example (%s)' % (exmaple_gid, ))
    
    return geneParser

def GeneIntron(GAPer, geneParser, geneID):
    class NoExonError(Exception):
        pass
    def Intron_From_Exom(exon_str, start, end, strand):
        """ Get Intron Area From Exon String
        Example:
            exon_str = "489361-489710,485040-485208,476887-476945,476738-476882"
            start = 476738
            end = 489710
            Intron_From_Exom(exon_str, start, end, strand="-")
        """
        import copy
        genomeCoorIntron = []
        exonList = [ [int(exon_region.split('-')[0]), int(exon_region.split('-')[1])] for exon_region in exon_str.split(',') ]
        if len(exonList) == 0: raise NoExonError
        exonList.sort(key=lambda x: x[0], reverse=False)
        if start < exonList[0][0]:
            genomeCoorIntron.append( [start, exonList[0][0]-1] )
        for idx in range(len(exonList)-1):
            exon_region = exonList[idx]
            genomeCoorIntron.append( [exonList[idx][1]+1, exonList[idx+1][0]-1] )
        if end > exonList[-1][1]:
            genomeCoorIntron.append( [exonList[-1][1]+1, end] )
        GeneCoorIntron = copy.deepcopy(genomeCoorIntron)
        if strand == '+':
            for idx in range(len(GeneCoorIntron)):
                GeneCoorIntron[idx][0] = GeneCoorIntron[idx][0] - start + 1
                GeneCoorIntron[idx][1] = GeneCoorIntron[idx][1] - start + 1
        elif strand == '-':
            for idx in range(len(GeneCoorIntron)):
                GeneCoorIntron[idx][0] = end - GeneCoorIntron[idx][0] + 1
                GeneCoorIntron[idx][1] = end - GeneCoorIntron[idx][1] + 1
                GeneCoorIntron[idx][0], GeneCoorIntron[idx][1] = GeneCoorIntron[idx][1], GeneCoorIntron[idx][0]
        genomeCoorIntron.sort(key=lambda x: x[0], reverse=True if strand == '-' else False)
        GeneCoorIntron.sort(key=lambda x: x[0])
        return genomeCoorIntron, GeneCoorIntron
    gene_start = geneParser[geneID]['start']
    gene_end = geneParser[geneID]['end']
    strand = geneParser[geneID]['strand']
    IsoformIntrons = {}
    for transID in geneParser[geneID]['transcript']:
        exon_str = GAPer[transID]['exon_str']
        trans_start = GAPer[transID]['start']
        trans_end = GAPer[transID]['end']
        try:
            #IsoformIntrons[transID] = Intron_From_Exom(exon_str, gene_start, gene_end, strand)
            IsoformIntrons[transID] = Intron_From_Exom(exon_str, trans_start, trans_end, strand)
        except NoExonError:
            sys.stderr.writelines("Warning: %s have no Exon_str, Skip it\n" % (transID, ))
    return IsoformIntrons


def GeneExon(GAPer, geneParser, geneID):
    class NoExonError(Exception):
        pass
    def Exom_from_str(exon_str, start, end, strand):
        """ Get Intron Area From Exon String
        Example:
            exon_str = "489361-489710,485040-485208,476887-476945,476738-476882"
            start = 476738
            end = 489710
            Intron_From_Exom(exon_str, start, end, strand="-")
        """
        import copy
        genomeCoorIntron = []
        exonList = [ [int(exon_region.split('-')[0]), int(exon_region.split('-')[1])] for exon_region in exon_str.split(',') ]
        if len(exonList) == 0: raise NoExonError
        exonList.sort(key=lambda x: x[0], reverse=False)
        for idx in range(len(exonList)):
            if strand == '+':
                gene_start, gene_end = exonList[idx][0]-start+1, exonList[idx][1]-start+1
            elif strand == '-':
                gene_end, gene_start = end-exonList[idx][0]+1, end-exonList[idx][1]+1
            genomeCoorIntron.append([gene_start, gene_end])
        genomeCoorIntron.sort(key=lambda x:x[0])
        return exonList, genomeCoorIntron
    gene_start = geneParser[geneID]['start']
    gene_end = geneParser[geneID]['end']
    strand = geneParser[geneID]['strand']
    IsoformExons = {}
    for transID in geneParser[geneID]['transcript']:
        exon_str = GAPer[transID]['exon_str']
        try:
            IsoformExons[transID] = Exom_from_str(exon_str, gene_start, gene_end, strand)
        except NoExonError:
            sys.stderr.writelines("Warning: %s have no Exon_str, Skip it\n" % (transID, ))
    return IsoformExons

def eval_start_stop_codon(GAPer, Sequence):
    start_minus_3 = 0; start_0 = 0; start_plus_3 = 0
    stop_minus_3 = 0; stop_0 = 0; stop_plus_3 = 0
    
    stop_list = ('TAA', 'TAG', 'TGA')
    
    for tid in GAPer:
        try:
            sequence = Sequence[tid]
        except KeyError:
            continue
        
        transFeature = GAPer[tid]
        if transFeature['gene_type'] not in ('mRNA', 'protein_coding'):
            continue
        
        cds_start = transFeature['cds_start']
        cds_end = transFeature['cds_end']
        trans_len = transFeature['trans_len']
        
        if cds_start > 10 and cds_end+10<trans_len:
            start_1 = sequence[cds_start-4:cds_start-1]
            start_2 = sequence[cds_start-1:cds_start+2]
            start_3 = sequence[cds_start+2:cds_start+5]
            stop_1 = sequence[cds_end-6:cds_end-3]
            stop_2 = sequence[cds_end-3:cds_end]
            stop_3 = sequence[cds_end:cds_end+3]
            if start_1 == 'ATG': start_minus_3 += 1
            if start_2 == 'ATG': start_0 += 1
            if start_3 == 'ATG': start_plus_3 += 1
            if stop_1 in stop_list: stop_minus_3 += 1
            if stop_2 in stop_list: stop_0 += 1
            if stop_3 in stop_list: stop_plus_3 += 1
    
    print("===========Check Start/Stop Codon Position===========")
    print( "%s\t%s\t%s" % (start_minus_3, start_0, start_plus_3) )
    print( "%s\t%s\t%s" % (stop_minus_3, stop_0, stop_plus_3) )
    if start_0 < start_minus_3 or start_0 < start_plus_3:
        sys.stderr.writelines("Unexpected Error: start codon must be the first 3 bases in CDS\n")
    if stop_minus_3 > stop_0 and stop_minus_3 > stop_plus_3:
        return -3
    elif stop_0 > stop_minus_3 and stop_0 > stop_plus_3:
        return 0
    elif stop_plus_3 > stop_minus_3 and stop_plus_3 > stop_0:
        return 3

def showRNAStructure(transFeature, sequence, bias=0):
    
    cds_start = transFeature['cds_start']
    cds_end = transFeature['cds_end']
    UTR5 = sequence[:cds_start-1]
    Start_codon = sequence[cds_start-1:cds_start+2]
    CDS = sequence[ cds_start+2:cds_end-3+bias ]
    Stop_codon = sequence[ cds_end-3+bias:cds_end+bias ]
    UTR3 = sequence[cds_end+bias:]
    
    if len(CDS)%3 != 0:
        sys.stderr.writelines("Warning: The length of the CDS is not an integer multiple of 3\n")
    
    formatCDS = ""
    i = 0; j = 0
    while i<len(CDS):
        if j % 2 == 0:
            formatCDS += Colors.f(CDS[i:i+3], fc='black', bc='yellow', ft='normal')
        else:
            formatCDS += Colors.f(CDS[i:i+3], fc='black', bc='lightyellow', ft='normal')
        i += 3
        j += 1
    
    UTR5 = Colors.f(UTR5, fc='green', bc='default', ft='normal')
    Start_codon = Colors.f(Start_codon, fc='lightred', bc='default', ft='bold')
    Stop_codon = Colors.f(Stop_codon, fc='lightred', bc='default', ft='bold') 
    UTR3 = Colors.f(UTR3, fc='green', bc='default', ft='normal')
    
    return UTR5+Start_codon+formatCDS+Stop_codon+UTR3

def labelRNAPosition(Parser, tid, region, bn=None, bw=None):
    """
    bn: Bin number (default: 50)
    bw: Bin width
    bn and bw can are mutate
    """
    import math
    
    ft = Parser.getTransFeature(tid)
    
    tLen = ft['trans_len']
    cds_start = ft['cds_start']
    cds_end = ft['cds_end']
    
    s,e = region
    assert 1 <= s <= e <= tLen
    
    gene_type = "other"
    if cds_start > 1 or cds_end < tLen:
        gene_type = "mRNA"
    if ft['gene_type'] in ('mRNA', 'protein_coding'):
        gene_type = "mRNA"
    
    if not bw:
        if not bn:
            bn = 50
        len_item = 1.0*tLen/bn
    else:
        len_item = bw
        bn = int(math.ceil(1.0*tLen/len_item))
    
    start_block = min(round(s/len_item),bn-1)
    end_block = min(round(e/len_item),bn-1)
    
    cds_start_block = min(round(cds_start/len_item),bn-1)
    cds_end_block = min(round(cds_end/len_item),bn-1)
    
    string = ""
    if gene_type == 'mRNA':
        
        #### Correct CDS boundary
        if s > cds_end and start_block == cds_end_block:
            start_block = min(start_block+1, bn-1)
            if end_block < start_block: end_block = start_block
        
        if e > cds_end and end_block == cds_end_block:
            end_block = min(cds_end_block+1,bn-1)
        
        if e < cds_start and end_block == cds_start_block:
            end_block = max(end_block-1,0)
            if start_block > end_block: start_block = end_block
        
        if s < cds_start and start_block == cds_start_block:
            start_block = max(cds_start_block-1,0)
        
        #### Label
        for i in range(bn):
            if i<cds_start_block or i>cds_end_block:
                c = "-"
            else:
                c = "|"
            if start_block<=i<=end_block:
                c = Colors.f(c, fc='red', ft='bold')
                #c = "\x1b[1;31;40m" + c + "\x1b[0m"
            string += c
    else:
        for i in range(bn):
            c = "-"
            if start_block<=i<=end_block:
                c = Colors.f(c, fc='red', ft='bold')
                #c = "\x1b[1;31;40m" + c + "\x1b[0m"
            string += c
    
    return string


def getRNAPosition(Parser, tid, region):
    """
    return type: [ "not_mRNA", "5UTR", "span_5UTR_CDS", "CDS", "span_CDS_3UTR", "3UTR", "span_5UTR_CDS_3UTR", "INVALID" ]
    """
    
    ft = Parser.getTransFeature(tid)
    
    tLen = ft['trans_len']
    cds_start = ft['cds_start']
    cds_end = ft['cds_end']
    
    s,e = region
    assert 1 <= s <= e <= tLen
    
    gene_type = "other"
    if cds_start > 1 or cds_end < tLen:
        gene_type = "mRNA"
    if ft['gene_type'] in ('mRNA', 'protein_coding'):
        gene_type = "mRNA"
    
    if gene_type != 'mRNA': return "not_mRNA"
    if e<cds_start: return "5UTR"
    if s<cds_start<=e<=cds_end: return "span_5UTR_CDS"
    if cds_start<=s<=e<=cds_end: return "CDS"
    if cds_start<=s<=cds_end<e: return "span_CDS_3UTR"
    if cds_end<s: return "3UTR"
    if s<cds_start<cds_end<e: return "span_5UTR_CDS_3UTR"
    return "INVALID"

def getLenSortedTransDictForGenes(Parser, only=[]):
    """
    return a dict like:
        gene_id: [ trans_id1, trans_id1, trans_id3... ]
        trans_idn is sorted by their length
    only constraint RNA type such as
        only = ['mRNA', 'protein_coding', 'snRNA'...]
    """
    trans_dict = {}
    
    gene_parser = Parser.getGeneParser()
    for gid in gene_parser:
        trans_dict[gid] = []
        for tid in gene_parser[gid]['transcript']:
            ft = Parser.GAPer[tid]
            length = ft['trans_len']
            if only:
                gene_type = ft['gene_type']
                if gene_type in only:
                    trans_dict[gid].append( (tid,length) )
            else:
                trans_dict[gid].append( (tid,length) )
        if trans_dict[gid]:
            trans_dict[gid].sort(key=lambda x: x[1], reverse=True)
            trans_dict[gid] = [ it[0] for it in trans_dict[gid] ]
        else:
            del trans_dict[gid]
    
    return trans_dict

def getGeneCombinedIntronExon(Parser, geneID, verbose=True):
    """
    Parse gene intron/exon regions:
        1. Combine exons from all transcripts, we defined them as exon regions;
        2. Remove the exon regions from gene regions, we defined them as intron regions
    
    return: intron_regions, exon_regions
    """
    
    geneExons = Parser.getGeneExon(geneID)
    
    # 1. collect
    exons = []
    for tid in geneExons:
        exons += geneExons[tid][1]
    exons.sort(key=lambda x: x[0])
    if len(exons) == 0: return None, None
    
    # 2. combine
    new_exons = [ exons[0] ]
    for i in range(1, len(exons)):
        if new_exons[-1][1] >= exons[i][0]:
            new_exons[-1][1] = max(new_exons[-1][1], exons[i][1])
        else:
            new_exons.append(exons[i])
    
    # 3. get intron
    new_introns = []
    if new_exons[0][0] != 1:
        if verbose: print(geneExons)
        raise Exception("Exons is not started from 1")
    for i in range(1, len(new_exons)):
        new_introns.append( (new_exons[i-1][1]+1, new_exons[i][0]-1) )
    if new_exons[-1][1] != Parser.getGeneParser()[geneID]['length']:
        if verbose: print(geneExons)
        raise Exception("Exons is not ended in geneLength")
    
    # 4. return
    return new_introns, new_exons


class ParseTransClass(object):
    def __init__(self, genomeCoorBedFile, seqFileName='', showAttr=True, remove_tid_version=False, remove_gid_version=False):
        """
        genomeCoorBedFile   -- A *.genomeCoor.bed file produced by parseGTF.py
        seqFileName         -- Transcriptome fasta file produced by parseGTF.py
        showAttr            -- Show a example
        remove_tid_version  -- Remove version information. ENST000000022311.2 => ENST000000022311
        remove_gid_version  -- Remove version information. ENSG000000022311.2 => ENSG000000022311
        """
        self.seqFileName = seqFileName
        self.genomeCoorBedFile = genomeCoorBedFile
        self.GAPer = ParseGenomeCoorBedFile(self.genomeCoorBedFile, showAttr=showAttr, remove_tid_version=remove_tid_version, remove_gid_version=remove_gid_version)
        if self.seqFileName != '':
            self.Seq = readSeq(self.seqFileName, remove_tid_version=remove_tid_version)
            self.stop_bias = eval_start_stop_codon(self.GAPer, self.Seq)
        
        sys.stderr.writelines('Warning: Trans/Gene Coorninate System: [1_based_start, 1_based_end]\n')
    
    def addSeq(self, seqFileName, remove_tid_version=False):
        """
        seqFileName         -- Transcriptome fasta file produced by parseGTF.py
        remove_tid_version  -- Remove version information. ENST000000022311.2 => ENST000000022311
        """
        self.seqFileName = seqFileName
        self.Seq = readSeq(self.seqFileName, remove_tid_version=remove_tid_version)
        self.stop_bias = eval_start_stop_codon(self.GAPer, self.Seq)
    
    def getTransFeature(self, transID, showAttr=False, verbose=True):
        """
        transID             -- Transcript ID
        showAttr            -- Show an example
        verbose             -- Print the information when transcript not found
        
        Return a dictionary:
            chr             -- Genome chromosome
            strand          -- + or -
            start           -- Genome start
            end             -- Genome end
            
            gene_name       -- Gene symbol
            gene_id         -- Gene id
            gene_type       -- Gene type
            
            trans_len       -- Transcript length
            
            utr_5_start     -- Transcript-based start site of 5'UTR
            utr_5_end       -- Transcript-based end site of 5'UTR
            cds_start       -- Transcript-based start site of CDS
            cds_end         -- Transcript-based end site of CDS
            utr_3_start     -- Transcript-based start site of 3'UTR
            utr_3_end       -- Transcript-based end site of 3'UTR
            
            exon_str        -- Genome-based exon string
        """
        
        try:
            transFeature = self.GAPer[transID]
        except KeyError:
            if verbose: sys.stderr.writelines('Warning: no this transID: %s in GAPer\n' % (transID, ))
            raise KeyError
        
        if showAttr:
            showBeautifulExample(transFeature, title='%s information' % (transID, ))
        
        return transFeature
    
    def getGeneParser(self, showAttr=True):
        """
        showAttr            -- Show an example
        
        Return a dictionary of dictionaries:
            { geneID => { ... }, ... } includes
                chr         -- Genome chromosome
                strand      -- + or -
                start       -- Genome start
                end         -- Genome end
                
                gene_name   -- Gene symbol
                gene_type   -- All transcript types
                
                length      -- Gene length (end-start+1)
                
                transcript  -- All transcripts belong to this gene
        """
        
        if not hasattr(self, 'GeneParser'):
            self.GeneParser = biuldGeneParser(self.GAPer, showAttr=showAttr)
        
        return self.GeneParser
    
    def getGeneIntron(self, geneID):
        """
        geneID              -- Gene id
        
        Get introns of all transcripts from this gene:
            { transID => [[intron1_start, intron1_end], [intron2_start, intron2_end], [intron3_start, intron3_end]],... }
        """
        
        return GeneIntron(self.GAPer, self.getGeneParser(), geneID)
    
    def getGeneExon(self, geneID):
        """
        geneID              -- Gene id
        
        Get introns of all transcripts from this gene:
            { transID => [[exon1_start, exon1_end], [exon2_start, exon2_end], [exon3_start, exon3_end]...],... }
        """
        
        return GeneExon(self.GAPer, self.getGeneParser(), geneID)
    
    def getGeneBin(self):
        """
        Get bins for gene coordination conversion
        """
        if not hasattr(self, 'GeneBin'):
            self.GeneBin = CoorFunc.binize_Gene(self, bw=100000)
        
        return self.GeneBin
    
    def getTransBin(self):
        """
        Get bins for transcript coordination conversion
        """
        if not hasattr(self, 'TransBin'):
            self.TransBin = CoorFunc.binize_Trans(self, bw=100000)
        
        return self.TransBin
    
    def geneCoor2genomeCoor(self, geneID, pos):
        """
        geneID              -- Gene id
        pos                 -- Gene-based position
        
        Convert gene-based coordinates to genome-based coordinates
        Return: [chrID, chrPos, Strand]
        """
        return CoorFunc.geneCoor2genomeCoor(self, geneID, pos)
    
    def genomeCoor2geneCoor(self, chrID, start, end, strand):
        """
        chrID               -- Chromosome id
        start               -- Chromosome-based start position
        end                 -- Chromosome-based end position
        strand              -- Chromosome strand
        
        Convert genome-based coordinates to gene-based coordinates
        Return [ [chrID, chrStart, chrEnd, geneID, geneStart, geneEnd], ... ]
        """
        GeneBin = self.getGeneBin()
        return CoorFunc.genomeCoor2geneCoor(GeneBin, self, chrID, start, end, strand)
    
    def genomeCoor2transCoor(self, chrID, start, end, strand):
        """
        chrID               -- Chromosome id
        start               -- Chromosome-based start position
        end                 -- Chromosome-based end position
        strand              -- Chromosome strand
        
        Convert genome-based coordinates to transcript-based coordinates
        Return [ [chrID, chrStart, chrEnd, transID, transStart, transEnd], ... ]
        """
        TransBin = self.getTransBin()
        return CoorFunc.genomeCoor2transCoor(TransBin, self, chrID, start, end, strand)
    
    def transCoor2genomeCoor(self, transID, pos):
        """
        transID             -- Transcript id
        pos                 -- Transcript-based position
        
        Convert transcript-based coordinates to genome-based coordinates
        Return [chrID, chrPos, Strand]
        """
        return CoorFunc.transCoor2genomeCoor(self, transID, pos)
    
    def geneCoor2transCoor(self, geneID, start, end):
        """
        geneID              -- Gene id
        start               -- Gene-based start position
        end                 -- Gene-based end position
        
        Convert gene-based coordinates to transcript-based coordinates
        Return  [chrID, chrStart, chrEnd, transID, transStart, transEnd], ... ]
        """
        TransBin = self.getTransBin()
        return CoorFunc.geneCoor2transCoor(TransBin, self, geneID, start, end)
    
    def transCoor2geneCoor(self, transID, start, end):
        """
        transID             -- Transcript id
        start               -- Transcript-based start position
        end                 -- Transcript-based end position
        
        Convert transcript-based coordinates to gene-based coordinates
        return [ [chrID, chrStart, chrEnd, geneID, geneStart, geneEnd], ... ]
        """
        GeneBin = self.getGeneBin()
        return CoorFunc.transCoor2geneCoor(GeneBin, self, transID, Start, End)
    
    def getTransList(self):
        """
        Return list of all transcript id
        """
        return self.GAPer.keys()
    
    def getGeneList(self):
        """
        Return list of all gene id
        """
        return self.getGeneParser(showAttr=False).keys()
    
    def getmRNATransList(self):
        """
        Return list of all mRNA id
        """
        mRNATransList = []
        
        for tid in self.GAPer:
            if self.GAPer[tid]['gene_type'] in ('mRNA', 'protein_coding'):
                mRNATransList.append(tid)
        
        return mRNATransList
    
    def getmRNAGeneList(self):
        """
        Return list of all mRNA gene id
        """
        mRNATransList = self.getmRNATransList()
        mRNAGeneSet = set([ self.GAPer[tid]['gene_id'] for tid in mRNATransList ])
        return list(mRNAGeneSet)
    
    def getLenSortedTransDictForGenes(self, only=[]):
        """
        only                --  Contraint transcript types, such as mRNA, snRNA ...

        Return a dict:
            { geneID => [ transID1, transID2, ... ], ...}
            transcripts are sorted by length
        """
        return getLenSortedTransDictForGenes(self, only=[])
    
    def getTransByGeneID(self, geneID):
        """
        geneID              -- Gene id
        
        Return transcripts belong to specific gene
        """
        gene_parser = self.getGeneParser(showAttr=False)
        return gene_parser[geneID]['transcript']
    
    def getGeneByGeneName(self, geneName):
        """
        geneName            -- Gene symbol
        
        Return gene id belong with specific gene name
        """
        if not hasattr(self, 'GeneNameToGeneID'):
            GeneParser = self.getGeneParser(showAttr=False)
            self.GeneNameToGeneID = { GeneParser[gid]['gene_name']:gid for gid in GeneParser }
        
        return self.GeneNameToGeneID[geneName]
    
    def getTransByGeneName(self, geneName):
        """
        geneName            -- Gene symbol
        
        Return transcripts belong to specific gene
        """
        geneID = self.getGeneByGeneName(geneName)
        return self.getGeneParser()[geneID]['transcript']
    
    def getTransByGeneType(self, geneType):
        """
        geneType            -- Gene type
        
        Return a list of transcripts belong to specific gene type
        """
        trans_list = []
        
        for tid in self.GAPer:
            if self.GAPer[tid]['gene_type'] == geneType:
                trans_list.append(tid)
        
        return trans_list
    
    def getTransSeq(self, transID):
        """
        transID             -- Transcript id
        
        Return transcript sequence
        """
        if self.seqFileName == "":
            sys.stderr.writelines('Not define Seq! Use addSeq() to add Seq\n')
            raise KeyError
        
        return self.Seq[transID]
    
    def getGeneCombinedIntronExon(self, geneID, verbose=True):
        """
        geneID          -- Gene id
        verbose         -- Print error information when error occured
        
        Parse gene intron/exon regions:
            1. Combine exons from all transcripts, they are defined as exon regions
            2. Remove the exon regions from gene regions, they are defined as intron regions
        
        return: [intron_regions, exon_regions]
        """
        
        return getGeneCombinedIntronExon(self, geneID=geneID, verbose=verbose)
    
    def showRNAStructure(self, transID):
        """
        transID         -- Transcript id
        
        Return string of mRNA sequence with color labeled its UTR/CDS/Codons
        """
        if self.seqFileName != "":
            ft = self.getTransFeature(transID)
            seq = self.getTransSeq(transID)
            return showRNAStructure(ft, seq, bias=self.stop_bias)
        else:
            sys.stderr.writelines('Not define Seq! Use addSeq() to add Seq\n')
            raise KeyError
    
    def labelRNAPosition(self, transID, region, bn=None, bw=None):
        """
        transID         -- Transcript id
        region          -- [start, end]
        bn              -- Bin number (default: 50)
        bw              -- Bin width
        
        return RNA structure and label the region
        -------|||||||||||||||||--------------------------
         5'UTR        CDS                3'UTR
        """
        return labelRNAPosition(self, transID, region, bn=bn, bw=bw)
    
    def getRNAPosition(self, transID, region):
        """
        transID         -- Transcript id
        region          -- [start, end]
        
        Return any one of:
            -> not_mRNA
            -> 5UTR
            -> span_5UTR_CDS
            -> CDS
            -> span_CDS_3UTR
            -> 3UTR
            -> span_5UTR_CDS_3UTR
            -> INVALID
        """
        return getRNAPosition(self, transID, region)

