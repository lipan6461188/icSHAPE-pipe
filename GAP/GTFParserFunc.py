#-*- coding:utf-8 -*-

import sys, os, random

###############################
## Shared Functions
###############################

def load_getstatusoutput():
    import subprocess
    if 'getoutput' in dir(subprocess):
        return subprocess.getstatusoutput
    else:
        import commands
        return commands.getstatusoutput

getstatusoutput = load_getstatusoutput()

def format_UTR(exonString, cdsString, strand, verbose=True):
    def correctList(List, strand):
        if strand == '+':    
            List.sort(key=lambda x: x[0], reverse=False)
            i = 0
            while i<len(List)-1:
                if List[i][1] == List[i+1][0]:
                    List[i][1] = List[i+1][1]
                    del List[i+1]
                else:
                    i += 1
        
        if strand == '-':
            List.sort(key=lambda x: x[0], reverse=True)
            i = 0
            while i<len(List)-1:
                if List[i][0] == List[i+1][1]:
                    List[i][0] = List[i+1][0]
                else:
                    i += 1
    
    exonList = []
    for region_str in exonString.split(','):
        item = region_str.split('-')
        exonList.append( [int(item[0]), int(item[1])] )
    
    cdsList = []
    for region_str in cdsString.split(','):
        item = region_str.split('-')
        cdsList.append( [int(item[0]), int(item[1])] )
    
    correctList(exonList, strand)
    correctList(cdsList, strand)
    
    utr_5 = []; utr_3 = []
    
    if strand == '+':
        for (exon_start, exon_end) in exonList:
            if exon_end < cdsList[0][0]:
                utr_5.append( (exon_start, exon_end) )
            elif (exon_start < cdsList[0][0] and exon_end == cdsList[0][1]) or (exon_start < cdsList[0][0] and len(cdsList) == 1):
                utr_5.append( (exon_start, cdsList[0][0]-1) )
                break
            elif exon_start == cdsList[0][0]:
                break
            else:
                if verbose:
                    sys.stderr.writelines('Impossible Event 1 '+"\n")
                    sys.stderr.writelines('exonString: '+exonString+"\n")
                    sys.stderr.writelines('cdsString: '+cdsString+"\n")
                    sys.stderr.writelines(strand+"\n")
                raise Exception("Impossible Event 1")
                break
        for (exon_start, exon_end) in exonList[::-1]:
            if exon_start > cdsList[-1][1]:
                utr_3.append( (exon_start, exon_end) )
            elif (cdsList[-1][1] < exon_end and exon_start == cdsList[-1][0]) or (cdsList[-1][1] < exon_end and len(cdsList) == 1):
                utr_3.append( (cdsList[-1][1]+1, exon_end) )
                break
            elif exon_end == cdsList[-1][1]:
                break
            else:
                if verbose:
                    sys.stderr.writelines('Impossible Event 2'+"\n")
                    sys.stderr.writelines('exonString: '+exonString+"\n")
                    sys.stderr.writelines('cdsString: '+cdsString+"\n")
                    sys.stderr.writelines(strand+"\n")
                raise Exception("Impossible Event 2")
                break
        utr_3.reverse()
    else:
        for (exon_start, exon_end) in exonList:
            if exon_start > cdsList[0][1]:
                utr_5.append( (exon_start, exon_end) )
            elif (cdsList[0][1] < exon_end and exon_start == cdsList[0][0]) or (cdsList[0][1] < exon_end and len(cdsList) == 1):
                utr_5.append( (cdsList[0][1]+1, exon_end) )
                break
            elif exon_end == cdsList[0][1]:
                break
            else:
                if verbose:
                    sys.stderr.writelines('Impossible Event 3'+"\n")
                    sys.stderr.writelines('exonString: '+exonString+"\n")
                    sys.stderr.writelines('cdsString: '+cdsString+"\n")
                    sys.stderr.writelines(strand+"\n")
                raise Exception("Impossible Event 3")
                break
        for (exon_start, exon_end) in exonList[::-1]:
            if exon_end < cdsList[-1][0]:
                utr_3.append( (exon_start, exon_end) )
            elif (exon_start < cdsList[-1][0] and exon_end == cdsList[-1][1]) or (exon_start < cdsList[-1][0] and len(cdsList) == 1):
                utr_3.append( (exon_start, cdsList[-1][0]-1) )
                break
            elif exon_start == cdsList[-1][0]:
                break
            else:
                if verbose:
                    sys.stderr.writelines('Impossible Event 4'+"\n")
                    sys.stderr.writelines('exonString: '+exonString+"\n")
                    sys.stderr.writelines('cdsString: '+cdsString+"\n")
                    sys.stderr.writelines(strand+"\n")
                raise Exception("Impossible Event 4")
                break
        utr_3.reverse()
    
    return utr_5, utr_3

def writeTranscriptome(genomeCoorFileName, genomeFaFile, outTranscriptomeFaFile, verbose=True, showAttr=False):
    import GAP
    import SeqFetch
    
    Parser = GAP.init(genomeCoorFileName, showAttr=showAttr)
    geneParser = Parser.getGeneParser(showAttr=showAttr)
    
    SeqFetcher = SeqFetch.SeqClass(genomeFaFile)
    
    OUT = open(outTranscriptomeFaFile, 'w')
    
    no_chr_set = set()
    for geneID in geneParser:
        for tid in geneParser[geneID]['transcript']:
            ft = Parser.getTransFeature(tid)
            geneType = ft['gene_type']
            geneName = ft['gene_name']
            seqAnno = "%s|%s|%s" % (geneID, geneName, geneType)
            
            strand = ft['strand']
            ChrID = ft['chr']
            
            if not SeqFetcher.has(ChrID):
                if ChrID not in no_chr_set:
                    no_chr_set.add(ChrID)
                    if verbose:
                        sys.stderr.writelines('Warning: %s not in %s' % (ChrID, genomeFaFile)+"\n")
                continue
            
            exonsList = Parser.getGeneExon(geneID)[tid][0]
            if strand == '-': exonsList.sort(key=lambda x: x[0], reverse=True)
            
            RNA_seq = ""
            for exon in exonsList:
                try:
                    RNA_seq += SeqFetcher.fetch(ChrID, exon[0]-1, exon[1], strand)
                except KeyError:
                    if verbose:
                        sys.stderr.writelines('Warning: KeyError -> %s\t%d\t%d\t%s' % (ChrID, exon[0]-1, exon[1], strand)+"\n")
                    continue
            
            OUT.writelines('>%s\t%s\n%s\n' % (tid, seqAnno, SeqFetch.cutSeq(RNA_seq)))
    
    OUT.close()

def is_scaffold(chrID):
    if chrID.startswith('chr') or chrID.startswith('NC_') or len(chrID)<=6:
        return False
    else:
        return True

###############################
## Ensembl GTF Functions
###############################

def parse_gtf_attributes(attributesString):
    attributes = {}
    for Tuple in attributesString.strip(';').split('; '):
        try:
            Key_Value = Tuple.strip(' ').split(' ')
            Key = Key_Value[0]
            Value = ' '.join(Key_Value[1:])
            attributes[ Key ] = Value.strip('"')
        except ValueError:
            print('Warning: '+Tuple+' cannot be unpacked, skip it')
    return attributes

def read_ensembl_gtf(gtfFileName, rem_scaffold=False):
    """
    gtfFileName             -- Ensembl/Gencode GTF file
    rem_scaffold            -- Remove scaffolds. scaffolds are defined as 
                               those chromosomes with id length > 6 and not startswith chr and NC_
    """
    gtf_container = {'gene':{}, 'RNA':{}, 'CDS':{}, 'exon':{}}
    
    lineCount = 0
    for line in open(gtfFileName):
        lineCount += 1
        if lineCount % 100000 == 0: print('\tRead %d Lines...' % (lineCount,))
        if line[0] == '#': continue
        
        (Chr, Source, Type, Chr_Start, Chr_End, Score, Strand, Phase, Attributes) = line.strip().split('\t')
        if rem_scaffold and is_scaffold(Chr):
            continue
        
        attributes = parse_gtf_attributes(Attributes)
        attributes['chr'] = Chr; attributes['start'] = Chr_Start; attributes['end'] = Chr_End; attributes['strand'] = Strand
        
        if Type == 'gene':
            if 'gene_name' not in attributes:
                attributes['gene_name'] = 'NULL'
            gtf_container[ 'gene' ][ attributes['gene_id'] ] = attributes
        elif Type == 'transcript':
            if 'gene_name' not in attributes:
                attributes['gene_name'] = 'NULL'
            if 'transcript_biotype' not in attributes and 'transcript_type' not in attributes:
                attributes['trans_type'] = 'NULL'
            else:
                if 'transcript_biotype' in attributes:
                    attributes['trans_type'] = attributes['transcript_biotype']
                    del attributes['transcript_biotype']
                else:
                    attributes['trans_type'] = attributes['transcript_type']
                    del attributes['transcript_type']
            gtf_container[ 'RNA' ][ attributes['transcript_id'] ] = attributes
        elif Type == 'exon':
            gtf_container[ 'exon' ][ attributes['transcript_id'] ] = gtf_container[ 'exon' ].get( attributes['transcript_id'], []) + [attributes]
        elif Type == 'CDS':
            gtf_container[ 'CDS' ][ attributes['transcript_id'] ] = gtf_container[ 'CDS' ].get( attributes['transcript_id'], []) + [attributes]
    
    return gtf_container

def write_gtf_genomeCoor_bed(gtf_container, genomeCoorFileName):
    
    tmpFileName = genomeCoorFileName + ".tmp"
    
    TMP = open(tmpFileName, 'w')
    for rna_ID in gtf_container['RNA']:
        transItem = get_format_Ensembl_gtf_RNA_item(gtf_container, rna_ID)
        if transItem:
            string = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (transItem[0], transItem[1], transItem[2], transItem[3], transItem[4], transItem[5], transItem[6], transItem[7])
            if len(transItem) == 9:
                string += "\t" + transItem[8]
            TMP.writelines(string+"\n")
    
    TMP.close()
    
    return_code, output = getstatusoutput( "sort -k 1,1 -k2n,3n %s > %s" % (tmpFileName, genomeCoorFileName) )
    if return_code != 0:
        print('Error: cat and sort step is broken. Please check the file: '+tmpFileName)
    else:
        os.remove(tmpFileName)
        print('Success!')

def get_format_Ensembl_gtf_RNA_item(gtf_container, rna_ID, verbose=True):
    
    RNA = gtf_container['RNA'][rna_ID]
    
    chrID = RNA['chr']
    GeneType = RNA['trans_type']
    GeneID = RNA['gene_id']
    TransID = rna_ID
    
    (Start, End, Strand, GeneID, TransID, GeneType) = (RNA['start'], RNA['end'], RNA['strand'], RNA['gene_name']+'='+GeneID, TransID, GeneType)
    
    exons = gtf_container['exon'][rna_ID]
    exonsList = [ exon['start']+'-'+exon['end'] for exon in exons ]
    exonString = ','.join(exonsList)
    
    format_list = [ chrID, Start, End, Strand, GeneID, TransID, GeneType, exonString ]
    
    # UTR
    if rna_ID in gtf_container['CDS']:
        CDSs = gtf_container['CDS'][rna_ID]
        cdsList = [ cds['start']+'-'+cds['end'] for cds in CDSs ]
        cdsString = ','.join(cdsList)
        
        try:
            (utr_5_list, utr_3_list) = format_UTR(exonString, cdsString, Strand, verbose=verbose)
        
        except:
            if verbose:
                sys.stderr.writelines('Skip this transcript: '+rna_ID+"\n")
                sys.stderr.writelines('exonString: '+exonString+"\n")
                sys.stderr.writelines('cdsString: '+cdsString+"\n")
                sys.stderr.writelines('strand: '+Strand+"\n")
                sys.stderr.writelines('\n\n')
            return None
        
        utr_5_string = ','.join(["%s-%s"%(it[0], it[1]) for it in utr_5_list])
        utr_3_string = ','.join(["%s-%s"%(it[0], it[1]) for it in utr_3_list])
        utrString = (utr_5_string+','+utr_3_string).strip(',')
        return format_list + [utrString]
    
    else:
        return format_list


###############################
## NCBI GFF3 Functions
###############################

def parse_gff3_attributes(attributesString):
    attributes = {}
    for Tuple in attributesString.split(';'):
        (Key, Value) = Tuple.split('=')
        if Key == 'Dbxref':
            for miniTuple in Value.split(','):
                miniKey_miniValue = miniTuple.split(':')
                if miniKey_miniValue[0] in attributes:
                     attributes[ miniKey_miniValue[0]+'_2' ] = miniKey_miniValue[-1]
                else:
                    attributes[ miniKey_miniValue[0] ] = miniKey_miniValue[-1]
        else:
            attributes[ Key ] = Value
    return attributes

def read_ncbi_gff3(gff3FileName, rem_scaffold=False, raw_chrID=False):
    """
    gtfFileName             -- NCBI GFF3 file
    rem_scaffold            -- Remove scaffolds. scaffolds are defined as 
                               those chromosomes with id length > 6 and not startswith chr and NC_
    raw_chrID               -- Don't convert NC_* to chrXX code
    """
    gff3_container = {'region': {}, 'others':{}, 'RNA':{}, 'exon':{}, 'CDS':{}}
    RNA_ids = set()
    mRNA_ids = set()
    
    lineCount = 0
    for line in open(gff3FileName):
        lineCount += 1
        if lineCount % 100000 == 0: print('\tRead %d Lines...' % (lineCount,))
        if line[0] == '#': continue
        
        (Chr, Source, Type, Chr_Start, Chr_End, Score, Strand, Phase, Attributes) = line.strip().split('\t')
        
        if not raw_chrID:
            if rem_scaffold and is_scaffold(Chr):
                continue
        
        attributes = parse_gff3_attributes(Attributes)
        attributes['chr'] = Chr; attributes['start'] = Chr_Start; attributes['end'] = Chr_End; attributes['strand'] = Strand
        
        if Type == 'region':
            if attributes['start'] == '1':
                gff3_container[ 'region' ][ Chr ] = attributes
        
        elif Type == 'exon':
            RNA_id = attributes['Parent']
            gff3_container[ 'exon' ][ RNA_id ] = gff3_container[ 'exon' ].get( RNA_id, []) + [attributes]
            RNA_ids.add(RNA_id)
        
        elif Type == 'CDS':
            RNA_id = attributes['Parent']
            gff3_container[ 'CDS' ][ RNA_id ] = gff3_container[ 'CDS' ].get( RNA_id, []) + [attributes]
            mRNA_ids.add(RNA_id)
        
        else:
            attributes['Type'] = Type
            item_id = attributes['ID']
            gff3_container['others'][ item_id ] = attributes
    
    known_gene_types = ("lnc_RNA", "misc_RNA", "miRNA", "pseudogene", "rRNA", "snoRNA", "snRNA", "tRNA", "precursor_RNA", "antisense_RNA", "RNase_P_RNA", "telomerase_RNA", "RNase_MRP_RNA", "protein_coding")
    
    for RNA_id in RNA_ids:
        RNA = gff3_container['others'][ RNA_id ]
        if 'gene' not in RNA:
            if 'Name' in RNA:
                RNA['gene'] = RNA['Name']
            else:
                RNA['gene'] = 'NULL'
        
        if 'GeneID' not in RNA: RNA['GeneID'] = 'NULL'
        
        if RNA_id in mRNA_ids:
            RNA[ 'gene_type' ] = 'mRNA'
        elif RNA[ 'Type' ] in known_gene_types:
            RNA[ 'gene_type' ] = RNA[ 'Type' ]
        elif "gene_biotype" in  RNA and RNA[ 'gene_biotype' ] in known_gene_types:
            RNA[ 'gene_type' ] = RNA[ 'gene_biotype' ]
        elif "ncrna_class" in RNA and RNA[ 'ncrna_class' ] in known_gene_types:
            RNA[ 'gene_type' ] = RNA[ 'ncrna_class' ]
        elif "gbkey" in RNA and RNA[ 'gbkey' ] in known_gene_types:
            RNA[ 'gene_type' ] = RNA[ 'gbkey' ]
        else:
            if 'Parent' in RNA:
                parent = gff3_container[ 'others' ][ RNA['Parent'] ]
                if 'gene_biotype' in parent:
                    RNA[ 'gene_type' ] = parent[ 'gene_biotype' ]
                else:
                    RNA[ 'gene_type' ] = RNA['Type']
            else:
                RNA[ 'gene_type' ] = RNA['Type']
        
        gff3_container['RNA'][RNA_id] = RNA
        
    return gff3_container

def write_gff3_genomeCoor_bed(gff3_container, genomeCoorFileName, raw_chrID=False):
    
    tmpFileName = genomeCoorFileName + ".tmp"
    
    NCToChr = build_NC_To_chr_dict(gff3_container, raw_chrID=raw_chrID)
    rnaIDToTransID = build_rnaID_To_transID(gff3_container)
    
    TMP = open(tmpFileName, 'w')
    for rna_ID in gff3_container['RNA']:
        transItem = get_format_NCBI_gff3_RNA_item(gff3_container, rna_ID, NCToChr, rnaIDToTransID)
        if transItem:
            string = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (transItem[0], transItem[1], transItem[2], transItem[3], transItem[4], transItem[5], transItem[6], transItem[7])
            if len(transItem) == 9:
                string += "\t" + transItem[8]
            TMP.writelines(string+"\n")
    
    TMP.close()
    
    return_code, output = getstatusoutput( "sort -k 1,1 -k2n,3n %s > %s" % (tmpFileName, genomeCoorFileName) )
    if return_code != 0:
        print('Error: cat and sort step is broken. Please check the file: '+tmpFileName)
    else:
        os.remove(tmpFileName)
        print('Success!')

def get_format_NCBI_gff3_RNA_item(gff3_container, rna_ID, NCToChr, rnaIDTotransID, verbose=True):
    """
        return [Start, End, Strand, GeneID, TransID, GeneType]
    """
    
    RNA = gff3_container['RNA'][rna_ID]
    
    NCID = RNA['chr']
    if NCID.startswith('NC'): 
        chrID = NCToChr[NCID]
    else:
        chrID = NCID
    
    GeneType = RNA['gene_type']
    TransID = rnaIDTotransID[rna_ID]
    
    (Start, End, Strand, GeneID, TransID, GeneType) = (RNA['start'], RNA['end'], RNA['strand'], RNA['gene']+'='+RNA['GeneID'], TransID, GeneType)
    
    # Exon
    exons = gff3_container['exon'][rna_ID]
    exonsList = [ exon['start']+'-'+exon['end'] for exon in exons ]
    exonString = ','.join(exonsList)
    
    format_list = [ chrID, Start, End, Strand, GeneID, TransID, GeneType, exonString ]
    
    # UTR
    if rna_ID in gff3_container['CDS']:
        CDSs = gff3_container['CDS'][rna_ID]
        cdsList = [ cds['start']+'-'+cds['end'] for cds in CDSs ]
        cdsString = ','.join(cdsList)
        
        try:
            (utr_5_list, utr_3_list) = format_UTR(exonString, cdsString, Strand, verbose=verbose)
        
        except:
            if verbose:
                sys.stderr.writelines('Skip this transcript: '+rna_ID+"\n")
                sys.stderr.writelines('exonString: '+exonString+"\n")
                sys.stderr.writelines('cdsString: '+cdsString+"\n")
                sys.stderr.writelines('strand: '+Strand+"\n")
                sys.stderr.writelines('\n\n')
            return None
        
        utr_5_string = ','.join(["%s-%s"%(it[0], it[1]) for it in utr_5_list])
        utr_3_string = ','.join(["%s-%s"%(it[0], it[1]) for it in utr_3_list])
        utrString = (utr_5_string+','+utr_3_string).strip(',')
        return format_list + [utrString]
    
    else:
        return format_list

def build_rnaID_To_transID(gff3_container):
    rnaIDToTransID = {}
    
    for rna_ID in gff3_container['RNA']:
        RNA = gff3_container['RNA'][rna_ID]
        try:
            TransID = RNA['transcript_id']
        except KeyError:
            TransID = rna_ID
        rnaIDToTransID[rna_ID] = TransID
    
    return rnaIDToTransID

def build_chr_To_NC_dict(gff3_container, raw_chrID=False):
    NC_To_chr = build_NC_To_chr_dict(gff3_container, raw_chrID=raw_chrID)
    
    chrToNC = {}
    for NC,chrID in NC_To_chr.items():
        chrToNC[chrID] = NC
    
    return chrToNC

def build_NC_To_chr_dict(gff3_container, raw_chrID=False):
    NCToChr = {}
    regions = gff3_container['region']
    
    for ChrID in regions:
        if raw_chrID:
            NCToChr[ChrID] = ChrID
        elif ChrID.startswith('NC'):    
            if 'chromosome' in regions[ChrID]:
                NCToChr[ChrID] =  'chr'+regions[ChrID]['chromosome']
            elif regions[ChrID]['genome'] == 'mitochondrion':
                NCToChr[ChrID] = 'chrM'
    
    return NCToChr


###########################
## GenomeCoor To TransCoor
###########################

def norm_exons(string):
    tuples = list(); raw_tuple = list()
    arr = string.split(',')
    for indx in range(len(arr)):
        aarr = arr[indx].split('-')
        if indx == 0:  tuples.append( (1, abs(int(aarr[1])-int(aarr[0]))+1 ) )
        else: tuples.append( (tuples[-1][1]+1, abs(int(aarr[1])-int(aarr[0]))+tuples[-1][1]+1 ) )
        raw_tuple.append( (int(aarr[0]), int(aarr[1])) )
    return tuples, raw_tuple

def norm_utr(exon_str, utr_str, strand):
    norm_tuple, raw_tuple = norm_exons(exon_str)
    utr_tuple = list()
    arr = utr_str.split(',')
    for indx in range(len(arr)):
        aarr = [ int(i) for i in arr[indx].split('-') ]
        flag = False
        for i in range(len(raw_tuple)):
            if strand == '-':
                if raw_tuple[i][0] == aarr[0] <= aarr[1] == raw_tuple[i][1]:
                    utr_tuple.append( ( norm_tuple[i][0], norm_tuple[i][1] ) ); flag = True; break
                elif raw_tuple[i][0] == aarr[0] <= aarr[1] < raw_tuple[i][1]:
                    utr_tuple.append( ( raw_tuple[i][1]-aarr[1]+norm_tuple[i][0], norm_tuple[i][1] ) ); flag = True; break
                elif raw_tuple[i][0] < aarr[0] <= aarr[1] == raw_tuple[i][1]:
                    utr_tuple.append( ( norm_tuple[i][0], raw_tuple[i][1]-aarr[0]+norm_tuple[i][0] ) ); flag = True; break
            if strand == '+':
                if raw_tuple[i][0] == aarr[0] <= aarr[1] == raw_tuple[i][1]:
                    utr_tuple.append( ( norm_tuple[i][0], norm_tuple[i][1] ) ); flag = True; break
                elif raw_tuple[i][0] == aarr[0] <= aarr[1] < raw_tuple[i][1]:
                    utr_tuple.append( ( norm_tuple[i][0], norm_tuple[i][0]+aarr[1]-aarr[0] ) ); flag = True; break
                elif raw_tuple[i][0] < aarr[0] <= aarr[1] == raw_tuple[i][1]:
                    utr_tuple.append( ( norm_tuple[i][1]-(aarr[1]-aarr[0]), norm_tuple[i][1] ) ); flag = True; break
        if not flag: print('Unexpected Result')
    return norm_tuple, utr_tuple

def genomeCoorBed_To_transCoorBed(genomeCoorFileName, transCoorFileName):
    IN = open(genomeCoorFileName)
    OUT = open(transCoorFileName, 'w')
    line = IN.readline()
    while line:
        arr = line.strip().split('\t')
        exonList = norm_exons(arr[7])[0]
        Len = exonList[-1][1]
        exonString = ','.join([str(item[0])+'-'+str(item[1]) for item in exonList])
        transCoorString = arr[5]+'\t'+arr[4]+'\t'+arr[6]+'\t'+str(Len)+'\t'+exonString
        if len(arr) == 9:
            utrList = norm_utr(arr[7], arr[8], arr[3])[1]
            utrString = ','.join([str(item[0])+'-'+str(item[1]) for item in utrList])
            transCoorString += '\t'+utrString
        OUT.writelines( transCoorString+"\n" )
        line = IN.readline()
    OUT.close()



