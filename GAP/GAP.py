#-*- coding:utf-8 -*-

def init(genomeCoorBedFile, seqFn='', showAttr=True, rem_tVersion=False, rem_gVersion=False):
    """
    genomeCoorBedFile   -- A *.genomeCoor.bed file produced by parseGTF.py
    seqFn               -- Transcriptome fasta file produced by parseGTF.py
    showAttr            -- Show an example
    rem_tVersion        -- Remove version information. ENST000000022311.2 => ENST000000022311
    rem_gVersion        -- Remove version information. ENSG000000022311.2 => ENSG000000022311
    
    Return a ParseTransClass object to parse genome
    """
    import ParseTrans
    return ParseTrans.ParseTransClass(genomeCoorBedFile, seqFileName=seqFn, showAttr=showAttr, 
                                        remove_tid_version=rem_tVersion, remove_gid_version=rem_gVersion)

def initGTF(AnnotationGTF, source, genomeFile='', showAttr=True, rem_tVersion=False, rem_gVersion=False, rem_scaffold=False, raw_chrID=False, verbose=False):
    """
    AnnotationGTF       -- Ensembl/Gencode GTF file or NCBI GFF3 file
    genomeFile          -- Genome file
    source              -- Gencode/Ensembl/NCBI
    showAttr            -- Show an example
    rem_tVersion        -- Remove version information. ENST000000022311.2 => ENST000000022311
    rem_gVersion        -- Remove version information. ENSG000000022311.2 => ENSG000000022311
    rem_scaffold        -- Remove scaffolds. scaffolds are defined as 
                           those chromosomes with id length > 6 and not startswith chr and NC_
    --rawchr            -- Use raw chromosome ID, don't convert NC_* to chrXX. Only useful for NCBI source GFF3
    verbose             -- Show process information
    
    Return a ParseTransClass object to parse genome
    """
    import GTFParserFunc
    
    if source == 'Gencode' or source == 'Ensembl':
        handle = GTFParserFunc.read_ensembl_gtf(AnnotationGTF, rem_scaffold=rem_scaffold)
        Parser = __build_parser(handle, genomeFn=genomeFile, source=source, showAttr=showAttr, rem_tVersion=rem_tVersion, rem_gVersion=rem_gVersion, verbose=verbose)
    elif source == 'NCBI':
        handle = GTFParserFunc.read_ncbi_gff3(AnnotationGTF, rem_scaffold=rem_scaffold, raw_chrID=raw_chrID)
        Parser = __build_parser(handle, genomeFn=genomeFile, source=source, showAttr=showAttr, rem_tVersion=rem_tVersion, rem_gVersion=rem_gVersion, verbose=verbose)
    else:
        sys.stderr.writelines("Error: source must be Gencode or NCBI\n")
    
    return Parser


def __build_parser(GenomeHandler, genomeFn='', source='Gencode', showAttr=True, rem_tVersion=False, rem_gVersion=False, verbose=False):
    import random, os, ParseTrans
    import GTFParserFunc
    
    rID = random.randint(10000, 99999)
    tmp_genomeCoor_file = "/tmp/tmp_%s_genomeCoor.bed" % (rID, )
    
    if source == 'Gencode' or source == 'Ensembl':
        GTFParserFunc.write_gtf_genomeCoor_bed(GenomeHandler, tmp_genomeCoor_file)
    elif source == 'NCBI':
        GTFParserFunc.write_gff3_genomeCoor_bed(GenomeHandler, tmp_genomeCoor_file, raw_chrID=raw_chrID)
    else:
        sys.stderr.writelines("Error!\n")
        return None
    
    tmp_transcriptome_file = ""
    if genomeFn:
        tmp_transcriptome_file = "/tmp/tmp_%s_transcriptome.fasta" % (rID, )
        GTFParserFunc.writeTranscriptome(tmp_genomeCoor_file, genomeFn, tmp_transcriptome_file, verbose=verbose, showAttr=showAttr)
    
    Parser = ParseTrans.ParseTransClass(tmp_genomeCoor_file, seqFileName=tmp_transcriptome_file, showAttr=showAttr, remove_tid_version=rem_tVersion, remove_gid_version=rem_gVersion)
    
    os.remove(tmp_genomeCoor_file)
    if tmp_transcriptome_file:
        os.remove(tmp_transcriptome_file)
    
    return Parser

