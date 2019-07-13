# GAP

<b>G</b>enome <b>A</b>nnotation <b>P</b>aser (GAP) is a tool for parsing Gencode/Ensembl GTF and NCBI GFF3 genome annotation files.

## Installation

### Prerequisites

Python 2.7 or Python 3

pysam: https://pypi.org/project/pysam/


### Installing

Add these lines to your `~/.bash_profile`, `[GAP]` is the <b>absolute</b> path of GAP.

```bash
export PYTHONPATH=[GAP]:$PYTHONPATH
export PATH=[GAP]:$PATH
```

## Main functions

* Fetch transcritome file from genome with GFF/GTF annotation;
* Conversion between three coordinate systems: Gene/Transcript/Genome;
* Get transcription information (chrID, gene type, length, CDS region...) given a transcription ID or gene name;
* Show mRNA structure
* ...

## Download data

> Gencode GTF/Genome: https://www.gencodegenes.org <br>
> Enesembl GTF/Genome: https://ensembl.org/info/data/ftp/index.html <br>
> NCBI GFF3/Genome: ftp://ftp.ncbi.nlm.nih.gov/genomes

## Quick start

#### 1. Convert a GTF file to *.genomeCoor.bed file
```shell
parseGTF.py -g chicken.gtf -o chicken -s ensembl --genome chicken_ensembl.fa
```
The genome file `--genome chicken_ensembl.fa` is optional.

It will produce three files: 

* `chicken.genomeCoor.bed`    -- A simple version of the genome-based annotation file
* `chicken.transCoor.bed`		-- A simple version of the transcript-based annotation file
* `chicken_transcriptome.fa`  -- Transcriptome file

#### 2. Read the annotation

```python
import GAP
chicken_parser = GAP.init("chicken.genomeCoor.bed", "chicken_transcriptome.fa")
```

Another way is to read the GTF file directly.

```python
chicken_parser = GAP.initGTF("chicken.gtf", genomeFile="chicken_ensembl.fa", source='Ensembl')
```

#### 3. Get all mRNAs

```python
mRNAs = chicken_parser.getmRNATransList(); print len(mRNAs)
# 30252
print mRNAs[:5]
# ['ENSGALT00000005443', 'ENSGALT00000005447', 'ENSGALT00000013873', 'ENSGALT00000001353', 'ENSGALT00000001352']
```
It shows that chicken has 30252 mRNA transcripts.

#### 4. Get GAPDH gene id and transcripts

```python
GAPDH_gID = chicken_parser.getGeneByGeneName("GAPDH"); print GAPDH_gID
# ENSGALG00000014442

GAPDH_transcripts = chicken_parser.getTransByGeneID(GAPDH_gID); print GAPDH_transcripts
# ['ENSGALT00000046744', 'ENSGALT00000086833', 'ENSGALT00000023323', 'ENSGALT00000086032', 'ENSGALT00000074237', 'ENSGALT00000090208', 'ENSGALT00000051222', 'ENSGALT00000054080', 'ENSGALT00000089752', 'ENSGALT00000085687']

## Print all transcript gene type and length
for tid in GAPDH_transcripts:
	trans_feature = chicken_parser.getTransFeature(tid)
	print tid, trans_feature['gene_type'], trans_feature['trans_len']

# ENSGALT00000046744 protein_coding 1076
# ENSGALT00000086833 protein_coding 1122
# ENSGALT00000023323 protein_coding 1288
# ENSGALT00000086032 protein_coding 1302
# ENSGALT00000074237 protein_coding 1091
# ENSGALT00000090208 protein_coding 670
# ENSGALT00000051222 protein_coding 1179
# ENSGALT00000054080 protein_coding 1498
# ENSGALT00000089752 protein_coding 434
# ENSGALT00000085687 protein_coding 991
```

#### 5. Get UTR and CDS region of longest transcript of GAPDH

```python
## Method 1
ft = chicken_parser.getTransFeature("ENSGALT00000054080")
cds_start, cds_end = ft['cds_start'], ft['cds_end']
GAPDH = chicken_parser.getTransSeq("ENSGALT00000054080")

UTR_5 = GAPDH[:cds_start-1]
CDS = GAPDH[cds_start-1:cds_end]
UTR_3 = GAPDH[cds_end:]

## Method 2
print chicken_parser.showRNAStructure("ENSGALT00000054080")
```
<img src="doc/img/showRNAStructure.png">

#### 6. Get genome coordination of GAPDH start codon

```python
chrID, chrPos, strand = chicken_parser.transCoor2genomeCoor("ENSGALT00000054080", cds_start)
print chrID, chrPos, strand
# ['1', 76953317, '+']
```
It shows that GAPDH start codon is located at 76953317 of positive strand of chromosome 1

#### 7. Find and label all GGAC motifs in GAPDH

```python
## Get Sequence
seq = chicken_parser.getTransSeq("ENSGALT00000054080")

## Collect motif sites
locs = []
start = 0
while 1:
    start = seq.find("GGAC", start+1)
    if start == -1: break
    locs.append(start)

## Label
for loc in locs:
    print str(loc)+"\t"+chicken_parser.labelRNAPosition("ENSGALT00000054080", [loc,loc])
```

<img src="doc/img/labelRNAPosition.png">

## Subjects

### Init a GAP object

> There are two ways to init a GAP object: <br />
> 1. Convert the GTF or GFF3 file to a *.genomeCoor.bed file and read it;<br />
> 2. Read the GTF or GFF3 directly.

#### Method 1 -- Parse GTF/GFF3 file

```bash
parseGTF.py -g chicken.gtf -o chicken -s ensembl --genome chicken_ensembl.fa
```
The `--genome` is optional.

```python
import GAP
chicken_parser = GAP.init("chicken.genomeCoor.bed", "chicken_transcriptome.fa")
```

#### Method 2 -- Read GTF/GFF3 directly

```python
import GAP
chicken_parser = GAP.initGTF("chicken.gtf", genomeFile="chicken_ensembl.fa", source='Ensembl')
```

<ul>
<li> <a href="#parseGTF">parseGTF.py</a> </li>
<li> <a href="#combineNCBIGenome">combineNCBIGenome.py</a> </li>
<li> <a href="#GAPinit">GAP.init</a> </li>
<li> <a href="#GAPinitGTF">GAP.initGTF</a> </li>
<li> <a href="#addSeq">addSeq</a> </li>
</ul>

### Get transcript or gene information
<ul>
<li> <a href="#getTransFeature">getTransFeature</a> </li>
<li> <a href="#getGeneParser">getGeneParser</a> </li>
<li> <a href="#getGeneIntron">getGeneIntron</a> </li>
<li> <a href="#getGeneExon">getGeneExon</a> </li>
<li> <a href="#getGeneCombinedIntronExon">getGeneCombinedIntronExon</a> </li>
</ul>

### Coordinate transformation

<img src="doc/img/coordinate_system.png">

<img src="doc/img/coordinate_system2.png">

<ul>
<li> <a href="#geneCoor2genomeCoor">geneCoor2genomeCoor</a> </li>
<li> <a href="#genomeCoor2geneCoor">genomeCoor2geneCoor</a> </li>
<li> <a href="#genomeCoor2transCoor">genomeCoor2transCoor</a> </li>
<li> <a href="#transCoor2genomeCoor">transCoor2genomeCoor</a> </li>
<li> <a href="#transCoor2geneCoor">transCoor2geneCoor</a> </li>
<li> <a href="#geneCoor2transCoor">geneCoor2transCoor</a> </li>
</ul>

### Get transcript list or gene list

<ul>
<li> <a href="#getTransList">getTransList</a> </li>
<li> <a href="#getGeneList">getGeneList</a> </li>
<li> <a href="#getmRNATransList">getmRNATransList</a> </li>
<li> <a href="#getmRNAGeneList">getmRNAGeneList</a> </li>
<li> <a href="#getLenSortedTransDictForGenes">getLenSortedTransDictForGenes</a> </li>
<li> <a href="#getTransByGeneID">getTransByGeneID</a> </li>
<li> <a href="#getGeneByGeneName">getGeneByGeneName</a> </li>
<li> <a href="#getTransByGeneName">getTransByGeneName</a> </li>
<li> <a href="#getTransByGeneType">getTransByGeneType</a> </li>
</ul>

### Get transcript sequence and mRNA structure

<ul>
<li> <a href="#getTransSeq">getTransSeq</a> </li>
<li> <a href="#showRNAStructure">showRNAStructure</a> </li>
<li> <a href="#labelRNAPosition">labelRNAPosition</a> </li>
<li> <a href="#getRNAPosition">getRNAPosition</a> </li>
</ul>

## Excutable files

<h4 id="parseGTF"> <span style="color:DarkGreen;font-family:menlo, sans-serif">parseGTF.py</h4>

##### <span style="color:darkred">Features</span>
<pre>
Parse GTF/GFF3 file and produce *.genomeCoor.bed, *.transCoor.bed and *_transcriptome.fa files
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
-g         genome annotation, Gencode/Ensembl GTF file or NCBI GFF3
-o         output file prefix
-s         [ensembl|gencode|ncbi] data source

optional:
--genome   fetch transcriptome from genome file, produce a prefix_transcriptome.fa file
--noscaffold	Remove scaffolds. Scaffolds are defined as those chromosomes with id length > 6 and not startswith chr and NC_
--rawchr	Use raw chromosome ID, don't convert NC_* to chrXX. Only useful for NCBI source GFF3
</pre>

##### <span style="color:darkred">Output</span>
<pre>
*.genomeCoor.bed, *.transCoor.bed and *_transcriptome.fa
</pre>

<h4 id="combineNCBIGenome"> <span style="color:DarkGreen;font-family:menlo, sans-serif">combineNCBIGenome.py</h4>

##### <span style="color:darkred">Features</span>
<pre>
Combine single chromosomes download from NCBI. Convert the NC_* code to chr* code
NC_006088.5 => chr1
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
./combineNCBIGenome.py chr1.fasta chr2.fasta chr3.fasta... > outFile.fa
</pre>

##### <span style="color:darkred">Return</span>
<pre>
A combined genome file with chr* code
</pre>

## Functions

<h4 id="GAPinit"> <span style="color:DarkGreen;font-family:menlo, sans-serif">GAP.init(genomeCoorBedFile, seqFn="", showAttr=True, rem_tVersion=False, rem_gVersion=False)</span> </h4>

##### <span style="color:darkred">Features</span>
<pre>
Init a GAP object from *.genomeCoor.bed file
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
genomeCoorBedFile   -- A *.genomeCoor.bed file produced by parseGTF.py
seqFn               -- Transcriptome fasta file produced by parseGTF.py
showAttr            -- Show an example
rem_tVersion        -- Remove version information. ENST000000022311.2 => ENST000000022311
rem_gVersion        -- Remove version information. ENSG000000022311.2 => ENSG000000022311
</pre>

##### <span style="color:darkred">Return</span>
<pre>
GAP object
</pre>


<h4 id="GAPinitGTF"> <span style="color:DarkGreen;font-family:menlo, sans-serif">initGTF(AnnotationGTF, source, genomeFile="", showAttr=True, rem_tVersion=False, rem_gVersion=False, verbose=False)</span> </h4>


##### <span style="color:darkred">Features</span>
<pre>
Init a GAP object from GTF/GFF3 file
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
AnnotationGTF       -- Ensembl/Gencode GTF file or NCBI GFF3 file
genomeFile          -- Genome file
source              -- Gencode/Ensembl/NCBI
showAttr            -- Show an example
rem_tVersion        -- Remove version information. ENST000000022311.2 => ENST000000022311
rem_gVersion        -- Remove version information. ENSG000000022311.2 => ENSG000000022311
verbose             -- Show process information
</pre>

##### <span style="color:darkred">Return</span>
<pre>
GAP object
</pre>


<h4 id="addSeq"> <span style="color:DarkGreen;font-family:menlo, sans-serif">addSeq(seqFileName, remove_tid_version=False)</span> </h4>

##### <span style="color:darkred">Features</span>
<pre>
If the sequence is not provided when init a GAP object, add the sequence to it.
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
seqFileName         -- Transcriptome fasta file produced by parseGTF.py
remove_tid_version  -- Remove version information. ENST000000022311.2 => ENST000000022311
</pre>

##### <span style="color:darkred">Return</span>
<pre>
No
</pre>

<h4 id="getTransFeature"> <span style="color:DarkGreen;font-family:menlo, sans-serif">getTransFeature(transID, showAttr=False, verbose=True)</span> </h4>

##### <span style="color:darkred">Features</span>
<pre>
Get transcript features given the transcript id
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
transID             -- Transcript ID
showAttr            -- Show an example
verbose             -- Print the warning information when transcript not found
</pre>

##### <span style="color:darkred">Return</span>
<pre>
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
</pre>

<h4 id="getGeneParser"> <span style="color:DarkGreen;font-family:menlo, sans-serif">getGeneParser(showAttr=True)</span> </h4>

##### <span style="color:darkred">Features</span>
<pre>
Get gene parser with gene informations
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
showAttr            -- Show an example
</pre>

##### <span style="color:darkred">Return</span>
<pre>
Return a list of dictionaries:
    { geneID => { ... }, ... } includes
        chr         -- Genome chromosome
        strand      -- + or -
        start       -- Genome start
        end         -- Genome end
        
        gene_name   -- Gene symbol
        gene_type   -- All transcript types
        
        length      -- Gene length (end-start+1)
        
        transcript  -- All transcripts belong to this gene
</pre>

<h4 id="getGeneIntron"> <span style="color:DarkGreen;font-family:menlo, sans-serif">getGeneIntron(geneID)</span> </h4>

##### <span style="color:darkred">Features</span>
<pre>
Get genome-based intron corrdinates of given gene
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
geneID              -- Gene id
</pre>

##### <span style="color:darkred">Return</span>
<pre>
Get introns of all transcripts of this gene:
            { transID => [[intron1_start, intron1_end], [intron2_start, intron2_end], [intron3_start, intron3_end]],... }
</pre>

<h4 id="getGeneExon"> <span style="color:DarkGreen;font-family:menlo, sans-serif">getGeneExon(geneID)</span> </h4>

##### <span style="color:darkred">Features</span>
<pre>
Get genome-based exon corrdinates of given gene
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
geneID              -- Gene id
</pre>

##### <span style="color:darkred">Return</span>
<pre>
Get exons of all transcripts of this gene:
            { transID => [[exon1_start, exon1_end], [exon2_start, exon2_end], [exon3_start, exon3_end]...],... }
</pre>

<h4 id="geneCoor2genomeCoor"> <span style="color:DarkGreen;font-family:menlo, sans-serif">geneCoor2genomeCoor(geneID, pos)</span> </h4>

##### <span style="color:darkred">Features</span>
<pre>
Convert gene-based coordinates to genome-based coordinates
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
geneID              -- Gene id
pos                 -- Gene-based position
</pre>

##### <span style="color:darkred">Return</span>
<pre>
Return: [chrID, chrPos, Strand]
</pre>

<h4 id="genomeCoor2geneCoor"> <span style="color:DarkGreen;font-family:menlo, sans-serif">genomeCoor2geneCoor(chrID, start, end, strand)</span> </h4>

##### <span style="color:darkred">Features</span>
<pre>
Convert genome-based coordinates to gene-based coordinates
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
chrID               -- Chromosome id
start               -- Chromosome-based start position
end                 -- Chromosome-based end position
strand              -- Chromosome strand
</pre>

##### <span style="color:darkred">Return</span>
<pre>
Return [ [chrID, chrStart, chrEnd, geneID, geneStart, geneEnd], ... ]
</pre>

<h4 id="genomeCoor2transCoor"> <span style="color:DarkGreen;font-family:menlo, sans-serif">genomeCoor2transCoor(chrID, start, end, strand)</span> </h4>

##### <span style="color:darkred">Features</span>
<pre>
Convert genome-based coordinates to transcript-based coordinates
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
chrID               -- Chromosome id
start               -- Chromosome-based start position
end                 -- Chromosome-based end position
strand              -- Chromosome strand
</pre>

##### <span style="color:darkred">Return</span>
<pre>
Return [ [chrID, chrStart, chrEnd, transID, transStart, transEnd], ... ]
</pre>

<h4 id="transCoor2genomeCoor"> <span style="color:DarkGreen;font-family:menlo, sans-serif">transCoor2genomeCoor(transID, pos)</span> </h4>

##### <span style="color:darkred">Features</span>
<pre>
Convert transcript-based coordinates to genome-based coordinates
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
transID             -- Transcript id
pos                 -- Transcript-based position
</pre>

##### <span style="color:darkred">Return</span>
<pre>
Return [chrID, chrPos, Strand]
</pre>

<h4 id="geneCoor2transCoor"> <span style="color:DarkGreen;font-family:menlo, sans-serif">geneCoor2transCoor(geneID, start, end)</span> </h4>

##### <span style="color:darkred">Features</span>
<pre>
Convert gene-based coordinates to transcript-based coordinates
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
geneID              -- Gene id
start               -- Gene-based start position
end                 -- Gene-based end position
</pre>

##### <span style="color:darkred">Return</span>
<pre>
Return  [chrID, chrStart, chrEnd, transID, transStart, transEnd], ... ]
</pre>

<h4 id="transCoor2geneCoor"> <span style="color:DarkGreen;font-family:menlo, sans-serif">transCoor2geneCoor(transID, start, end)</span> </h4>

##### <span style="color:darkred">Features</span>
<pre>
Convert transcript-based coordinates to gene-based coordinates
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
transID             -- Transcript id
start               -- Transcript-based start position
end                 -- Transcript-based end position
</pre>

##### <span style="color:darkred">Return</span>
<pre>
return [ [chrID, chrStart, chrEnd, geneID, geneStart, geneEnd], ... ]
</pre>

<h4 id="getTransList"> <span style="color:DarkGreen;font-family:menlo, sans-serif">getTransList()</span> </h4>

##### <span style="color:darkred">Features</span>
<pre>
Return list of all transcript id
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
No
</pre>

##### <span style="color:darkred">Return</span>
<pre>
return [ transID1, transID2, ... ]
</pre>

<h4 id="getGeneList"> <span style="color:DarkGreen;font-family:menlo, sans-serif">getGeneList()</span> </h4>

##### <span style="color:darkred">Features</span>
<pre>
Return list of all gene id
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
No
</pre>

##### <span style="color:darkred">Return</span>
<pre>
return [ geneID1, geneID2, ... ]
</pre>

<h4 id="getmRNATransList"> <span style="color:DarkGreen;font-family:menlo, sans-serif">getmRNATransList()</span> </h4>

##### <span style="color:darkred">Features</span>
<pre>
Return list of all mRNA id
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
No
</pre>

##### <span style="color:darkred">Return</span>
<pre>
return [ mRNAID1, mRNAID2, ... ]
</pre>

<h4 id="getmRNAGeneList"> <span style="color:DarkGreen;font-family:menlo, sans-serif">getmRNAGeneList()</span> </h4>

##### <span style="color:darkred">Features</span>
<pre>
Return list of all mRNA gene id
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
No
</pre>

##### <span style="color:darkred">Return</span>
<pre>
return [ mRNAGeneID1, mRNAGeneID2, ... ]
</pre>

<h4 id="getLenSortedTransDictForGenes"> <span style="color:DarkGreen;font-family:menlo, sans-serif">getLenSortedTransDictForGenes(only=[])</span> </h4>

##### <span style="color:darkred">Features</span>
<pre>
Return a dictionary of geneID to sorted transID list
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
only                --  Contraint transcript types, such as mRNA, snRNA ...
</pre>

##### <span style="color:darkred">Return</span>
<pre>
return { geneID => [ transID1, transID2, ... ], ...}
transcripts are sorted by length
</pre>

<h4 id="getTransByGeneID"> <span style="color:DarkGreen;font-family:menlo, sans-serif">getTransByGeneID(geneID)</span> </h4>

##### <span style="color:darkred">Features</span>
<pre>
Return transcripts belong to specific gene
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
geneID              -- Gene id
</pre>

##### <span style="color:darkred">Return</span>
<pre>
return [ transID1, transID2... ]
</pre>

<h4 id="getGeneByGeneName"> <span style="color:DarkGreen;font-family:menlo, sans-serif">getGeneByGeneName(geneName)</span> </h4>

##### <span style="color:darkred">Features</span>
<pre>
Return gene id with specific gene name
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
geneName            -- Gene symbol
</pre>

##### <span style="color:darkred">Return</span>
<pre>
return geneID
</pre>

<h4 id="getTransByGeneName"> <span style="color:DarkGreen;font-family:menlo, sans-serif">getTransByGeneName(geneName)</span> </h4>

##### <span style="color:darkred">Features</span>
<pre>
Return transcripts belong to specific gene
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
geneName            -- Gene symbol
</pre>

##### <span style="color:darkred">Return</span>
<pre>
return [ transID1, transID2... ]
</pre>

<h4 id="getTransByGeneType"> <span style="color:DarkGreen;font-family:menlo, sans-serif">getTransByGeneType(geneType)</span> </h4>

##### <span style="color:darkred">Features</span>
<pre>
 Return a list of transcripts belong to specific gene type
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
geneType            -- Gene type
</pre>

##### <span style="color:darkred">Return</span>
<pre>
return [ transID1, transID2... ]
</pre>

<h4 id="getTransSeq"> <span style="color:DarkGreen;font-family:menlo, sans-serif">getTransSeq(transID)</span> </h4>

##### <span style="color:darkred">Features</span>
<pre>
 Return transcript sequence
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
transID             -- Transcript id
</pre>

##### <span style="color:darkred">Return</span>
<pre>
return sequence
</pre>

<h4 id="getGeneCombinedIntronExon"> <span style="color:DarkGreen;font-family:menlo, sans-serif">getGeneCombinedIntronExon(geneID, verbose=True)</span> </h4>

##### <span style="color:darkred">Features</span>
<pre>
Parse gene intron/exon regions:
    1. Combine exons from all transcripts, they are defined as exon regions
    2. Remove the exon regions from gene regions, they are defined as intron regions
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
geneID          -- Gene id
verbose         -- Print error information when error occured
</pre>

##### <span style="color:darkred">Return</span>
<pre>
return: [intron_regions, exon_regions]
</pre>

<h4 id="showRNAStructure"> <span style="color:DarkGreen;font-family:menlo, sans-serif">showRNAStructure(transID)</span> </h4>

##### <span style="color:darkred">Features</span>
<pre>
Return string of mRNA sequence with color labeled its UTR/CDS/Codons
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
transID         -- Transcript id
</pre>

##### <span style="color:darkred">Return</span>
<pre>
return: colored sequence
</pre>

<h4 id="labelRNAPosition"> <span style="color:DarkGreen;font-family:menlo, sans-serif">labelRNAPosition(transID, region, bn=None, bw=None)</span> </h4>

##### <span style="color:darkred">Features</span>
<pre>
Return string of mRNA sequence with color labeled its UTR/CDS/Codons

        return RNA structure and label the region
        -------|||||||||||||||||--------------------------
         5'UTR        CDS                3'UTR
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
transID         -- Transcript id
region          -- [start, end]
bn              -- Bin number (default: 50)
bw              -- Bin width
bw and bn cannot be specified at the same time
</pre>

##### <span style="color:darkred">Return</span>
<pre>
return: -------|||||||||||||||||--------------------------
</pre>

<h4 id="getRNAPosition"> <span style="color:DarkGreen;font-family:menlo, sans-serif">getRNAPosition(transID, region)</span> </h4>

##### <span style="color:darkred">Features</span>
<pre>
Get the position of region located in mRNA
</pre>

##### <span style="color:darkred">Parameters</span>
<pre>
transID         -- Transcript id
region          -- [start, end]
</pre>

##### <span style="color:darkred">Return</span>
<pre>
Return any one of:
    -> not_mRNA
    -> 5UTR
    -> span_5UTR_CDS
    -> CDS
    -> span_CDS_3UTR
    -> 3UTR
    -> span_5UTR_CDS_3UTR
    -> INVALID
</pre>

## Notice

Complete documentation is in doc/Introduction.html


## Authors

* **Li Pan** - *Programmer* - [Zhanglab](http://zhanglab.life.tsinghua.edu.cn)


