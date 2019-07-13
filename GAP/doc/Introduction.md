# GAP

<b>G</b>enome <b>A</b>nnotation <b>P</b>aser (<a href="https://github.com/lipan6461188/GAP" style="text-decoration:none;color:blue;">GAP</a>) is a tool for parsing Gencode/Ensembl GTF and NCBI GFF3 genome annotation files.

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
<img src="img/showRNAStructure.png">

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

<img src="img/labelRNAPosition.png">

<style>
.hov:hover {
	background-color: Tomato;
	color: white;
	cursor:pointer;
}

.hov {
	display:inline-block;
	border:1px solid black;
	padding:6px;
	border-radius:8px;
	font-family:Consolas,Monaco,Lucida Console,Liberation Mono,DejaVu Sans Mono, sans-serif;
	float:left;
	margin-left: 10px;
	margin-right: 10px;
	margin-bottom: 10px;
	width: 250px;
	text-align: center;
	text-decoration: none;
	color: darkblue;
}

.clearfix:after {
	content: "";
	clear: both;
	display: table;
}
</style>

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

<div class="clearfix">
<a class="hov" href="#parseGTF">parseGTF.py</a>
<a class="hov" href="#combineNCBIGenome">combineNCBIGenome.py</a>
<a class="hov" href="#GAPinit">GAP.init</a>
<a class="hov" href="#GAPinitGTF">GAP.initGTF</a>
<a class="hov" href="#addSeq">addSeq</a>
</div>

### Get transcript or gene information
<div class="clearfix">
<a class="hov" href="#getTransFeature">getTransFeature</a>
<a class="hov" href="#getGeneParser">getGeneParser</a>
<a class="hov" href="#getGeneIntron">getGeneIntron</a>
<a class="hov" href="#getGeneExon">getGeneExon</a>
<a class="hov" href="#getGeneCombinedIntronExon">getGeneCombinedIntronExon</a>
</div>

### Coordinate transformation

<img src="img/coordinate_system.png">

<img src="img/coordinate_system2.png">

<div class="clearfix">
<a class="hov" href="#geneCoor2genomeCoor">geneCoor2genomeCoor</a>
<a class="hov" href="#genomeCoor2geneCoor">genomeCoor2geneCoor</a>
<a class="hov" href="#genomeCoor2transCoor">genomeCoor2transCoor</a>
<a class="hov" href="#transCoor2genomeCoor">transCoor2genomeCoor</a>
<a class="hov" href="#transCoor2geneCoor">transCoor2geneCoor</a>
<a class="hov" href="#geneCoor2transCoor">geneCoor2transCoor</a>
</div>

### Get transcript list or gene list

<div class="clearfix">
<a class="hov" href="#getTransList">getTransList</a>
<a class="hov" href="#getGeneList">getGeneList</a>
<a class="hov" href="#getmRNATransList">getmRNATransList</a>
<a class="hov" href="#getmRNAGeneList">getmRNAGeneList</a>
<a class="hov" href="#getLenSortedTransDictForGenes">getLenSortedTransDictForGenes</a>
<a class="hov" href="#getTransByGeneID">getTransByGeneID</a>
<a class="hov" href="#getGeneByGeneName">getGeneByGeneName</a>
<a class="hov" href="#getTransByGeneName">getTransByGeneName</a>
<a class="hov" href="#getTransByGeneType">getTransByGeneType</a>
</div>

### Get transcript sequence and mRNA structure

<div class="clearfix">
<a class="hov" href="#getTransSeq">getTransSeq</a>
<a class="hov" href="#showRNAStructure">showRNAStructure</a>
<a class="hov" href="#labelRNAPosition">labelRNAPosition</a>
<a class="hov" href="#getRNAPosition">getRNAPosition</a>
</div>

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

## Demos

### 1. Coordinate transformation

```python
import GAP

hg19_parser = GAP.init("hg19.genomeCoor.bed")

hg19_parser.geneCoor2genomeCoor("ENSG00000227232.5_2", 10)
# ['chr1', 29561, '-']

hg19_parser.geneCoor2transCoor("ENSG00000227232.5_2", 1, 11)
# [['chr1', 29560, 29570, 'ENST00000488147.1_1', 1, 11]]

hg19_parser.transCoor2genomeCoor("ENST00000456328.2_1", 1)
# ['chr1', 11869, '+']

hg19_parser.transCoor2geneCoor("ENST00000456328.2_1", 1, 11)
# [['chr1', 11869, 11879, 'ENSG00000223972.5_2', 1, 11]]

hg19_parser.genomeCoor2geneCoor("chr1", 11869, 11869+10, '+')
# [['chr1', 11869, 11879, 'ENSG00000223972.5_2', 1, 11]]

hg19_parser.genomeCoor2transCoor("chr1", 11869, 11869+10, '+')
# [['chr1', 11869, 11879, 'ENST00000456328.2_1', 1, 11]]

hg19_parser.genomeCoor2transCoor("chr14", 11869, 11869+10, '+')
# A chromosome or TransID/GeneID not in genomefile will raise a KeyError

hg19_parser.transCoor2geneCoor("ENST00000456328.2_1", 1, 1100000)
# out of Gene range will raise CoorFunc.out_of_range
```

### 2. Compare CDS/UTR length of Gencode/RefSeq annotation

```python
import GAP

def statistic_UTR_dist(genomeCoorFile):
    parser = GAP.init(genomeCoorFile)
    utr5_len_dict = {}
    utr3_len_dict = {}
    cds_len_dict = {}
    for transID in parser.getmRNATransList():
        RNA = parser.getTransFeature(transID)
        
        cds_start, cds_end = RNA['cds_start'], RNA['cds_end']
        
        length = RNA['trans_len']
        
        utr5_len = cds_start - 1
        utr3_len = length - cds_end
        cds_len = cds_end - cds_start + 1
        
        gene_name = RNA['gene_name']
        
        if utr5_len > 0 and utr3_len > 0:
            utr5_len_dict[gene_name] = utr5_len_dict.get(gene_name, []) + [utr5_len]
            utr3_len_dict[gene_name] = utr3_len_dict.get(gene_name, []) + [utr3_len]
            cds_len_dict[gene_name] = cds_len_dict.get(gene_name, []) + [cds_len]
    
    return utr5_len_dict, utr3_len_dict, cds_len_dict

hg38_utr5_gencode, hg38_utr3_gencode, hg38_cds_gencode = statistic_UTR_dist("hg38_gencode.genomeCoor.bed")
hg38_utr5_ncbi,    hg38_utr3_ncbi,    hg38_cds_ncbi    = statistic_UTR_dist("hg38_ncbi.genomeCoor.bed")

mm10_utr5_gencode, mm10_utr3_gencode, mm10_cds_gencode = statistic_UTR_dist("mm10_gencode.genomeCoor.bed")
mm10_utr5_ncbi,    mm10_utr3_ncbi,    mm10_cds_ncbi    = statistic_UTR_dist("mm10_ncbi.genomeCoor.bed")

def plot_single_anno(UTR5_Dict, UTR3_Dict, CDS_Dict):
    import seaborn as sns
    import pandas as pd
    import matplotlib.pyplot as plt
    
    dataSet = []
    for gene_name in UTR5_Dict:
        dataSet += [ [item, "UTR5"] for item in UTR5_Dict[gene_name] ]
    
    for gene_name in UTR3_Dict:
        dataSet += [ [item, "UTR3"] for item in UTR3_Dict[gene_name] ]
    
    for gene_name in CDS_Dict:
        dataSet += [ [item, "CDS"] for item in CDS_Dict[gene_name] ]
    
    dataSet = pd.DataFrame(dataSet, columns=['length', 'type'])
    
    # show some statistics
    print "median of UTR5: ", dataSet.loc[dataSet.type=="UTR5", "length"].median()
    print "median of UTR3: ", dataSet.loc[dataSet.type=="UTR3", "length"].median()
    print "median of CDS: ", dataSet.loc[dataSet.type=="CDS", "length"].median()
    
    # plot
    sns.violinplot(data=dataSet, x='type', y='length')
    plt.xlabel("mRNA Part", fontsize="large")
    plt.ylabel("Length Distribution", fontsize="large")
    plt.ylim(-2000,5000)
    plt.show()

plot_single_anno(hg38_utr5_gencode, hg38_utr3_gencode, hg38_cds_gencode)
plot_single_anno(hg38_utr5_ncbi, hg38_utr3_ncbi, hg38_cds_ncbi)

plot_single_anno(mm10_utr5_gencode, mm10_utr3_gencode, mm10_cds_gencode)
plot_single_anno(mm10_utr5_ncbi, mm10_utr3_ncbi, mm10_cds_ncbi)

def compair_anno(Gencode, RefSeq, title):
    import numpy as np
    import seaborn as sns
    import pandas as pd
    import matplotlib.pyplot as plt
    
    dataSet = []
    for gene_name in set(Gencode)&set(RefSeq):
        dataSet.append( [np.median(Gencode[gene_name]), "Gencode"] )
        dataSet.append( [np.median(RefSeq[gene_name]), "RefSeq"] )
    
    dataSet = pd.DataFrame(dataSet, columns=['length', 'type'])
    
    # show some statistics
    print "median of Gencode: ", dataSet.loc[dataSet.type=="Gencode", "length"].median()
    print "median of RefSeq: ", dataSet.loc[dataSet.type=="RefSeq", "length"].median()
    
    # plot
    sns.violinplot(data=dataSet, x='type', y='length')
    plt.xlabel(title, fontsize="large")
    plt.ylabel("Length Distribution", fontsize="large")
    plt.ylim(-500,2000)
    plt.show()


compair_anno(hg38_utr5_gencode, hg38_utr5_ncbi, "5UTR")
compair_anno(hg38_utr3_gencode, hg38_utr3_ncbi, "3UTR")
compair_anno(hg38_cds_gencode,  hg38_cds_ncbi,  "CDS")
```

## Authors

* **Li Pan** - *Programmer* - [Zhanglab](http://zhanglab.life.tsinghua.edu.cn)


