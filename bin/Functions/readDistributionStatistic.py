#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import numpy
import getopt
import os
import version
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from multiprocessing import Pool


"""

F1="/Share/home/zhangqf8/zhuoer/raw/GSE45803/SRR815612.fastq.gz"
F2="/Share/home/zhangqf8/test/Structure_Database/1.uniqReads/yeast_vivo_rep1.fq"
F3="/Share/home/zhangqf8/test/Structure_Database/2.trimReads/yeast_vivo_rep1.fq.gz"
F4="/Share/home/zhangqf8/test/Structure_Database/3.rem_rRNA/yeast_vivo_rep1.fq"
F5="/Share/home/zhangqf8/test/Structure_Database/4.map_smallRNA/yeast_vivo_rep1.fq"
F6="/Share/home/zhangqf8/test/Structure_Database/5.map_genome/yeast_vivo_rep1.sorted.bam"

N1="/Share/home/zhangqf8/zhuoer/raw/GSE45803/SRR815613.fastq.gz"
N2="/Share/home/zhangqf8/test/Structure_Database/1.uniqReads/yeast_vivo_rep2.fq"
N3="/Share/home/zhangqf8/test/Structure_Database/2.trimReads/yeast_vivo_rep2.fq.gz"
N4="/Share/home/zhangqf8/test/Structure_Database/3.rem_rRNA/yeast_vivo_rep2.fq"
N5="/Share/home/zhangqf8/test/Structure_Database/4.map_smallRNA/yeast_vivo_rep2.fq"
N6="/Share/home/zhangqf8/test/Structure_Database/5.map_genome/yeast_vivo_rep2.sorted.bam"

Q1="/Share/home/zhangqf8/zhuoer/raw/GSE45803/SRR815614.fastq.gz"
Q2="/Share/home/zhangqf8/test/Structure_Database/1.uniqReads/yeast_vivo_rep3.fq"
Q3="/Share/home/zhangqf8/test/Structure_Database/2.trimReads/yeast_vivo_rep3.fq.gz"
Q4="/Share/home/zhangqf8/test/Structure_Database/3.rem_rRNA/yeast_vivo_rep3.fq"
Q5="/Share/home/zhangqf8/test/Structure_Database/4.map_smallRNA/yeast_vivo_rep3.fq"
Q6="/Share/home/zhangqf8/test/Structure_Database/5.map_genome/yeast_vivo_rep3.sorted.bam"


sys.argv = ['bin', '-@', ",".join([F1,F2,F3,F4,F5,F6]), \
    '-@', ",".join([N1,N2,N3,N4,N5,N6]), '-@', ",".join([Q1,Q2,Q3,Q4,Q5,Q6]), \
    '--processtag', \
    'raw,unique,trimmed,map_rRNA,smallRNA,mapGenome', '--sampletag', 'yeast_vivo_rep1,yeast_vivo_rep2,yeast_vivo_rep3', \
    '-o', '/Share/home/zhangqf8/test/Structure_Database/test']


readDistributionStatistic \
    -@ raw1,p1,p2,p3,p4 \
    -@ raw1,p1,p2,p3,p4 \
    -@ raw1,p1,p2,p3,p4 \
    --sampletag vivo1,vivo2,vivo3 \
    --processtag raw,collapse,trim,remove_rRNA,final
"""



Usage = """
readDistributionStatistic - Count fastq or sam to statistic reads distribution
==============================================================================
    This script counts the number of reads remaining in each step of one 
or more of the original fastq files after a series of processing, and finally 
generates a report. Multiple -@ are allowed in the parameter, each of which 
is a series of processing of the original fastq file.
    -@ raw_fastq,rest_of_reads_step1,rest_of_reads_step2,...,final_map_sam
    All rest_of_reads_stepX can be sam/bam/fastq file, but only final_map_sam can 
contains un-mapped reads
    --processtag gives a tag name for each process procedure
    --sampletag gives a sample name for each -@ sample

\x1b[1mUSAGE:\x1b[0m
  %s -@ ... -@ ... --sampletag ... --processtag ... -o outFolder

\x1b[1mHELP:\x1b[0m
  -@                <String,String...>
                        fastq/sam/bam files separated by comma
  --sampletag       <String,String...>
                        tag for each file
  --processtag      <String,String...>
                        tag for each process
  -o                <String>
                        A directory to output the report

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
    params = { 
        'outFolder': None,
        'sampleTags': None,
        'processTags': None,
        'filepipe': []
    }
    opts, args = getopt.getopt(sys.argv[1:], 'h@:o:', ['sampletag=','processtag='])
    for op, value in opts:
        if op == '-h':
            print(Usage)
            sys.exit(-1)
        # Basic Parameters
        elif op == '-@':
            files = value.strip().split(',')
            params['filepipe'].append(files)
        elif op == '-o':
            params['outFolder'] = os.path.abspath(value).rstrip('/') + '/'
        elif op == '--sampletag':
            tags = value.strip().split(',')
            params['sampleTags'] = tags
        elif op == '--processtag':
            tags = value.strip().split(',')
            params['processTags'] = tags
        
        else:
            sys.stderr.writelines( "Error: unrecognized parameter: "+op+"\n" )
            print(Usage)
            sys.exit(-1)
    
    # check
    if not params['outFolder'] or not params['sampleTags'] or not params['processTags'] or not params['filepipe']:
        sys.stderr.writelines( "Error: Specify -@ --sampletag --processtag -o\n" )
        print(Usage)
        exit(-1)
    
    fileLen = len(params['sampleTags'])
    if fileLen != len(params['filepipe']):
        sys.stderr.writelines( "Error: different number of -@ (%d) and --sampletag (%d) \n" % (fileLen, len(params['filepipe'])) )
        print(Usage)
        exit(-1)
    
    processLen = len(params['processTags'])
    if processLen != len(params['filepipe'][0]):
        sys.stderr.writelines( "Error: different number of -@ files (%d) and --processTags (%d) \n" % (processLen, len(params['filepipe'][0])) )
        print(Usage)
        exit(-1)
    
    for file in params['filepipe']:
        if processLen != len(file):
            sys.stderr.writelines( "Error: different number of -@ files (%d) and --processTags (%d) \n" % (processLen, len(file)) )
            print(Usage)
            exit(-1)
    
    return params

def count_fq(inFile):
    
    if inFile.endswith('.gz'):
        count = getoutput("gzip -d %s -c | wc -l" % (inFile, ))
        count = int(count)
    else:
        count = getoutput("wc -l %s | awk '{print $1}'" % (inFile, ))
        count = int(count)
    
    assert count%4 == 0
    return count/4

def count_sam(inFile):
    
    CMD = "samtools view %s | awk 'BEGIN{map=0;unmap=0;}{if($3==\"*\"){unmap+=1}else{map+=1}}END{print map\"\t\"unmap}'"
    
    mapcount, unmapcount = getoutput(CMD % (inFile, )).strip().split('\t')
    return int(mapcount), int(unmapcount)

def stackedBarPlot( stackedBars, stackedLabels, barLabels, stackedColors=None ):
    """
    stackedBars = [ [29, 10, 21], [24, 11, 33] ]
    stackedLabels = ['stack1', 'stack2', 'stack3'] 
    barLabels = ['bar1', 'bar2']
    stackedColors = ['red', 'blue', 'green']
    stackedBarPlot( stackedBars, stackedLabels, barLabels, stackedColors)
    """
    
    import matplotlib.pyplot as plt
    
    def getStackedBarList(stackedBars, stackedLabels, barLabels):
        barList = []
        for bar,bar_label in zip(stackedBars,barLabels):
            for bar_item,stack_label in zip(bar,stackedLabels):
                barList.append( (bar_item, stack_label, bar_label) )
        return barList
    
    assert len(stackedBars) == len(barLabels)
    assert  len(stackedBars[0]) == len(stackedLabels)
    
    if stackedColors:
        assert len(stackedLabels) == len(stackedColors)
    else:
        stackedColors = sns.color_palette("hls", len(stackedLabels))
    
    stackedBarList = getStackedBarList(stackedBars, stackedLabels, barLabels)
    
    last_y = [0]*len(barLabels)
    for i,stack_label in enumerate( stackedLabels ):
        sub_dict = {it[2]:it[0] for it in stackedBarList if it[1]==stack_label}
        
        y = []
        for bar_label in barLabels:
            y.append(sub_dict[bar_label])
        
        plt.bar(range(1,len(barLabels)+1), y, color=stackedColors[i], bottom=last_y, label=stack_label, width=0.3, linewidth=0)
        last_y = [ y_i+y_j for y_i,y_j in zip(y, last_y) ]
    
    plt.legend(frameon=False)
    plt.xticks( range(1,len(barLabels)+1), barLabels )
    plt.xticks(rotation=45)


def count_pipe(pipefiles):
    readcount_list = []
    for i,file in enumerate(pipefiles):
        print("count: "+file)
        if file.endswith('fastq') or file.endswith('.fq') or file.endswith('fastq.gz') or file.endswith('fq.gz'):
            readcount = count_fq(file)
            readcount_list.append(readcount)
        elif file.endswith('.sam') or file.endswith('.bam'):
            mapcount,unmapcount = count_sam(file)
            if unmapcount == 0:
                readcount_list.append( mapcount )
            else:
                if i == len(pipefiles)-1:
                    readcount_list.append( (mapcount,unmapcount) )
                else:
                    readcount_list.append( mapcount+unmapcount )
    return readcount_list




params = init()
if not os.path.exists(params['outFolder']):
    os.mkdir(params['outFolder'])

if not os.path.exists(params['outFolder']+'img'):
    os.mkdir(params['outFolder']+'img')


p = Pool(len(params['filepipe']))
ReadCount = p.map(count_pipe, params['filepipe'])

has_final_sam = True
for data in ReadCount:
    if not isinstance(data[-1], tuple) and not isinstance(data[-1], list):
        has_final_sam = False

if has_final_sam:
    for data in ReadCount:
        if data[-2] != sum(data[-1]):
            sys.stderr.writelines("Warning: number of final sam file (%d) != number of read before map (%d)\n" % (data[-2], sum(data[-1])))

bar_data = []
for data in ReadCount:
    tmp_bar_data = []
    for i in range(1,len(data)-1):
        tmp_bar_data.append( data[i-1]-data[i] )
    if has_final_sam:
        tmp_bar_data.append( data[-1][0] )
        tmp_bar_data.append( data[-1][1] )
    else:
        tmp_bar_data.append( data[-1] )
    bar_data.append(tmp_bar_data)

processTags = params['processTags'][1:]
sampleTags = params['sampleTags']

if has_final_sam:
    processTags += ['un-mapped']

colors = ['#ffc107','#00ffff','#ffdddd','#8bc34a','#ffffcc','#cddc39','#2196f3','#87ceeb','#607d8b','#00bcd4']


width = int(len(processTags)/2 + 3)
plt.figure(figsize=(width, 5))
stackedBarPlot( bar_data, processTags, sampleTags, colors[:len(processTags)])
plt.xlabel("Sample")
plt.ylabel("Number of reads")
plt.tight_layout()
plt.savefig(params['outFolder']+'img/barplot.png')
plt.close()


table_head = "               <tr>\n                   <th>Sample name</th>\n"
for processname in [params['processTags'][0]]+processTags:
    table_head += "                   <th>{}</th>\n".format(processname)

table_head += "               </tr>\n"

table_content = ""
for samplename,data in zip(sampleTags,ReadCount):
    total = data[0]
    table_content += "               <tr>\n"
    table_content += "                   <td>"+samplename+"</td>\n"
    table_content += "                   <td>%s (%.2f%%)</td>\n" % ("{:,}".format(int(data[0])), 100.0)
    for value in data[1:-1]:
        table_content += "                   <td>%s (%.2f%%)</td>\n" % ("{:,}".format(int(value)), 100.0*value/total)
    if has_final_sam:
        table_content += "                   <td>%s (%.2f%%)</td>\n" % ("{:,}".format(int(data[-1][0])), 100.0*data[-1][0]/total)
        table_content += "                   <td>%s (%.2f%%)</td>\n" % ("{:,}".format(int(data[-1][1])), 100.0*data[-1][1]/total)
    else:
        table_content += "                   <td>%s (%.2f%%)</td>\n" % ("{:,}".format(int(data[-1])), 100.0*data[-1]/total)
    table_content += "               </tr>\n"

table = """
        <table class="w3-table-all" style="width:100%;overflow-x: scroll;">
{}
        </table>
""".format( table_head + table_content )

#print(table)

html = """
<!DOCTYPE html>
<html>
<head>
<title>Read distribution</title>
<meta http-equiv=Content-Type content="text/html;charset=utf-8">
<meta name="viewport" content="width=device-width;initial-scale=1.0">
<link rel="stylesheet" type="text/css" href="https://www.w3schools.com/w3css/4/w3.css">

</head>
<body>

<div id="id01" class="w3-modal" style="cursor:pointer">
  <div class="w3-modal-content w3-card-4" style="cursor:default">
    <img src="img/barplot.png" alt="Reads statistics" width="100%" />
  </div>
</div>

<div class="w3-content" style="max-width: 1000px">
    <div class="w3-panel w3-pink w3-center"><h3>Read distribution statistics</h3></div>
    <p><b>command</b></p>
    <p class="w3-border w3-lightgray" style="overflow:scroll;width:500px">
    {0}
    </p>

    <div class="w3-container w3-section">
        <h3>1. Reads distribution <b>(click to enlarge)</b> </h3>
        <p>The figure below shows the number of remaining reads for each processing step.</p>
        <img src="img/barplot.png" alt="Reads statistics" onclick="document.getElementById('id01').style.display='block'" />
    </div>

    <div class="w3-container w3-section">
        <h3>2. Reads left for each step </b></h3>
        {1}
    </div>
</div>
<script>
var modal = document.getElementById('id01');
window.onclick = function(event) {{
  if (event.target == modal) {{
    modal.style.display = "none";
  }}
}}
</script>
</body>
</html>
"""

command = " ".join(sys.argv)

OUT = open(params['outFolder']+'report.html', 'w')
OUT.writelines(html.format(command, table))
OUT.close()

