
import pandas as pd
import seaborn as sns
import os, sys
import matplotlib.pyplot as plt
import matplotlib
import getopt
import version
font = {'size': 12}
matplotlib.rc('font', **font)
matplotlib.use('Agg')

Usage = """
plotGenomeRTRepCor - replicate correlation, base ratio and number of RT
=======================================================================
\x1b[1mUSAGE:\x1b[0m 
  %s -i input.gtab -o directory
\x1b[1mHELP:\x1b[0m
  -i                    <String>
                                Input countRT file record the RT (prduced by countRT)
  -o                    <String>
                                Directory to output HTML report

  More options:
  --minBD               <Int>
                                Minimun basedensity (default: 100)
  --winSize             <Int>
                                Window size for each replicate calculation (default: 100)
  --bases               <char,char,...>
                                Bases to calculate Correlation (default: All)

\x1b[1mWARNING:\x1b[0m
    Basedensity is appended after the RT column

\x1b[1mVERSION:\x1b[0m
    %s

\x1b[1mAUTHOR:\x1b[0m
    Li Pan

""" % (sys.argv[0], version.Version)


def init():
    params = { 'inFile': None, 'outFolder': None, 'minBD': 100, 'winSize': 100, 'bases': None }
    opts, args = getopt.getopt(sys.argv[1:], 'hi:o:', ['minBD=', 'winSize=', 'bases='])
    for op, value in opts:
        if op == '-h':
            print(Usage)
            sys.exit(-1)
        # Basic Parameters
        elif op == '-i':
            params['inFile'] = os.path.abspath(value)
        elif op == '-o':
            params['outFolder'] = os.path.abspath(value).rstrip('/') + '/'
        elif op == '--minBD':
            params['minBD'] = int(value)
        elif op == '--winSize':
            params['winSize'] = int(value)
        elif op == '--bases':
            params['bases'] = value.strip().split(',')
            for b in params['bases']:
                if len(b) != 1:
                    sys.stderr.writelines("Error: --bases should be bases seperated with comma. such as --bases A,T,C,G\n")
                    exit(-1)
        else:
            sys.stderr.writelines("Error: unrecognized parameter: "+op+"\n")
            print(Usage)
            sys.exit(-1)
    
    # check
    if (not params['inFile']) or (not params['outFolder']):
          sys.stderr.writelines("Error: Please specify -i and -o"+"\n")
          print(Usage)
          sys.exit(-1)
    
    return params


def read_gTab_head(IN):
    last_pos = IN.tell()
    
    gTab_head = {}
    
    rt_files = []
    bd_files = []
    
    line = IN.readline()
    while line:
        if not line.startswith('@'):
            ### Check column nums
            data = line.strip().split()
            if len(data) != gTab_head['ColNum']:
                sys.stderr.writelines("Error: actual column number (%s) != labeled number (%s)" % (len(data), gTab_head['ColNum'])+"\n")
                exit(-1)
            
            IN.seek(last_pos)
            break
        tag, num = line.strip()[1:].split()
        
        if tag.startswith('RT_'):
            rt_files.append( tag[3:] )
        elif tag.startswith('BD_'):
            bd_files.append( tag[3:] )
        
        gTab_head[tag] = int(num)
        last_pos = IN.tell()
        line = IN.readline()    
    
    if 'ChrID' not in gTab_head:
        sys.stderr.writelines("Error: ChrID not in head\n")
        exit(-1)
    elif 'ColNum' not in gTab_head:
        sys.stderr.writelines("Error: ColNum not in head\n")
        exit(-1)
    elif 'Strand' not in gTab_head:
        sys.stderr.writelines("Error: Strand not in head\n")
        exit(-1)
    elif 'ChrPos' not in gTab_head:
        sys.stderr.writelines("Error: ChrPos not in head\n")
        exit(-1)
    
    #rt_files.sort()
    #bd_files.sort()
    
    if sorted(rt_files) != sorted(bd_files):
        print("Error: rt files is different with bd files")
        print("rt_files: "+" ".join(rt_files))
        print("bd_files: "+" ".join(bd_files))
        exit(-1)
    
    if len(rt_files) == 0:
        sys.stderr.writelines("Error: it seems not a .gtab file produced by countRT")
        exit(-1)
    
    rt_bd_cols = { file:(gTab_head['RT_'+file], gTab_head['BD_'+file]) for file in rt_files }
    
    return gTab_head,rt_bd_cols,rt_files

"""
def parse_head_as_rtcount(gTab_head):
    rt_files = []
    bd_files = []
    for key in gTab_head:
        if key.startswith('RT_'):
            rt_files.append( key.lstrip('RT_') )
        elif key.startswith('BD_'):
            bd_files.append( key.lstrip('BD_') )
    rt_files.sort()
    bd_files.sort()
    if rt_files != bd_files:
        print("Error: rt files is different with bd files")
        exit(-1)
    if len(rt_files) == 0:
        sys.stderr.writelines("Error: it seems not a .gtab file produced by countRT")
        exit(-1)
    rt_bd_cols = { file:(gTab_head['RT_'+file], gTab_head['BD_'+file]) for file in rt_files }
    return rt_bd_cols
"""

def calcRTReplicateCorrelation(IN, gTab_head, rt_bd_cols, files, bases=None, minBD=100, windowsize=100):
    import scipy
    import scipy.stats
    
    base_col = gTab_head.get('Base', 0)
    
    file_base = {}
    
    for file in files:
        file_base[file] = { 'A':0, 'T':0, 'C':0, 'G':0 }
    
    file_cor = {}
    rt_list = {}
    for i in range(len(files)):
        for j in range(i+1, len(files)):
            key = files[i]+files[j]
            rt_list[key] = []
            file_cor[key] = []    
    
    for line in IN:
        data = line.strip().split()
        tmp_rt_list = []
        tmp_bd_list = []
        for file in files:
            rt = int(data[ rt_bd_cols[file][0]-1 ])
            bd = int(data[ rt_bd_cols[file][1]-1 ])
            if base_col:
                base = data[ base_col-1 ].upper()
                if base == 'U':
                    base = 'T'
                if base != 'N':
                    file_base[file][base] += rt
                tmp_rt_list.append( rt )
                tmp_bd_list.append( bd )
        
        if bases and base not in bases:
            continue
        
        for i in range(len(files)):
            fi = files[i]
            for j in range(i+1, len(files)):
                fj = files[j]
                if tmp_bd_list[i]>=minBD and tmp_bd_list[j]>=minBD:
                    key = fi+fj
                    rt_list[key].append( [tmp_rt_list[i], tmp_rt_list[j]] )
                    if len(rt_list[key]) == windowsize:
                        x = [ it[0] for it in rt_list[key] ]
                        y = [ it[1] for it in rt_list[key] ]
                        r, p = scipy.stats.pearsonr(x, y)
                        if p < 0.05:
                            file_cor[key].append( r )
                        rt_list[key] = []
    
    return file_cor, file_base, files

def count_base_ratio(file_base):
    file_base_ratio = {}
    baseset = set()
    for file in file_base:
        total = sum( list(file_base[file].values()) )
        if total == 0: return False, False, False
        file_base_ratio[file] = {}
        for base in file_base[file]:
            file_base_ratio[file][base] = 1.0*file_base[file][base] / (total+1)
            baseset.add(base)
    
    # check bases
    baselist = list(baseset)
    base_ratio_list = []
    files = list(file_base_ratio.keys())
    for file in files:
        rl = [ file_base_ratio[file][base] for base in baselist ]
        base_ratio_list.append(rl)
    return base_ratio_list, files, baselist

def plot_base_pie(base_ratio_list, files, baselist, outprefix):
    colors = ['gold', 'yellowgreen', 'lightcoral', 'lightskyblue'];
    files = pure_file_name(files)
    for i in range(len(files)):
        base_ratio = base_ratio_list[i]
        explode = (0, 0, 0, 0);
        plt.figure(figsize=(4,4))
        baseratiolabels = [ "%s(%d%%)"%(b,100*r) for b,r in zip(baselist,base_ratio) ]
        plt.pie(base_ratio, labels=baseratiolabels, colors=colors, startangle=140, shadow=False, explode=explode, rotatelabels=False);
        plt.axis('equal');
        plt.tight_layout()
        plt.savefig(outprefix+files[i]+".png")
        plt.close()

def init_list_rect(rowNum, colNum, init_value=0):
    rect = []
    
    for i in range(rowNum):
        row = []
        for j in range(colNum):
            row.append(init_value)
        rect.append(row)
    
    return rect

def calc_averaged_correlation(file_cor, files):
    prect = init_list_rect(len(files), len(files), init_value=0)
    i = 0
    vmin=1
    vmax=0
    for i in range(len(files)):
        fi = files[i]
        for j in range(i+1, len(files)):
            fj = files[j]
            key = fi + fj
            if i == j:
                prect[i][j] = prect[j][i] = 1
            else:
                #print(file_cor.keys())
                v = sum(file_cor[key])/len(file_cor[key])
                prect[i][j] = prect[j][i] = v
                if v > vmax: vmax = v
                if v < vmin: vmin = v
            j += 1
        i += 1
    return prect, vmin, vmax

def pure_file_name(raw_file_list):
    pure_file_list = []
    for file in raw_file_list:
        if '.' in file:
            pure_file_list.append( ".".join(file.split('.')[:-1]) )
        else:
            pure_file_list.append( file )
    return pure_file_list

def hex2rgb(hex):
    return tuple(int(hex[i:i+2], 16) for i in (0, 2, 4))

def plot_seqdepth(file_base, files, outfn):
    purefiles = pure_file_name(files)
    depth = []
    for file,purefile in zip(files,purefiles):
        depth.append( [sum(file_base[file].values()), purefile] )
    depth = pd.DataFrame(depth, columns=['Number of rt','file'])
    width = 1 + len(files)
    plt.figure(figsize=(width,5))
    colors = ['#f44336', '#2196f3', '#3f51b5', '#9c27b0', '#009688', '#673ab7', '#ff5722', '#e91e63', '#ff9800'] * 2
    #colors = [ hex2rgb(hex[1:]) for hex in colors ]
    sns.barplot(data=depth, x='file', y='Number of rt', palette=colors[:len(depth)])
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(outfn)
    plt.close()

def plot_cor_heatmap(file_cor, files, outfn):
    if len(files) == 1:
        return False
    prect, vmin, vmax = calc_averaged_correlation(file_cor, files)
    prect = pd.DataFrame(prect, columns=pure_file_name(files), index=pure_file_name(files))
    mask = prect.copy()
    mask.loc[:,:] = False
    for i in range(mask.shape[0]):
        mask.iloc[i,i] = True
    
    size = 4 + int(len(files)/2)
    plt.figure(figsize=(size,size))
    sns.heatmap(data=prect, cmap=sns.color_palette("BuGn"), vmin=vmin, vmax=vmax, mask=mask, annot=True, fmt='.3f', cbar=False)
    plt.tight_layout()
    plt.savefig(outfn)
    plt.close()
    return True


params = init()
if not os.path.exists(params['outFolder']):
    os.mkdir(params['outFolder'])

if not os.path.exists(params['outFolder']+'img'):
    os.mkdir(params['outFolder']+'img')

## Read file and calculate
IN = open(params['inFile'], 'r')
gTab_head,rt_bd_cols,files = read_gTab_head(IN)

has_genome = ('Base' in gTab_head)
if params['bases'] and (not has_genome):
    sys.stderr.writelines("Error: --genome should be provided when --bases specified\n")
    exit(-1)

file_cor, file_base, files = calcRTReplicateCorrelation(IN, gTab_head, rt_bd_cols, files, bases=params['bases'], minBD=params['minBD'], windowsize=params['winSize'])
IN.close()

plot_seqdepth(file_base, files, params['outFolder']+'img/seqdepth.png')

if has_genome:
    base_ratio_list, files2, baselist = count_base_ratio(file_base)
    if base_ratio_list != False:
        plot_base_pie(base_ratio_list, files2, baselist, params['outFolder']+'img/')

heatmap_ploted = plot_cor_heatmap(file_cor, files, params['outFolder']+'img/cor_heatmap.png')

html_head = """
<!DOCTYPE html>
<html>
<head>
<title>{0} replicate report</title>
<meta http-equiv=Content-Type content="text/html;charset=utf-8">
<meta name="viewport" content="width=device-width;initial-scale=1.0">
<link rel="stylesheet" type="text/css" href="https://www.w3schools.com/w3css/4/w3.css">
<style>
.hover_full {{
    width: 500px;
}}
.hover_full:hover{{
    width:1000px;
}}
</style>
</head>
<body>
<div class="w3-content" style="max-width: 1000px">
    <div class="w3-panel w3-purple w3-center"><h3>{0}</h3></div>
    <p>command: {1}</p>
"""

html_tail = """
</div>
</body>
</html>
"""

html_body_1 = """
    <div class="w3-container w3-section">
        <h3>1. Number of RT <b>(Hover to enlarge the figure)</b></h3>
        <p>The figure below counts the total number of RTs in each sample.</p>
        <img src="img/seqdepth.png" alt="Number of RT" class="hover_full" />
    </div>
"""

html_body_2 = """
    <div class="w3-container w3-section">
        <h3>2. Ratio of RT for each base</h3>
        <p>The figure below shows the RT percentage of each base in each sample.</p>
{}
    </div>
"""


html_body_3 = """
    <div class="w3-container w3-section">
        <h3>3. Correlation heatmap {} <b>(Hover to enlarge the figure)</b></h3>
        <p>The figure below shows the correlation between the samples. Filter out the bases on the genome with coverage greater than the threshold, slide with a sliding window and calculate the Pearson correlation coefficient for each window</p>
        <img src="img/cor_heatmap.png" alt="Correlation heatmap" class="hover_full" />
    </div>
"""

html_body_2_template = """
        <div style="float:left;max-width:300px;padding:5px;">
            <div class="w3-center" style="width:100%;height:20%">
                <span>{0}</span>
            </div>
            <div style="width:100%;height:80%">
                <img src="img/{0}.png" alt="{0}" width="100%"/>
            </div>
        </div>
"""

pure_file = params['inFile'].split('/')[-1]
command = " ".join(sys.argv)

OUT = open(params['outFolder']+'report.html', 'w')
OUT.writelines(html_head.format(pure_file, command))
OUT.writelines(html_body_1)

if has_genome:
    html_2_body = ""
    for file in pure_file_name(files):
        html_2_body += html_body_2_template.format(file)
    OUT.writelines(html_body_2 .format(html_2_body))
else:
    info = "        <p style=\"color:#f44336\">Base information is not included in the input file, specify -genome in countRT</p>"
    OUT.writelines(html_body_2.format(info))

if heatmap_ploted:
    if params['bases']:
        info = "(only {0} bases)".format(",".join(params['bases']))
    else:
        info = ""
    OUT.writelines(html_body_3.format(info))


OUT.writelines(html_tail)
OUT.close()

