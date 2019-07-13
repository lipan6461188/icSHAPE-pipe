
import pandas as pd
import seaborn as sns
import os, sys
import matplotlib.pyplot as plt
import matplotlib
import getopt
import random
import version
font = {'size': 12}
matplotlib.rc('font', **font)
matplotlib.use('Agg')

Usage = """
plotGenomeSHAPEdist - base ratio and shape distribution
=======================================================
\x1b[1mUSAGE:\x1b[0m 
  %s -i input.gtab -o directory
\x1b[1mHELP:\x1b[0m
  -i                    <String>
                                Input genome shape file record the RT (prduced by calcSHAPE)
  -o                    <String>
                                Directory to output HTML report

  More options:
  --minBD               <Int>
                                Minimun basedensity (default: 100)

  Warning:
        @Base tag is required

\x1b[1mWARNING:\x1b[0m
    Basedensity is appended after the RT column

\x1b[1mVERSION:\x1b[0m
    %s

\x1b[1mAUTHOR:\x1b[0m
    Li Pan

""" % (sys.argv[0], version.Version)


def init():
    params = { 'inFile': None, 'outFolder': None, 'minBD': 100 }
    opts, args = getopt.getopt(sys.argv[1:], 'hi:o:', ['minBD='])
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
    elif 'Base' not in gTab_head:
        sys.stderr.writelines("Error: Base not in head\n")
        exit(-1)
    elif 'N_RT' not in gTab_head:
        sys.stderr.writelines("Error: N_RT not in head\n")
        exit(-1)
    elif 'N_BD' not in gTab_head:
        sys.stderr.writelines("Error: N_BD not in head\n")
        exit(-1)
    elif 'Shape' not in gTab_head:
        sys.stderr.writelines("Error: Shape not in head\n")
        exit(-1)
    
    return gTab_head

def count_Base_RTSHAPE(IN, gTab_head, minBD=100):
    
    base_col = gTab_head['Base']
    n_rt_col = gTab_head['N_RT']
    n_bd_col = gTab_head['N_BD']
    d_rt_col = gTab_head.get('D_RT', 0)
    shape_col = gTab_head['Shape']
    
    base_n_rt = { 'A':0, 'T':0, 'C':0, 'G':0 }
    base_d_rt = { 'A':0, 'T':0, 'C':0, 'G':0 }
    base_shape = { 'A':[], 'T':[], 'C':[], 'G':[] }
        
    for line in IN:
        data = line.strip().split()
        
        nrt = int(data[ n_rt_col-1 ])
        nbd = int(data[ n_bd_col-2 ])
        
        base = data[ base_col-1 ].upper()
        if base == 'U':
            base = 'T'
        if base == 'N':
            continue
        
        base_n_rt[base] += nrt
        if d_rt_col:
            drt = int(data[ d_rt_col-1 ])
            base_d_rt[base] += drt
        shape = data[ shape_col-1 ]
        if shape != '-1':
            base_shape[base].append( float(shape) )
    
    return base_n_rt, base_d_rt, base_shape

#####################################

def adjacent_values(vals, q1, q3):
    import numpy
    vals = sorted(vals)
    
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = numpy.clip(upper_adjacent_value, q3, vals[-1])
    
    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = numpy.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value

def set_axis_style(ax, labels):
    import numpy
    
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(numpy.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)

def violin(ax, data_list, labels, colors=None, rem_ext=0):
    """
    fig, axs = plt.subplots(nrows=5, ncols=1, figsize=(6, 12), sharey=True)
    
    axs[0].set_title('smartSHAPE 1 ng')
    axs[0].set_ylabel('smart SHAPE score')
    data = [ [],[],[],[] ]
    violin(axs[0], data, ['A', 'T', 'C', 'G'])
    
    fig.tight_layout()
    fig.show()
    """
    import numpy
    
    if colors == None:
        colors = ['#D43F3A'] * len(data_list)
    else:
        assert len(colors) == len(data_list)
        colors = colors[::-1]
    
    if rem_ext:
        assert 0.0 <= rem_ext <= 0.5
        import copy
        data_list = copy.deepcopy(data_list)
        for idx in range(len(data_list)):
            data_list[idx].sort()
            remNum = int(len(data_list[idx]) * rem_ext)
            start = remNum; end = len(data_list[idx]) - remNum
            data_list[idx] = data_list[idx][start:end]
    
    parts = ax.violinplot(data_list, showmeans=False, showmedians=False, showextrema=False)
    
    for pc in parts['bodies']:
        pc.set_facecolor(colors.pop())
        pc.set_edgecolor('black')
        pc.set_alpha(1)
    
    quartile1 = []; medians = []; quartile3 = []
    for data in data_list:
        quartile1.append( numpy.percentile(data, 25) )
        medians.append( numpy.percentile(data, 50) )
        quartile3.append( numpy.percentile(data, 75) )
    
    whiskers = numpy.array([
        adjacent_values(sorted_array, q1, q3)
        for sorted_array, q1, q3 in zip(data_list, quartile1, quartile3)])
    
    whiskersMin, whiskersMax = whiskers[:, 0], whiskers[:, 1]
    
    inds = numpy.arange(1, len(medians) + 1)
    ax.scatter(inds, medians, marker='o', color='white', s=20, zorder=3)
    ax.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=8)
    ax.vlines(inds, whiskersMin, whiskersMax, color='k', linestyle='-', lw=1)
    
    set_axis_style(ax, labels)

#####################################

def plot_base_pie(base_ratio_dict, outfn):
    colors = ['gold', 'yellowgreen', 'lightcoral', 'lightskyblue']
    baselist = sorted(list(base_ratio_dict.keys()))
    base_num = [ base_ratio_dict[b] for b in baselist ]
    total = sum(base_num)
    base_ratio = [ 1.0*it/total for it in base_num ]
    explode = (0, 0, 0, 0)
    plt.figure(figsize=(4,4))
    baseratiolabels = [ "%s(%d%%)"%(b,100*r) for b,r in zip(baselist,base_ratio) ]
    plt.pie(base_ratio, labels=baseratiolabels, colors=colors, startangle=140, shadow=False, explode=explode, rotatelabels=False)
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig(outfn)
    plt.close()

def plot_base_shape_violin(base_shape, outfn):
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6, 6), sharey=True)
    colors = ['gold', 'yellowgreen', 'lightcoral', 'lightskyblue']
    baselist = sorted(list(base_shape.keys()))
    
    axs.set_ylabel('SHAPE score')
    data = [ base_shape[b] for b in baselist ]
    i = 0
    sample_data = []
    for i in range(len(data)):
        if len(data[i]) > 100000:
            sample_data.append( random.sample(data[i], 100000) )
        elif len(data[i]) == 0:
            sample_data.append( [0.3]*100 )
        else:
            sample_data.append( data[i] )
    violin(axs, sample_data, baselist, colors=colors)
    
    fig.tight_layout()
    fig.savefig(outfn)

params = init()
if not os.path.exists(params['outFolder']):
    os.mkdir(params['outFolder'])

if not os.path.exists(params['outFolder']+'img'):
    os.mkdir(params['outFolder']+'img')

## Read file and calculate
IN = open(params['inFile'], 'r')
gTab_head = read_gTab_head(IN)

base_n_rt, base_d_rt, base_shape = count_Base_RTSHAPE(IN, gTab_head, minBD=params['minBD'])
IN.close()

plot_base_pie(base_n_rt, params['outFolder']+'img/nrt.png')
if sum(list(base_d_rt.values())):
    plot_base_pie(base_d_rt, params['outFolder']+'img/drt.png')

plot_base_shape_violin(base_shape, params['outFolder']+'img/shapevioline.png')

html_head = """
<!DOCTYPE html>
<html>
<head>
<title>{0} reactivity distribution</title>
<meta http-equiv=Content-Type content="text/html;charset=utf-8">
<meta name="viewport" content="width=device-width;initial-scale=1.0">
<link rel="stylesheet" type="text/css" href="https://www.w3schools.com/w3css/4/w3.css">
</head>
<body>
<div class="w3-content" style="max-width: 1000px">
    <div class="w3-panel w3-pink w3-center"><h3>{0}</h3></div>
    <p>command: {1}</p>
    <div class="w3-container w3-section">
        <h3>1. Ratio of RT for each base</h3>
        <p>The figure below shows the RT percentage of each base in each sample.</p>
        <div style="float:left;max-width:300px;padding:5px;">
            <div class="w3-center" style="width:100%;height:20%">
                <span>RT Ratio in treatment samples</span>
            </div>
            <div style="width:100%;height:80%">
                <img src="img/nrt.png" alt="Treatment Base RT Ratio" width="100%"/>
            </div>
        </div>
"""

html_tail = """
    </div>
    <div class="w3-container w3-section">
        <h3>2. Reactivity distribution</h3>
        <p>The figure below shows the distribution of SHAPE scores for each base.</p>
        <img src="img/shapevioline.png" alt="Reactivity distribution" style="max-width:500px;" />
    </div>
</div>
</body>
</html>
"""

html_medium = """
        <div style="float:left;max-width:300px;padding:5px;">
            <div class="w3-center" style="width:100%;height:20%">
                <span>RT Ratio in control samples</span>
            </div>
            <div style="width:100%;height:80%">
                <img src="img/drt.png" alt="Control Base RT Ratio" width="100%"/>
            </div>
        </div>
"""

pure_file = params['inFile'].split('/')[-1]
command = " ".join(sys.argv)

OUT = open(params['outFolder']+'report.html', 'w')
OUT.writelines(html_head.format(pure_file, command))

if base_d_rt['A']:
    OUT.writelines(html_medium)

OUT.writelines(html_tail)
OUT.close()

