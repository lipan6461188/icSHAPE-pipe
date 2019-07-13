#!/usr/bin/env bash
set -e
#trap 'exit 130' INT
#trap "exit 1" TERM

BOLD_s="\x1b[1;39;49m"
RED_s="\x1b[0;31;107m"
END="\x1b[0m"

function abspath(){ echo ${1%/*}; }
function absfn(){ echo ${1##*/}; }
function fnroot() { echo ${1%.*}; }
function fnrep() { echo $(fnroot $(absfn $1)); }

export -f abspath
export -f absfn
export -f fnroot
export -f fnrep

function print_usage()
{
    echo -e "${BOLD_s}Usage:${END}"
    echo -e "    $proc_name [-D D1.fastq,D2.fastq...] -N N1.fastq,N2.fastq... -f genome.fa [-i STAR_genome_index]"
    echo -e "               -g genome.gtf -s ncbi/ensembl -r rRNA.fa -a adaptor.fa -o outdir -p 20 -l 13 [-q Z-ZQF] [-v] [-c]"
    echo -e "${BOLD_s}Options:${END}"
    echo -e "    -D     -- DMSO samples (optional)"
    echo -e "    -N     -- NAI samples"
    echo -e "    -f     -- reference genome (.fasta file)"
    echo -e "    -i     -- STAR index of reference genome (optional)"
    echo -e "    -g     -- GTF/GFF3 annotation"
    echo -e "    -s     -- annotation source: ncbi/ensembl"
    echo -e "    -r     -- abundant sequence (such as rRNA,tRNA...)"
    echo -e "    -a     -- adaptor sequence"
    echo -e "    -o     -- output folder"
    echo -e "    -p     -- how many threads to use at most"
    echo -e "    -l     -- how many bases to trim at 5'"
    echo -e "    -q     -- LSF queue (default: Z-ZQF)"
    echo -e "    -v     -- verbose mode, output command but out execute"
    echo -e "    -c     -- The input reads are cleaned reads. readCollapse,trim will be skipped"
    echo -e ""
    echo -e "${BOLD_s}Version:${END}"
    echo -e "    2019-06-03"
    echo -e ""
    echo -e "${BOLD_s}Author:${END}"
    echo -e "    Li Pan (lip16@mails.tsinghua.edu.cn)"
    echo -e "\n"
}

function exeInPath()
{
    exe=$1
    if [ "$(which $exe)" == "" ]; then
        echo "false"
    else
        echo "true"
    fi
}

function job_id() { echo $2 | gawk 'match($0, /<([0-9]+)>/, id){print id[1]}'; }

export -f exeInPath
export -f job_id

declare proc_name=$0;
declare command="${@: -1}";

declare BinDir=$(abspath ${proc_name})

declare DMSOFiles=()
declare NAIFiles=()
declare GenomeFasta=""
declare GenomeIndex=""
declare GenomeGTF=""
declare GTFSource=""
declare rRNAFasta=""
declare Adaptor=""
declare OutDir=""
declare Cores=""
declare LeadTrim=""
declare verbose=0
declare isCleaned="false"

declare LSFqueue="Z-ZQF"
declare Threads=20


while getopts ":D:N:f:i:g:s:r:a:o:p:l:q:cv" opt; do 
  case ${opt} in
    D )
        IFS=',' read -r -a DMSOFiles <<< "${OPTARG}";
      ;;
    N )
        IFS=',' read -r -a NAIFiles <<< "${OPTARG}";
      ;;
    f )
        GenomeFasta=${OPTARG};
      ;;
    i )
        GenomeIndex=${OPTARG};
      ;;
    g )
        GenomeGTF=${OPTARG};
      ;;
    s )
        GTFSource=${OPTARG};
      ;;
    r )
        rRNAFasta=${OPTARG};
      ;;
    a )
        Adaptor=${OPTARG};
      ;;
    o )
        OutDir=$(echo ${OPTARG} | sed 's/\///');
      ;;
    p )
        Cores=${OPTARG};
      ;;
    l )
        LeadTrim=${OPTARG};
      ;;
    q )
        LSFqueue=${OPTARG};
      ;;
    v )
        verbose=1;
      ;;
    c )
        isCleaned="true";
      ;;
    \? )
        print_usage;
      ;;
  esac
done
shift $((OPTIND -1))

if [[ ${#NAIFiles[@]} -eq 0 ]]; then
    echo -e "${RED_s}Error: Specify -N${END}";
    print_usage;
    exit;
elif [[ "${GenomeFasta}" == "" ]]; then
    echo -e "${RED_s}Error: Specify -f${END}";
    print_usage;
    exit;
elif [[ "${GenomeGTF}" == "" ]]; then
    echo -e "${RED_s}Error: Specify -g${END}";
    print_usage;
    exit;
elif [[ "${GTFSource}" == "" ]]; then
    echo -e "${RED_s}Error: Specify -s${END}";
    print_usage;
    exit;
elif [[ "${rRNAFasta}" == "" ]]; then
    echo -e "${RED_s}Error: Specify -r${END}";
    print_usage;
    exit;
elif [[ "${OutDir}" == "" ]]; then
    echo -e "${RED_s}Error: Specify -o${END}";
    print_usage;
    exit;
elif [[ "${Cores}" == "" ]]; then
    echo -e "${RED_s}Error: Specify -p${END}";
    print_usage;
    exit;
elif [[ "${isCleaned}" == "false" ]] && [[ "${Adaptor}" == "" ]]; then
    echo -e "${RED_s}Error: Specify -a${END}";
    print_usage;
    exit;
elif [[ "${isCleaned}" == "false" ]] && [[ "${LeadTrim}" == "" ]]; then
    echo -e "${RED_s}Error: Specify -l${END}";
    print_usage;
    exit;
fi

declare errorFn=$OutDir/error
declare logFn=$OutDir/log

####################
## Detection
####################

if [ ! -f ${BinDir}/fetchSmallRNA.py ]; then
    echo -e "${RED_s}Error: ${BinDir}/fetchSmallRNA.py not exist${END}"
    exit;
fi

if [ ! -f ${BinDir}/falen ]; then
    echo -e "${RED_s}Error: ${BinDir}/falen not exist${END}"
    exit;
fi

if [ "$(exeInPath bwait)" == "false" ]; then
    echo -e ${RED_s}"Error: bwait not in your PATH"${END};
    exit;
fi

if [ "$(exeInPath bsub)" == "false" ]; then
    echo -e ${RED_s}"Error: bsub not in your PATH"${END};
    exit;
fi

for abssample in ${DMSOFiles[@]};
do
    if [ ! -f ${abssample} ]; then
        echo -e ${RED_s}"Error: ${abssample} not exist"${END};
        exit;
    fi
done

for abssample in ${NAIFiles[@]};
do
    if [ ! -f ${abssample} ]; then
        echo -e ${RED_s}"Error: ${abssample} not exist"${END};
        exit;
    fi
done

if [ ! -f ${GenomeFasta} ]; then
    echo -e ${RED_s}"Error: ${GenomeFasta} not exist"${END};
    exit;
fi  

if [ ! -f ${GenomeGTF} ]; then
    echo -e ${RED_s}"Error: ${GenomeGTF} not exist"${END};
    exit;
fi 

if [ ! -f ${rRNAFasta} ]; then
    echo -e ${RED_s}"Error: ${rRNAFasta} not exist"${END};
    exit;
fi 

if [ "${isCleaned}" == "false" ] && [ ! -f ${Adaptor} ]; then
    echo -e ${RED_s}"Error: ${Adaptor} not exist"${END};
    exit;
fi

####################
## Create directory
####################

CMD1=""
if [ "${GenomeIndex}" == "" ]; then
    CMD1="mkdir -p ${OutDir}/index"; echo -e "  >>  ${BOLD_s}${CMD1}${END}";
fi

CMD2="mkdir -p ${OutDir}/rRNA_index"; echo -e "  >>  ${BOLD_s}${CMD2}${END}";
CMD3="mkdir -p ${OutDir}/smallRNA"; echo -e "  >>  ${BOLD_s}${CMD3}${END}";
CMD4="mkdir -p ${OutDir}/GTF"; echo -e "  >>  ${BOLD_s}${CMD4}${END}";
CMD5="mkdir -p ${OutDir}/1.readCollapse"; echo -e "  >>  ${BOLD_s}${CMD5}${END}";
CMD6="mkdir -p ${OutDir}/2.trim"; echo -e "  >>  ${BOLD_s}${CMD6}${END}";
CMD7="mkdir -p ${OutDir}/3.rem_rRNA"; echo -e "  >>  ${BOLD_s}${CMD7}${END}";
CMD8="mkdir -p ${OutDir}/4.mapSmallRNA"; echo -e "  >>  ${BOLD_s}${CMD8}${END}";
CMD9="mkdir -p ${OutDir}/5.mapGenome"; echo -e "  >>  ${BOLD_s}${CMD9}${END}";
CMD10="mkdir -p ${OutDir}/6.calcFPKM"; echo -e "  >>  ${BOLD_s}${CMD10}${END}";
CMD11="mkdir -p ${OutDir}/7.sam2tab"; echo -e "  >>  ${BOLD_s}${CMD11}${END}";
CMD12="mkdir -p ${OutDir}/8.calcGenomeSHAPE"; echo -e "  >>  ${BOLD_s}${CMD12}${END}";
CMD13="mkdir -p ${OutDir}/9.transSHAPE"; echo -e "  >>  ${BOLD_s}${CMD13}${END}";
CMD14="mkdir -p ${OutDir}/10.bedGraph"; echo -e "  >>  ${BOLD_s}${CMD14}${END}";
CMD15="mkdir -p ${OutDir}/11.quality_control"; echo -e "  >>  ${BOLD_s}${CMD15}${END}";

if [ $verbose -eq 0 ]; then
    if [ "${CMD1}" != "" ]; then
        eval ${CMD1}
    fi
    eval ${CMD2}
    eval ${CMD3}
    eval ${CMD4}
    eval ${CMD5}
    eval ${CMD6}
    eval ${CMD7}
    eval ${CMD8}
    eval ${CMD9}
    eval ${CMD10}
    eval ${CMD11}
    eval ${CMD12}
    eval ${CMD13}
    eval ${CMD14}
    eval ${CMD15}
fi


####################
## Preparation 
####################

## Build index for rRNA
CMD="bowtie2-build --quiet ${rRNAFasta} ${OutDir}/rRNA_index/rRNA"
echo -e "  >>  ${BOLD_s}${CMD}${END}";
if [ $verbose -eq 0 ]; then eval ${CMD}; fi

CMD="${BinDir}/falen ${rRNAFasta} > ${OutDir}/rRNA_index/rRNA.len"
echo -e "  >>  ${BOLD_s}${CMD}${END}";
if [ $verbose -eq 0 ]; then eval ${CMD}; fi
#${BinDir}/falen ${rRNAFasta} > ${OutDir}/rRNA_index/rRNA.len

## Build index for genome

jobIDs=()
declare i=0
if [ "${GenomeIndex}" == "" ]; then
    GenomeIndex=${OutDir}/index
    CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} -n ${Cores} \"icSHAPE-pipe starbuild -i ${GenomeFasta} -o ${GenomeIndex}/ --gtf ${GenomeGTF} -p ${Cores}\""
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then 
        info=$(eval ${CMD})
        jobIDs+=("$(job_id $info)")
        i=$[$i+1]
    fi
fi

for jobid in ${jobIDs[@]}; 
do
    CMD="bwait $jobid"
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then
        echo -e "starbuild bwait $jobid";
        eval ${CMD}
    fi
done

## Covert GTF file

jobIDs=()
declare i=0

CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} \"icSHAPE-pipe parseGTF -g ${GenomeGTF} -o ${OutDir}/GTF/Anno -s ${GTFSource} --genome ${GenomeFasta}\""
echo -e "  >>  ${BOLD_s}${CMD}${END}";
if [ $verbose -eq 0 ]; then 
    info=$(eval ${CMD})
    jobIDs+=("$(job_id $info)")
    i=$[$i+1]
fi

for jobid in ${jobIDs[@]}; 
do
    CMD="bwait $jobid"
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then
        echo -e "parseGTF bwait $jobid";
        eval ${CMD}
    fi
done

## prepare smallRNA
CMD="${BinDir}/fetchSmallRNA.py ${OutDir}/GTF/Anno_transcriptome.fa ${OutDir}/smallRNA/smallRNA.fa 200"
echo -e "  >>  ${BOLD_s}${CMD}${END}";
if [ $verbose -eq 0 ]; then eval ${CMD}; fi

CMD="bowtie2-build --quiet ${OutDir}/smallRNA/smallRNA.fa ${OutDir}/smallRNA/smallRNA"
echo -e "  >>  ${BOLD_s}${CMD}${END}";
if [ $verbose -eq 0 ]; then eval ${CMD}; fi

CMD="${BinDir}/falen ${OutDir}/smallRNA/smallRNA.fa > ${OutDir}/smallRNA/smallRNA.len"
echo -e "  >>  ${BOLD_s}${CMD}${END}";
if [ $verbose -eq 0 ]; then eval ${CMD}; fi

####################
## Data process
####################

######## 1. readCollapse

if [ "${isCleaned}" == "false" ]; then
    jobIDs=()
    declare i=0

    for abssample in ${DMSOFiles[@]};
    do
        oFn=${OutDir}/1.readCollapse/$(fnrep ${abssample}).fastq
        CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} \"icSHAPE-pipe readcollapse -U ${abssample} -o ${oFn} --simplify\""
        echo -e "  >>  ${BOLD_s}${CMD}${END}";
        if [ $verbose -eq 0 ]; then 
            info=$(eval ${CMD})
            jobIDs+=("$(job_id $info)")
            i=$[$i+1]
        fi
    done

    for abssample in ${NAIFiles[@]};
    do
        oFn=${OutDir}/1.readCollapse/$(fnrep ${abssample}).fastq
        CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} \"icSHAPE-pipe readcollapse -U ${abssample} -o ${oFn} --simplify\""
        echo -e "  >>  ${BOLD_s}${CMD}${END}";
        if [ $verbose -eq 0 ]; then 
            info=$(eval ${CMD})
            jobIDs+=("$(job_id $info)")
            i=$[$i+1]
        fi
    done

    for jobid in ${jobIDs[@]}; 
    do
        CMD="bwait $jobid"
        echo -e "  >>  ${BOLD_s}${CMD}${END}";
        if [ $verbose -eq 0 ]; then
            echo -e "readcollapse bwait $jobid";
            eval ${CMD}
        fi
    done
fi

######## 2. Trim reads

if [ "${isCleaned}" == "false" ]; then
    jobIDs=()
    declare i=0

    for abssample in ${DMSOFiles[@]};
    do
        iFn=${OutDir}/1.readCollapse/$(fnrep ${abssample}).fastq
        oFn=${OutDir}/2.trim/$(fnrep ${abssample}).fastq
        CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} -n ${Cores} \"icSHAPE-pipe trim -i ${iFn} -o ${oFn} -l ${LeadTrim} -a ${Adaptor} -p ${Cores} -m 25\""
        echo -e "  >>  ${BOLD_s}${CMD}${END}";
        if [ $verbose -eq 0 ]; then 
            info=$(eval ${CMD})
            jobIDs+=("$(job_id $info)")
            i=$[$i+1]
        fi
    done

    for abssample in ${NAIFiles[@]};
    do
        iFn=${OutDir}/1.readCollapse/$(fnrep ${abssample}).fastq
        oFn=${OutDir}/2.trim/$(fnrep ${abssample}).fastq
        CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} -n ${Cores} \"icSHAPE-pipe trim -i ${iFn} -o ${oFn} -l ${LeadTrim} -a ${Adaptor} -p ${Cores} -m 25\""
        echo -e "  >>  ${BOLD_s}${CMD}${END}";
        if [ $verbose -eq 0 ]; then 
            info=$(eval ${CMD})
            jobIDs+=("$(job_id $info)")
            i=$[$i+1]
        fi
    done

    for jobid in ${jobIDs[@]}; 
    do
        CMD="bwait $jobid"
        echo -e "  >>  ${BOLD_s}${CMD}${END}";
        if [ $verbose -eq 0 ]; then
            echo -e "trim bwait $jobid";
            eval ${CMD}
        fi
    done
fi

######## 3.0 link

if [ "${isCleaned}" == "true" ]; then
    for abssample in ${DMSOFiles[@]};
    do
        oFn=${OutDir}/2.trim/$(fnrep ${abssample}).fastq
        CMD="ln -s $(realpath ${abssample}) ${oFn}"
        echo -e "  >>  ${BOLD_s}${CMD}${END}";
        if [ $verbose -eq 0 ]; then eval ${CMD}; fi
    done

    for abssample in ${NAIFiles[@]};
    do
        oFn=${OutDir}/2.trim/$(fnrep ${abssample}).fastq
        CMD="ln -s $(realpath ${abssample}) ${oFn}"
        echo -e "  >>  ${BOLD_s}${CMD}${END}";
        if [ $verbose -eq 0 ]; then eval ${CMD}; fi
    done
fi

######## 3. Remove rRNA

jobIDs=()
declare i=0

for abssample in ${DMSOFiles[@]};
do
    iFn=${OutDir}/2.trim/$(fnrep ${abssample}).fastq
    oFn=${OutDir}/3.rem_rRNA/$(fnrep ${abssample}).fastq
    oSAM=${OutDir}/3.rem_rRNA/$(fnrep ${abssample}).sam
    CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} -n ${Cores} \"icSHAPE-pipe cleanFq -i ${iFn} -o ${oFn} -x ${OutDir}/rRNA_index/rRNA -p ${Cores} --mode Local --sam ${oSAM}\""
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then 
        info=$(eval ${CMD})
        jobIDs+=("$(job_id $info)")
        i=$[$i+1]
    fi
done

for abssample in ${NAIFiles[@]};
do
    iFn=${OutDir}/2.trim/$(fnrep ${abssample}).fastq
    oFn=${OutDir}/3.rem_rRNA/$(fnrep ${abssample}).fastq
    oSAM=${OutDir}/3.rem_rRNA/$(fnrep ${abssample}).sam
    CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} -n ${Cores} \"icSHAPE-pipe cleanFq -i ${iFn} -o ${oFn} -x ${OutDir}/rRNA_index/rRNA -p ${Cores} --mode Local --sam ${oSAM}\""
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then 
        info=$(eval ${CMD})
        jobIDs+=("$(job_id $info)")
        i=$[$i+1]
    fi
done

for jobid in ${jobIDs[@]}; 
do
    CMD="bwait $jobid"
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then
        echo -e "remove_abundant_RNA bwait $jobid";
        eval ${CMD}
    fi
done

######## 5. Map to smallRNA

jobIDs=()
declare i=0

for abssample in ${DMSOFiles[@]};
do
    iFn=${OutDir}/3.rem_rRNA/$(fnrep ${abssample}).fastq
    oFn=${OutDir}/4.mapSmallRNA/$(fnrep ${abssample}).fastq
    oSAM=${OutDir}/4.mapSmallRNA/$(fnrep ${abssample}).sam
    CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} -n ${Cores} \"icSHAPE-pipe cleanFq -i ${iFn} -o ${oFn} -x ${OutDir}/smallRNA/smallRNA -p ${Cores} --mode Local --sam ${oSAM}\""
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then 
        info=$(eval ${CMD})
        jobIDs+=("$(job_id $info)")
        i=$[$i+1]
    fi
done

for abssample in ${NAIFiles[@]};
do
    iFn=${OutDir}/3.rem_rRNA/$(fnrep ${abssample}).fastq
    oFn=${OutDir}/4.mapSmallRNA/$(fnrep ${abssample}).fastq
    oSAM=${OutDir}/4.mapSmallRNA/$(fnrep ${abssample}).sam
    CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} -n ${Cores} \"icSHAPE-pipe cleanFq -i ${iFn} -o ${oFn} -x ${OutDir}/smallRNA/smallRNA -p ${Cores} --mode Local --sam ${oSAM}\""
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then 
        info=$(eval ${CMD})
        jobIDs+=("$(job_id $info)")
        i=$[$i+1]
    fi
done

for jobid in ${jobIDs[@]}; 
do
    CMD="bwait $jobid"
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then
        echo -e "remove_smallRNA bwait $jobid";
        eval ${CMD}
    fi
done

######## 4. Map to genome

jobIDs=()
declare i=0

for abssample in ${DMSOFiles[@]};
do
    iFn=${OutDir}/4.mapSmallRNA/$(fnrep ${abssample}).fastq
    oFn=${OutDir}/5.mapGenome/$(fnrep ${abssample})
    CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} -n ${Cores} \"icSHAPE-pipe mapGenome -i ${iFn} -o ${oFn} -x ${GenomeIndex} -p ${Cores} --noMut5\""
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then 
        info=$(eval ${CMD})
        jobIDs+=("$(job_id $info)")
        i=$[$i+1]
    fi
done

for abssample in ${NAIFiles[@]};
do
    iFn=${OutDir}/4.mapSmallRNA/$(fnrep ${abssample}).fastq
    oFn=${OutDir}/5.mapGenome/$(fnrep ${abssample})
    CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} -n ${Cores} \"icSHAPE-pipe mapGenome -i ${iFn} -o ${oFn} -x ${GenomeIndex} -p ${Cores} --noMut5\""
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then 
        info=$(eval ${CMD})
        jobIDs+=("$(job_id $info)")
        i=$[$i+1]
    fi
done

for jobid in ${jobIDs[@]}; 
do
    CMD="bwait $jobid"
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then
        echo -e "map_Genome bwait $jobid";
        eval ${CMD}
    fi
done

######## 5. Calculate FPKM

jobIDs=()
declare i=0

declare RPKMCalcFn=()
if [[ ${#DMSOFiles[@]} -eq 0 ]]; then
    RPKMCalcFn=${NAIFiles[@]} 
else
    RPKMCalcFn=${DMSOFiles[@]} 
fi

for abssample in ${RPKMCalcFn[@]};
do
    iFn=${OutDir}/5.mapGenome/$(fnrep ${abssample}).sorted.bam
    oFn=${OutDir}/6.calcFPKM/$(fnrep ${abssample})
    CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} -n ${Cores} \"icSHAPE-pipe calcFPKM -i ${iFn} -o ${oFn} -G ${GenomeGTF} -p ${Cores}\""
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then 
        info=$(eval ${CMD})
        jobIDs+=("$(job_id $info)")
        i=$[$i+1]
    fi
done

for jobid in ${jobIDs[@]}; 
do
    CMD="bwait $jobid"
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then
        echo -e "calcFPKM bwait $jobid";
        eval ${CMD}
    fi
done

######## 6. Sam file to tab

# genome

jobIDs=()
declare i=0

for abssample in ${DMSOFiles[@]};
do
    iFn=${OutDir}/5.mapGenome/$(fnrep ${abssample}).sorted.bam
    oFn=${OutDir}/7.sam2tab/$(fnrep ${abssample}).tab
    CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} \"icSHAPE-pipe sam2tab -in ${iFn} -out ${oFn}\""
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then 
        info=$(eval ${CMD})
        jobIDs+=("$(job_id $info)")
        i=$[$i+1]
    fi
    
    iFn=${OutDir}/3.rem_rRNA/$(fnrep ${abssample}).sam
    oFn=${OutDir}/7.sam2tab/rRNA_$(fnrep ${abssample}).tab
    CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} \"icSHAPE-pipe sam2tab -in ${iFn} -out ${oFn}\""
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then 
        info=$(eval ${CMD})
        jobIDs+=("$(job_id $info)")
        i=$[$i+1]
    fi

    iFn=${OutDir}/4.mapSmallRNA/$(fnrep ${abssample}).sam
    oFn=${OutDir}/7.sam2tab/smallRNA_$(fnrep ${abssample}).tab
    CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} \"icSHAPE-pipe sam2tab -in ${iFn} -out ${oFn}\""
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then 
        info=$(eval ${CMD})
        jobIDs+=("$(job_id $info)")
        i=$[$i+1]
    fi
done

for abssample in ${NAIFiles[@]};
do
    iFn=${OutDir}/5.mapGenome/$(fnrep ${abssample}).sorted.bam
    oFn=${OutDir}/7.sam2tab/$(fnrep ${abssample}).tab
    CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} \"icSHAPE-pipe sam2tab -in ${iFn} -out ${oFn}\""
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then 
        info=$(eval ${CMD})
        jobIDs+=("$(job_id $info)")
        i=$[$i+1]
    fi

    iFn=${OutDir}/3.rem_rRNA/$(fnrep ${abssample}).sam
    oFn=${OutDir}/7.sam2tab/rRNA_$(fnrep ${abssample}).tab
    CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} \"icSHAPE-pipe sam2tab -in ${iFn} -out ${oFn}\""
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then 
        info=$(eval ${CMD})
        jobIDs+=("$(job_id $info)")
        i=$[$i+1]
    fi

    iFn=${OutDir}/4.mapSmallRNA/$(fnrep ${abssample}).sam
    oFn=${OutDir}/7.sam2tab/smallRNA_$(fnrep ${abssample}).tab
    CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} \"icSHAPE-pipe sam2tab -in ${iFn} -out ${oFn}\""
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then 
        info=$(eval ${CMD})
        jobIDs+=("$(job_id $info)")
        i=$[$i+1]
    fi
done

for jobid in ${jobIDs[@]}; 
do
    CMD="bwait $jobid"
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then
        echo -e "sam2tab bwait $jobid";
        eval ${CMD}
    fi
done


######## 7. Calculate SHAPE score

# genome

declare rRNA_D_str=""
declare rRNA_N_str=""
declare smallRNA_D_str=""
declare smallRNA_N_str=""
declare D_str=""
declare N_str=""

for abssample in ${DMSOFiles[@]};
do
    D_str=${OutDir}/7.sam2tab/$(fnrep ${abssample}).tab","${D_str}
    rRNA_D_str=${OutDir}/7.sam2tab/rRNA_$(fnrep ${abssample}).tab","${rRNA_D_str}
    smallRNA_D_str=${OutDir}/7.sam2tab/smallRNA_$(fnrep ${abssample}).tab","${smallRNA_D_str}
done
D_str=$(echo ${D_str} | sed 's/,$//')
rRNA_D_str=$(echo ${rRNA_D_str} | sed 's/,$//')
smallRNA_D_str=$(echo ${smallRNA_D_str} | sed 's/,$//')

if [[ $D_str != "" ]]; then D_str="-D "${D_str}; fi
if [[ $rRNA_D_str != "" ]]; then rRNA_D_str="-D "${rRNA_D_str}; fi
if [[ $smallRNA_D_str != "" ]]; then smallRNA_D_str="-D "${smallRNA_D_str}; fi

declare MODE=calcSHAPE
if [[ ${#DMSOFiles[@]} -eq 0 ]]; then MODE=calcSHAPENoCont; fi

for abssample in ${NAIFiles[@]};
do
    N_str=${OutDir}/7.sam2tab/$(fnrep ${abssample}).tab","${N_str}
    rRNA_N_str=${OutDir}/7.sam2tab/rRNA_$(fnrep ${abssample}).tab","${rRNA_N_str}
    smallRNA_N_str=${OutDir}/7.sam2tab/smallRNA_$(fnrep ${abssample}).tab","${smallRNA_N_str}
done
N_str=$(echo ${N_str} | sed 's/,$//')
rRNA_N_str=$(echo ${rRNA_N_str} | sed 's/,$//')
smallRNA_N_str=$(echo ${smallRNA_N_str} | sed 's/,$//')

jobIDs=()
declare i=0

CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} \"icSHAPE-pipe ${MODE} \
    ${D_str} \
    -N ${N_str} \
    -size ${GenomeIndex}/chrNameLength.txt \
    -ijf ${GenomeIndex}/sjdbList.fromGTF.out.tab \
    -genome ${GenomeFasta} \
    -bases A,T,C,G \
    -out ${OutDir}/8.calcGenomeSHAPE/shape.gTab\""

echo -e "  >>  ${BOLD_s}${CMD}${END}";
if [ $verbose -eq 0 ]; then 
    info=$(eval ${CMD})
    jobIDs+=("$(job_id $info)")
    i=$[$i+1]
fi

CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} \"icSHAPE-pipe ${MODE} \
    ${rRNA_D_str} \
    -N ${rRNA_N_str} \
    -size ${OutDir}/rRNA_index/rRNA.len \
    -genome ${rRNAFasta} \
    -bases A,T,C,G \
    -out ${OutDir}/8.calcGenomeSHAPE/rRNA_shape.gTab \
    -non-sliding\""

echo -e "  >>  ${BOLD_s}${CMD}${END}";
if [ $verbose -eq 0 ]; then 
    info=$(eval ${CMD})
    jobIDs+=("$(job_id $info)")
    i=$[$i+1]
fi


CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} \"icSHAPE-pipe ${MODE} \
    ${smallRNA_D_str} \
    -N ${smallRNA_N_str} \
    -size ${OutDir}/smallRNA/smallRNA.len \
    -genome ${OutDir}/smallRNA/smallRNA.fa \
    -bases A,T,C,G \
    -out ${OutDir}/8.calcGenomeSHAPE/smallRNA_shape.gTab \
    -non-sliding\""

echo -e "  >>  ${BOLD_s}${CMD}${END}";
if [ $verbose -eq 0 ]; then 
    info=$(eval ${CMD})
    jobIDs+=("$(job_id $info)")
    i=$[$i+1]
fi

for jobid in ${jobIDs[@]}; 
do
    CMD="bwait $jobid"
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then
        echo -e "calculate SHAPE score bwait $jobid";
        eval ${CMD}
    fi
done

######## 8. generate transcriptome-based SHAPE score and RTBD

declare rpkm_str=""

for abssample in ${RPKMCalcFn[@]};
do
    rpkm_str=${OutDir}/6.calcFPKM/$(fnrep ${abssample})/isoforms.fpkm_tracking","${rpkm_str}
done
rpkm_str=$(echo ${rpkm_str} | sed 's/,$//')

# genome

jobIDs=()
declare i=0

CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} -n ${Cores} \"icSHAPE-pipe genSHAPEToTransSHAPE \
    -i ${OutDir}/8.calcGenomeSHAPE/shape.gTab \
    -o ${OutDir}/9.transSHAPE/final.shape \
    -g ${OutDir}/GTF/Anno.genomeCoor.bed \
    -p ${Cores} \
    -r ${rpkm_str}\""

echo -e "  >>  ${BOLD_s}${CMD}${END}";
if [ $verbose -eq 0 ]; then 
    info=$(eval ${CMD})
    jobIDs+=("$(job_id $info)")
    i=$[$i+1]
fi

for jobid in ${jobIDs[@]}; 
do
    CMD="bwait $jobid"
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then
        echo -e "genomeSHAPE2transSHAPE bwait $jobid";
        eval ${CMD}
    fi
done

# rRNA

jobIDs=()
declare i=0

CMD="icSHAPE-pipe genSHAPEToTransSHAPE \
    -i ${OutDir}/8.calcGenomeSHAPE/rRNA_shape.gTab \
    -o ${OutDir}/9.transSHAPE/final.shape \
    -s ${OutDir}/rRNA_index/rRNA.len \
    -p 1 \
    --app"

echo -e "  >>  ${BOLD_s}${CMD}${END}";
if [ $verbose -eq 0 ]; then eval ${CMD}; fi

# smallRNA

CMD="icSHAPE-pipe genSHAPEToTransSHAPE \
    -i ${OutDir}/8.calcGenomeSHAPE/smallRNA_shape.gTab \
    -o ${OutDir}/9.transSHAPE/final.shape \
    -s ${OutDir}/smallRNA/smallRNA.len \
    -p 1 \
    --app"

echo -e "  >>  ${BOLD_s}${CMD}${END}";
if [ $verbose -eq 0 ]; then eval ${CMD}; fi

######## 10. Visualization

jobIDs=()
declare i=0

CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} -n ${Cores} \"icSHAPE-pipe genSHAPEToBedGraph -i ${OutDir}/8.calcGenomeSHAPE/shape.gTab -o ${OutDir}/10.bedGraph/ -c 200\""
echo -e "  >>  ${BOLD_s}${CMD}${END}";
if [ $verbose -eq 0 ]; then 
    info=$(eval ${CMD})
    jobIDs+=("$(job_id $info)")
    i=$[$i+1]
fi


##########################
##### Replicate control
##########################

######## 1. Reads distribution statistic


declare sample_str=""
declare tag_str=""
for abssample in ${DMSOFiles[@]};
do
    h=$(fnrep ${abssample})
    if [ "${isCleaned}" == "true" ]; then
        sample_str=${sample_str}" -@ ${abssample},${OutDir}/3.rem_rRNA/${h}.fastq,${OutDir}/4.mapSmallRNA/${h}.fastq,${OutDir}/5.mapGenome/${h}.unsorted.bam"
    else
        sample_str=${sample_str}" -@ ${abssample},${OutDir}/1.readCollapse/${h}.fastq,${OutDir}/2.trim/${h}.fastq,${OutDir}/3.rem_rRNA/${h}.fastq,${OutDir}/4.mapSmallRNA/${h}.fastq,${OutDir}/5.mapGenome/${h}.unsorted.bam"
    fi
    
    tag_str=${tag_str}","${h}
done

for abssample in ${NAIFiles[@]};
do
    h=$(fnrep ${abssample})
    if [ "${isCleaned}" == "true" ]; then
        sample_str=${sample_str}" -@ ${abssample},${OutDir}/3.rem_rRNA/${h}.fastq,${OutDir}/4.mapSmallRNA/${h}.fastq,${OutDir}/5.mapGenome/${h}.unsorted.bam"
    else
        sample_str=${sample_str}" -@ ${abssample},${OutDir}/1.readCollapse/${h}.fastq,${OutDir}/2.trim/${h}.fastq,${OutDir}/3.rem_rRNA/${h}.fastq,${OutDir}/4.mapSmallRNA/${h}.fastq,${OutDir}/5.mapGenome/${h}.unsorted.bam"
    fi

    tag_str=${tag_str}","${h}
done

tag_str=$(echo $tag_str | sed 's/^,//')

CMD=""
if [ "${isCleaned}" == "true" ]; then
    CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} \"icSHAPE-pipe readDistributionStatistic \
        ${sample_str} \
        --sampletag ${tag_str} \
        --processtag Raw,rRNA,smallRNA,mapGenome \
        -o ${OutDir}/11.quality_control/read_map_statistics\""
else
    CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} \"icSHAPE-pipe readDistributionStatistic \
        ${sample_str} \
        --sampletag ${tag_str} \
        --processtag Raw,Unique,Trim,rRNA,smallRNA,mapGenome \
        -o ${OutDir}/11.quality_control/read_map_statistics\""
fi

echo -e "  >>  ${BOLD_s}${CMD}${END}";
if [ $verbose -eq 0 ]; then 
    info=$(eval ${CMD})
    jobIDs+=("$(job_id $info)")
    i=$[$i+1]
fi


######## 2. Reads map statistic

for abssample in ${DMSOFiles[@]};
do
    h=$(fnrep ${abssample})
    CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} \"icSHAPE-pipe samStatistics -i ${OutDir}/5.mapGenome/${h}.unsorted.bam -o ${OutDir}/11.quality_control/${h}_SamMap -g ${OutDir}/GTF/Anno.genomeCoor.bed --fast\""
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then 
        info=$(eval ${CMD})
        jobIDs+=("$(job_id $info)")
        i=$[$i+1]
    fi
done

for abssample in ${NAIFiles[@]};
do
    h=$(fnrep ${abssample})
    CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} \"icSHAPE-pipe samStatistics -i ${OutDir}/5.mapGenome/${h}.unsorted.bam -o ${OutDir}/11.quality_control/${h}_SamMap -g ${OutDir}/GTF/Anno.genomeCoor.bed --fast\""
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then 
        info=$(eval ${CMD})
        jobIDs+=("$(job_id $info)")
        i=$[$i+1]
    fi
done

for jobid in ${jobIDs[@]}; 
do
    CMD="bwait $jobid"
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then
        echo -e "map statistics bwait $jobid";
        eval ${CMD}
    fi
done


######## 3. RT replicates


declare tab_str=""
for abssample in ${DMSOFiles[@]};
do
    h=$(fnrep ${abssample})
    tab_str=${tab_str}","${OutDir}"/7.sam2tab/"${h}".tab"
done

for abssample in ${NAIFiles[@]};
do
    h=$(fnrep ${abssample})
    tab_str=${tab_str}","${OutDir}"/7.sam2tab/"${h}".tab"
done

tab_str=$(echo $tab_str | sed 's/^,//')

jobIDs=()
declare i=0

CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} \"icSHAPE-pipe countRT \
    -in ${tab_str} \
    -ijf ${GenomeIndex}/sjdbList.fromGTF.out.tab \
    -size ${GenomeIndex}/chrNameLength.txt \
    -out ${OutDir}/11.quality_control/RT_count.gTab \
    -genome ${GenomeFasta}\""

echo -e "  >>  ${BOLD_s}${CMD}${END}";
if [ $verbose -eq 0 ]; then 
    info=$(eval ${CMD})
    jobIDs+=("$(job_id $info)")
    i=$[$i+1]
fi

for jobid in ${jobIDs[@]}; 
do
    CMD="bwait $jobid"
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then
        echo -e "countRT bwait $jobid";
        eval ${CMD}
    fi
done

jobIDs=()
declare i=0

CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} \"icSHAPE-pipe plotGenomeRTRepCor -i ${OutDir}/11.quality_control/RT_count.gTab -o ${OutDir}/11.quality_control/RTRepCoorelation\""
echo -e "  >>  ${BOLD_s}${CMD}${END}";
if [ $verbose -eq 0 ]; then 
    info=$(eval ${CMD})
    jobIDs+=("$(job_id $info)")
    i=$[$i+1]
fi

for jobid in ${jobIDs[@]}; 
do
    CMD="bwait $jobid"
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then
        echo -e "plotGenomeRTRepCor bwait $jobid";
        eval ${CMD}
    fi
done


######## 4. Plot genome-based SHAPE distribution

jobIDs=()
declare i=0

CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} \"icSHAPE-pipe plotGenomeSHAPEdist -i ${OutDir}/8.calcGenomeSHAPE/shape.gTab -o ${OutDir}/11.quality_control/GenomeSHAPEDist\""
echo -e "  >>  ${BOLD_s}${CMD}${END}";
if [ $verbose -eq 0 ]; then 
    info=$(eval ${CMD})
    jobIDs+=("$(job_id $info)")
    i=$[$i+1]
fi

######## 5. transcript-base SHAPE statistics

CMD="bsub -q ${LSFqueue} -e ${errorFn} -o ${logFn} \"icSHAPE-pipe transSHAPEStatistics -i ${OutDir}/9.transSHAPE/final.shape -g ${OutDir}/GTF/Anno.genomeCoor.bed -o ${OutDir}/11.quality_control/transSHAPE\""
echo -e "  >>  ${BOLD_s}${CMD}${END}";
if [ $verbose -eq 0 ]; then 
    info=$(eval ${CMD})
    jobIDs+=("$(job_id $info)")
    i=$[$i+1]
fi

for jobid in ${jobIDs[@]}; 
do
    CMD="bwait $jobid"
    echo -e "  >>  ${BOLD_s}${CMD}${END}";
    if [ $verbose -eq 0 ]; then
        echo -e "SHAPE distribution bwait $jobid";
        eval ${CMD}
    fi
done

echo -e "${RED_s}Success!${END}"

