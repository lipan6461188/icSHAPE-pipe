#ifndef SLIDING_SHAPE_H
#define SLIDING_SHAPE_H


#include <param.h>
#include <sam.h>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <sstream>
#include <map>
#include <deque>
#include <numeric>

#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#include "fasta.h"

using namespace std;
using namespace pan;

#define null -1.0

extern Color::Modifier RED;
extern Color::Modifier DEF;
extern Color::Modifier YELLOW;

enum NormSampleMethod { SMART, UPPER, QUARTILE, DECILE, VIGINTILE };
enum NormCalcMethod { MEDIAN, MEAN, PEAK };

// for icSHAPE: 
//     * SUBSTRACTION -- nai_rt-sub_fac*dmso_rt
//     * DIVIDING -- div_fac*nai_bd/dmso_bd
//     * COMPLEX -- div_fac*(nai_rt-sub_fac*dmso_rt)/dmso_bd
// for smart-SHAPE
//     * SUBSTRACTION -- nai_rt/nai_bd
//     * RT -- nai_rt
enum EnrichMethod { SUBSTRACTION, DIVIDING, COMPLEX, RT, LOG }; // RT is for smart-shape
//enum STRAND{ NEG=0, POS=1 };


struct Map_Record
{
    STRAND strand = POSITIVE;
    vector<Region> regions;  

    Map_Record(const vector<Region> &regions, const STRAND &strand): 
        strand(strand), regions(regions){ }
    Map_Record(const StringArray &line_data);
    Map_Record() = default;
};




struct icSHAPE_Param
{
    NormSampleMethod norm_sample_method = DECILE;
    NormCalcMethod norm_calc_method = MEDIAN;
    uINT norm_sample_factor = 2;

    bool no_sliding = false;

    EnrichMethod enrich_method = COMPLEX;
    float winsor_factor = 0.05;
    uINT winsor_scaling = 1;

    float sub_factor = 0.25;
    float div_factor = 10;
    float add_factor = 1;

    uINT out_min_cov = 50;
    uINT min_cov = 200;
    uINT wsize = 200;
    uINT wstep = 5;
};

struct smartSHAPE_Param
{
    bool no_sliding = false;
    
    EnrichMethod enrich_method = RT;

    float sink_index = 0;

    float sub_factor = 0.25;
    float div_factor = 10;
    float add_factor = 1;

    float winsor_factor = 0.05;

    uINT out_min_cov = 50;
    uINT min_cov = 200;
    uINT wsize = 200;
    uINT wstep = 5;

    uLONG BD_ext = 0;
};

struct Junction: public Region
{
    uLONG support = 0;

    Junction(const uLONG &first, const uLONG &second): Region(first, second){ }
    uLONG width() const { return second - first + 1; }
};


using JunctionArray = vector<Junction>;

/**** Junctions ****/

// read sjdbList.fromGTF.out.tab file (in STAR index directory)
void load_junctions(const string &file_name, MapStringT<JunctionArray> &junctions, const MapStringuLONG &chr_size);

// build junction index, left means build with junction start; right means build with junction end
void buildLeftJunctionMap(JunctionArray &junctions, map<uLONG, vector<Junction*>> &junc_map, const uLONG &binsize=100000);
void buildRightJunctionMap(JunctionArray &junctions, map<uLONG, vector<Junction*>> &junc_map, const uLONG &binsize=100000);

// count junction reads
void build_junction_support(const vector<Map_Record> &record_array, JunctionArray &junctions);

// preserve junctions with max supported reads when junctions overlap
void combine_junction(JunctionArray &junctions);

// check junction overlap and combine junction when overklap
void check_overlap(JunctionArray &junctions, const uLONG &size);

// write junctions into a file
void write_junctions(const MapStringT<JunctionArray> &junctions, const string &outFile);


/**** Read Tab file (from sam file) ****/

// read a single chromosome data, input file must be sorted!!
bool read_chr(vector<Map_Record> &record_array, ifstream *hander, string &chr_id);

/**** Read chromosome size file ****/

// load chromosome size file, it can be chrNameLength.txt from STAR index directory
void load_chr_size(const string &file_name, MapStringuLONG &chr_size);

/**** calculate RT and BD ****/

// calculate chromsome RT and BD in positive strand (input record_array must be in positive strand)
void calc_chr_BDRT_Pos(uIntArray &BD, uIntArray &RT, const vector<Map_Record> &record_array, JunctionArray &junctions,const uLONG &size, const uLONG &BD_ext=0, const uLONG &binsize=100000);

// calculate chromsome RT and BD in negative strand (input record_array must be in negative strand)
void calc_chr_BDRT_Neg(uIntArray &BD, uIntArray &RT, const vector<Map_Record> &record_array, JunctionArray &junctions, const uLONG &size, const uLONG &BD_ext=0, const uLONG &binsize=100000);


/**** Cutoff functions ****/

// calculate averaged score from a float array, if the valid score number is less than half of number, it will return null
float mean_score(const FloatArray &data_array, unsigned short &c_count);

// discriminate if a region is low expressed, is the number of base exceed minBD exceed bd.size()/5, then it is not a low exp region
bool low_exp_region(const deque<float> &bd, const uLONG &minBD=10);

// if bd/rt is below a cutoff then skip a calculate
bool valid_cov(const deque<float> &bd, const deque<float> &rt, const float &min_bd=50, const float &min_rt=1.0);

/**** icSHAPE fansion (with NAI+DMSO) ****/

// calculate SHAPE score with NAI and DMSO RT/BD in non-junction regions
// score must be init
void sliding_non_junction(  const uIntArray &NAI_RT, const uIntArray &DMSO_RT,
                            const uIntArray &NAI_BD, const uIntArray &DMSO_BD,
                            const JunctionArray &junctions, FloatArray &score, 
                            const string &chr_id, const STRAND &strand, 
                            ofstream &OUT, const icSHAPE_Param &param,
                            const vector<bool> &chr_mask, const string &chr_seq);

// sliding all junctions in junction regions with NAI + DMSO RT/BD
// score must be init
void sliding_junction( const uIntArray &NAI_RT, const uIntArray &DMSO_RT,
                        const uIntArray &NAI_BD, const uIntArray &DMSO_BD,
                        const JunctionArray &junctions, FloatArray &score, 
                        const string &chr_id, const STRAND &strand, 
                        ofstream &OUT, const icSHAPE_Param &param,
                        const vector<bool> &chr_mask, const string &chr_seq);

// sliding a single junction with NAI + DMSO RT/BD
void sliding_single_junction(const uIntArray &NAI_RT, const uIntArray &DMSO_RT,
                            const uIntArray &NAI_BD, const uIntArray &DMSO_BD,
                            const uLONG &start, const uLONG &end, FloatArray &score,
                            const string &chr_id, const STRAND &strand, 
                            ofstream &OUT, const icSHAPE_Param &param,
                            const vector<bool> &chr_mask, const string &chr_seq);

// calculate SHAPE score with NAI + DMSO RT/BD
void calculate_score(const deque<float> &nai_rt, const deque<float> &nai_bd, 
                    const deque<float> &dmso_rt, const deque<float> &dmso_bd, 
                    FloatArray &scores, const icSHAPE_Param &param);

void calculate_mask_score(const deque<float> &nai_rt, const deque<float> &nai_bd, 
                    const deque<float> &dmso_rt, const deque<float> &dmso_bd, 
                    FloatArray &scores, const deque<bool> &mask, const icSHAPE_Param &param);

// calculate enrichment with DMSO and NAI
void calcEnrich(const deque<float> &dmso_bd, const deque<float> &dmso_rt,
                const deque<float> &nai_bd, const deque<float> &nai_rt,
                const float &dmso_bd_sf, const float &dmso_rt_sf, 
                const float &nai_bd_sf, const float &nai_rt_sf, 
                FloatArray &score, const icSHAPE_Param &param);

/**** smart-SHAPE fansion (with NAI only) ****/

// calculate SHAPE score with NAI RT/BD in non-junction regions,
// score must be init
void sliding_non_junction(  const uIntArray &NAI_RT, const uIntArray &NAI_BD,
                            const JunctionArray &junctions, FloatArray &score, 
                            const string &chr_id, const STRAND &strand, 
                            ofstream &OUT, const smartSHAPE_Param &param,
                            const vector<bool> &chr_mask, const string &chr_seq);

// sliding all junctions in junction regions with NAI RT/BD
// score must be init
void sliding_junction( const uIntArray &NAI_RT, const uIntArray &NAI_BD,
                        const JunctionArray &junctions, FloatArray &score, 
                        const string &chr_id, const STRAND &strand, 
                        ofstream &OUT, const smartSHAPE_Param &param,
                        const vector<bool> &chr_mask, const string &chr_seq);

// sliding a single junction with NAI RT/BD
void sliding_single_junction(const uIntArray &NAI_RT, const uIntArray &NAI_BD,
                            const uLONG &start, const uLONG &end, FloatArray &score,
                            const string &chr_id, const STRAND &strand, 
                            ofstream &OUT, const smartSHAPE_Param &param,
                            const vector<bool> &chr_mask, const string &chr_seq);

// calculate smart-SHAPE score with NAI RT and BD
void calculate_score(const deque<float> &nai_rt, const deque<float> &nai_bd, FloatArray &scores, const smartSHAPE_Param &param);

void calculate_mask_score(const deque<float> &nai_rt, const deque<float> &nai_bd, FloatArray &scores, const deque<bool> &mask, const smartSHAPE_Param &param);

// calculate enrichment with only NAI
void calcEnrich(const deque<float> &nai_bd, const deque<float> &nai_rt, FloatArray &score, const smartSHAPE_Param &param);

// sink all SHAPE that means xi <- xi - sorted_scores[sink_index*len]
void sink_score(FloatArray &score, const float &sink_index);



/**** statistics functions ****/

// winsorization a data array that means all data are assigned to 0-1
void winsorization(FloatArray &score, const float &winsor_factor=0.05);

// Get winsor upper(U) and lower(L), and normalize each raw xi to (xi-L)/(U-L)
void winsorWindow(const FloatArray &score, const float &winsor_factor, float &winsorLower, float &winsorUpper);


// calculate calcScalingFactor when input a array
float calcScalingFactor(const deque<float> &data_array, const icSHAPE_Param &param);

// Sync tab files (from sam2tab), must be sorted
bool sync_chrs(vector<ifstream *> &i_vec, StringArray &chr_ids, vector< vector<Map_Record> > &record_array);

// Check if input handle opened
void check_input_handle(ifstream &IN, const string &fn);

// Check if output handle opened
void check_output_handle(ofstream &OUT, const string &fn);

/**** load sequence ****/

inline bool build_chr_mask(qFasta &seq_holder, 
                    const string &chrID, 
                    const STRAND &strand, 
                    const vector<char> &bases, 
                    vector<bool> &chr_mask,
                    string &chr_seq)
{
    if(not seq_holder.has_chr(chrID))
    {
        return false;
    }
    chr_mask.clear();
    uLONG chrLen = seq_holder.get_chr_len(chrID);
    chr_seq = seq_holder.get_chr_subbseq(chrID, 0, string::npos, strand);

    //cout << chr_seq.substr(0, 100) << endl << endl;
    //cout << chr_seq.substr(chr_seq.size()-100) << endl;

    if(strand == NEGATIVE)
        reverse(chr_seq.begin(),chr_seq.end());

    chr_mask.assign(chrLen+1, false);
    for(uLONG i=1; i<=chrLen; i++)
        if( find(bases.cbegin(), bases.cend(), chr_seq[i-1]) != bases.cend() )
            chr_mask[i] = true;

    return true;
}



struct gTab_line
{
    string chrID;
    char strand;
    uLONG pos;
    string left_content;
};

inline bool gTab_line_cmp(const gTab_line *p1, const gTab_line *p2)
{
    if(p1->chrID < p2->chrID)
    {
        return true;
    }else if(p1->chrID > p2->chrID)
    {
        return false;
    }else{
        if(p1->strand < p2->strand)
        {
            return true;
        }else if(p1->strand > p2->strand)
        {
            return false;
        }else{
            if(p1->pos < p2->pos)
            {
                return true;
            }else{
                return false;
            }
        }
    }
}

inline void sort_gTab(const string &inFn, const string &outFn)
{
    vector<gTab_line*> gtab_vector;
    string head;
    ifstream IN;
    IN.open(inFn, ifstream::in);
    string line;
    while (getline(IN, line))
    {
        if(line[0] == '@')
        {
            head += line + '\n';
            continue;
        }

        istringstream iss(line);
        gTab_line *pline = new gTab_line;
        iss >> pline->chrID >> pline->strand >> pline->pos;
        getline (iss, pline->left_content);
        gtab_vector.push_back(pline);
    }
    IN.close();

    sort(gtab_vector.begin(), gtab_vector.end(), gTab_line_cmp);

    ofstream OUT;
    OUT.open(outFn, ofstream::out);
    OUT << head;
    for(gTab_line *pline: gtab_vector)
    {
        OUT << pline->chrID << "\t" << pline->strand << "\t" << to_string(pline->pos)  << "\t" << pline->left_content << "\n";
        delete pline;
    }
    OUT.close();
}



/*
class SeqFilterIter
{
public:
    qFasta seq_holder;
    vector<char> filtbases;
    string chrID;
    STRAND strand;
    uLONG pos=0;

    uLONG chrLen=0;
    bool good=false;

public:
    SeqFilterIter(  const string &genome_fa_file,
                    const string &chrID,
                    const STRAND &strand,
                    const vector<char> &filtbases
                  ): chrID(chrID), strand(strand), filtbases(filtbases) {
        this->seq_holder.load_fasta_file(genome_fa_file);
        if(not this->seq_holder.has_chr(chrID))
        {
            this->good = false;
        }else{
            this->chrLen = this->seq_holder.get_chr_len(chrID);
            this->good = true;
        }
    }
    
    bool fetch_base_mask(uLONG start, uLONG len, vector<bool> &mask)
    {
        if(start < pos)
        {
            this->good = false;
            throw runtime_error("start < pos is illigal");
        }
        if(start+len>=this->chrLen)
        {
            this->good = false;
            throw runtime_error("start+len>=chrLen");
        }
        mask.clear();
        mask.assign(len, false);
        string subseq = get_chr_subbseq(this->chrID, start, len, this->strand);
        
        for(int i=0;i<len;i++)
            if( find(filtbases.cbegin(), filtbases.cend(), subseq[i]) != filtbases.cend() )
                mask[i] = true;
    }
}
*/


#endif
