#ifndef SAM_H
#define SAM_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <exception>
#include <string>
#include <cassert>
#include <algorithm>

#include "pan_type.h"
#include "string_split.h"
#include "exceptions.h"
#include "htslib.h"

namespace pan{

using std::string;
using std::vector;
using std::unordered_map;
using std::istream;
using std::ostream;
using std::istringstream;
using std::pair;

/* 
    Read Match Pattern:
        MATCH_LEFT: 30M30S
        MATCH_RIGHT: 30S30M
        MATCH_CENTER: 15S30M15S
        MATCH_TOTAL: 60M
        MATCH_SEG: 10M10S10M10S10M10S/15M30S15M
 */
enum MATCH_PATTERN{
    MATCH_LEFT, MATCH_RIGHT, MATCH_CENTER, MATCH_TOTAL, MATCH_SEG
};

struct Sam_Record;

// read information
bool read_is_mapped(const Sam_Record &read_record);
bool read_is_gapped(const Sam_Record &read_record);
bool read_is_reverse(const Sam_Record &read_record);
bool read_is_primary(const Sam_Record &read_record);
bool read_is_multimap(const vector<Sam_Record> &read_records);

/* a line of sam file */
struct Sam_Record
{
    string read_id;
    uINT flag;
    string chr_id;
    uLONG pos;
    uINT map_quanlity;
    string cigar;
    
    string read_id_next;
    uLONG pos_next;
    uINT temp_len;

    string read_seq;
    string read_quality;

    vector<string> attributes;

    string get_attr(const string &attr_name);

    char strand()const { return read_is_reverse(*this) ? '-' : '+'; }
};

using ReadPair = pair<Sam_Record, Sam_Record>;

/****** FUNCTION from Sam file *******/
bool read_a_sam_record(istream &IN, Sam_Record &read_record);                   // read a record from sam file
void write_a_sam_record(ostream &OUT, const Sam_Record &read_record);           // write a record to sam file
bool read_a_read_record(istream &IN, vector<Sam_Record> &read_records);         // read records for multi-map from sam file (slow)
void write_a_read_record(ostream &OUT, const vector<Sam_Record> &read_records); // write records for multi-map to sam file (slow)
void remove_read_mutimap(Sam_Record &read_record);                              // remove multimap flag from tag eg. 256 => 0
void add_read_mutimap(Sam_Record &read_record);                                 // add multimap flag to tag  eg. 0 => 256
void remove_read_reversed(Sam_Record &read_record);                             // remove reverse flag from tag eg. 16 => 0
void add_read_reversed(Sam_Record &read_record);                                // add multiple flag to tag eg. 0 => 16

/****** FUNCTION from Bam file *******/
bool read_a_sam_record(BGZF* fn_hd, bam_hdr_t *hdr, Sam_Record &read_record);   // read a record from bam file
bool read_a_read_record(BGZF* fn_hd, bam_hdr_t *hdr, vector<Sam_Record> &read_records, Sam_Record* &p_cache);
                                                                                // read records for multi-map from bam file (slow)

/****** FUNCTION from pair-end reads *******/
inline bool read_is_paired(const Sam_Record &read_record){ return read_record.flag & 1; }
inline bool read_is_mate_reversed(const Sam_Record &read_record){ return read_record.flag & 32; }
inline bool read_is_first_in_pair(const Sam_Record &read_record){ return read_record.flag & 64; }

// if two reads are paired
inline bool is_paired_reads(const Sam_Record &read_1, const Sam_Record &read_2)
{
    return (read_1.read_id == read_2.read_id) and (read_1.flag & 1) and (read_2.flag & 1) and ((read_1.flag & 64 and read_2.flag & 128) or (read_1.flag & 128 and read_2.flag & 64)) and (read_1.pos_next == read_2.pos);
}

ostream& operator<<(ostream &OUT, const Sam_Record &read_record);
ostream& operator<<(ostream &OUT, const vector<Sam_Record> &read_records);
ostream& operator<<(ostream &OUT, const ReadPair &read_pair);

// Compare two record
bool operator<(const Sam_Record& record_1, const Sam_Record& record_2);

struct Sam_Head
{
    MapStringT<uLONG> trans_len;

    // TODO: add more structured data structure of sam head
    StringArray head_list;
};

/****** FUNCTION from Sam file *******/
bool read_sam_head(istream &IN, Sam_Head &sam_head);
void write_sam_head(ostream &OUT, const Sam_Head &sam_head);

/****** FUNCTION from Bam file *******/
bool read_sam_head(bam_hdr_t *hdr, Sam_Head &sam_head);

ostream& operator<<(ostream &OUT, const Sam_Head &read_records);

// filter some reads
void filter_unmapped_record(vector<Sam_Record> &read_records);
void filter_ungapped_record(vector<Sam_Record> &read_records);
void filter_gapped_record(vector<Sam_Record> &read_records);
void filter_reversed_record(vector<Sam_Record> &read_records);
void filter_secondary_record(vector<Sam_Record> &read_records);

// filter reads from pair-end
//void filter_second_in_pair(vector<Sam_Record> &read_records);

//*********************************
//*********************************
//  Functions for juncted reads
//*********************************
//*********************************

// give a cigar string, return all splited items
void split_cigar(const string &cigar, vector<int> &cigarLen, vector<char> &cigarAlpha);
// give a cigar and a startpos of cigar, return all match region in the cigar
void get_global_match_region(const string &cigar, uLONG startPos, RegionArray &matchRegion);
// return all match region in the cigar
void get_local_match_region(const string &cigar, RegionArray &matchRegion);
// reverse cigar code
string reverse_cigar(const string &cigar);

// get larger soft-clip from 5' or 3' 
uINT read_max_softclip(const string &cigar);

/*  get cigar code of matched part 
    Example:
        cout << get_match_cigar("30S10M1D20M10S") << endl;
        // 10M1D20M
*/
string get_match_cigar(const string &raw_cigar);

/*  get cigar code of region part 
    Example:
        cout << get_region_cigar("30S10M30S", Region(11, 40)) << endl;
        // 20S10M10S
*/
string get_region_cigar(const string &raw_cigar, const Region &region);

string::size_type remove_left_D(string &raw_cigar);
string::size_type remove_right_D(string &raw_cigar);

/*  get the match pattern of a reads
    1. 60M          return: MATCH_TOTAL
    2. 10S50M       return: MATCH_RIGHT
    3. 50M10S       return: MATCH_LEFT
    4. 6S34M20S     return: MATCH_LEFT
    5. 10S30M10S    return: MATCH_CENTER
*/
MATCH_PATTERN cigar_match_pattern(const string &raw_cigar);



//  ================ Group Two Reads into One Read ================

enum GROUP_STATUS{
    GROUP_SUCCESS=0,    // success
    GROUP_DIFF_ID,      // map different read id
    GROUP_UNMMAPED,     // unmmaped read
    GROUP_DIFF_STRAND,  // map to different strand
    GROUP_DIFF_CHR,     // map to different chromosome
    GROUP_MULTI_MATCH,  // read map to too many regions
    GROUP_OVERLAP_GLOBAL, // two map to overlaped global region
    GROUP_ABNORMAL_MATCH_PATTERN, // two reads are not MATCH_LEFT/MATCH_RIGHT pair
    GROUP_INVALID_BORDER          // unexpected start-end border
};

//  ================ Group Two Reads ================

/*  group two reads into one read
    two read should:
        0. same read id
        1. map to same chromosome
        2. map to sense or anti-map
        3. valid to group(global overlap is not allowed, local overlap is allowed)

    Return: if success, return true
 */
GROUP_STATUS group_read_pair(const ReadPair &raw_read_pair, Sam_Record &grouped_sam_record);


//  ================ Trim Two Reads ================

enum TRIM_STATUS{
    TRIM_SUCCESS=0,    // success
    TRIM_DIFF_ID,      // map different read id
    TRIM_UNMMAPED,     // unmmaped read
    TRIM_MULTI_MATCH,  // read map to too many regions
    TRIM_OVERLAP_GLOBAL, // two map to overlaped global region
    TRIM_ABNORMAL_MATCH_PATTERN, // two reads are not MATCH_LEFT/MATCH_RIGHT pair
    TRIM_INVALID_BORDER,         // unexpected start-end border
    TRIM_OVERLAP_INCLUDE         // too much overlap
};

/*  trim two reads so that they can duplex
    two read should:
        0. same read id
        1. valid to group(global overlap is not allowed, local overlap is allowed)

    Return: if success, return true
 */
static pair<uLONG, uLONG> __unuseful_bias;
TRIM_STATUS trim_read_pair( const ReadPair &raw_read_pair, 
                            ReadPair &trimmed_read_pair, 
                            pair<uLONG, uLONG> &bias=__unuseful_bias);


/***************************************/
/***************************************/
//            Calculate RPKM 
//  TODO: check and finilize these functions...
/***************************************/
/***************************************/

struct RPKM_PARAM
{
    ostream *p_verbose_out = nullptr;

    bool preserve_multimap = true;
    bool preserve_reverse_map = false;
};

/*
    paired-end:
        E00477:208:HG7F5CCXY:4:1102:4239:11084  99      mate reverse strand & first in pair
        E00477:208:HG7F5CCXY:4:1102:4239:11084  147     read reverse strand & second in pair     
        E00477:208:HG7F5CCXY:4:1102:4239:11084  355     mate reverse strand & first in pair & not primary alignment
        E00477:208:HG7F5CCXY:4:1102:4239:11084  403     read reverse strand & second in pair & not primary alignment

        E00477:208:HG7F5CCXY:4:1102:5416:11084  163     mate reverse strand & second in pair
        E00477:208:HG7F5CCXY:4:1102:5416:11084  83      read reverse strand & first in pair
        E00477:208:HG7F5CCXY:4:1102:13595:11084 419     mate reverse strand & second in pair & not primary alignment
        E00477:208:HG7F5CCXY:4:1102:13595:11084 339     read reverse strand & first in pair & not primary alignment

        E00477:208:HG7F5CCXY:4:1101:1560:30140  73      read paired & mate unmapped & first in pair
        E00477:208:HG7F5CCXY:4:1101:1580:19768  137     read paired & mate unmapped & second in pair
        E00477:208:HG7F5CCXY:4:1101:1580:26624  89      read paired & mate unmapped & read reverse strand &  first in pair      
*/

void statistic_sam_abundance(const string &sam_file, 
        MapStringT<uLONG> &uniq_mapped_reads, 
        MapStringT<double> &multi_mapped_reads,
        MapStringT<uLONG> &ref_len,
        const RPKM_PARAM &param);

inline double rpkm_func(const double &total_mapped_reads, const double &map_trans_reads, const uLONG &trans_len)
{
    return map_trans_reads * 1000000000 / (trans_len * total_mapped_reads);
}

void calc_rpkm(const string &sam_file, MapStringT<double> &rpkm, const RPKM_PARAM &param);


}
#endif