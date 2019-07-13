

#ifndef ALIGN_H
#define ALIGN_H

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <unordered_map>
#include <fstream>
#include <memory>
#include <exception>
#include <algorithm>
#include <math.h>

#include "pan_type.h"
#include "sstructure.h"

namespace pan{

using std::string;
using std::vector;
using std::unordered_map;
using std::istream;
using std::ostream;
using std::istringstream;
using std::pair;
using std::ifstream;
using std::runtime_error;
using std::min;
using std::max;
using std::shared_ptr;
using std::cout;
using std::endl;

const uLONG UNDEFINED_COOR = -1UL;


struct Sto_Record
{
    string chr_id;
    uLONG align_length = 0;
    uLONG seq_length = 0;

    string align_seq;       /* seq with - */
    string seq;             /* seq without - */

    //istream &operator >>(istream &, species_align &);
    //operator bool(){ return species_name.empty() ? false : true; }
    Sto_Record &operator+=(const Sto_Record&);
};

//istream &operator>>(istream &IN, Sto_Record &record);
bool read_an_sto_record(istream &IN, Sto_Record &sto_record);
ostream &operator<<(ostream &OUT, Sto_Record &align);

using MapStringStoRecord = MapStringT<shared_ptr<Sto_Record>>;

class Multi_Align
{
public:
    Multi_Align(const string &sto_file):sto_file(sto_file) { read_in_align(); }
    Multi_Align(const MapStringString &genome) { load_in_align(genome); }

    ~Multi_Align(){ if(p_ss) delete p_ss; }

    uINT size() const { return capacity; }
    uLONG length() const { return align_length; }
    const StringArray& keys() const { return chr_ids; }
    inline const Sto_Record& get_sto_record(const string &chr_id) const;
    inline bool has(const string& chr_id) const;
    inline const string &get_annotation(const string &chr_id) const;
    inline const string &get_align_seq(const string &chr_id) const;
    inline const string &get_raw_seq(const string &chr_id) const;

    // input 0-based output 1-based
    uLONG raw_coor_to_align_coor(const string &chr_id, uLONG coor) const;
    // input 0-based output 1-based
    uLONG align_coor_to_raw_coor(const string &chr_id, uLONG coor) const;

    inline string get_chr_align_sub_seq(const string &chr_id, uLONG start, uLONG length=string::npos) const;
    inline string get_chr_align_raw_seq(const string &chr_id, uLONG start, uLONG length=string::npos) const;

private:
    string sto_file;
    uLONG align_length = 0;
    uINT capacity = 0;

    MapStringStoRecord alignments;
    MapStringT<uLONGArray> raw_to_align_coor;
    MapStringT<uLONGArray> align_to_raw_coor;
    MapStringString sequence_annotation;

    StringArray chr_ids;

    SStructure *p_ss = nullptr;

    static void read_sto_head(istream &IN, MapStringString &annotation);
    void read_in_align();
    void load_in_align(const MapStringString &genome);
    void check_align_length();
    void build_raw_to_align_coor();
    void build_align_to_raw_coor();
};



const Sto_Record& Multi_Align::get_sto_record(const string &chr_id) const
{
    return *alignments.at(chr_id);
}

bool Multi_Align::has(const string& chr_id)const { 
    return ( find(chr_ids.cbegin(), chr_ids.cend(), chr_id) == chr_ids.end() ) ? false : true; 
}

string Multi_Align::get_chr_align_sub_seq(const string &chr_id, uLONG start, uLONG length) const
{
    return alignments.at(chr_id)->align_seq.substr(start, length);
}

string Multi_Align::get_chr_align_raw_seq(const string &chr_id, uLONG start, uLONG length) const
{
    return alignments.at(chr_id)->seq.substr(start, length);
}

const string &Multi_Align::get_annotation(const string &chr_id) const
{
    return sequence_annotation.at(chr_id);
}

const string &Multi_Align::get_align_seq(const string &chr_id) const
{
    return alignments.at(chr_id)->align_seq;
}

const string &Multi_Align::get_raw_seq(const string &chr_id) const
{
    return alignments.at(chr_id)->seq;
}

}
#endif // ALIGN_H

