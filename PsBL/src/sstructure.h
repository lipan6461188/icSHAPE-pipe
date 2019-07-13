
#ifndef SSTRUCRURE_H
#define SSTRUCRURE_H

#include "pan_type.h"

#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <unordered_map>
#include <cstring>
#include <memory>
#include <iostream>
#include <stack>

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

struct Base_Pair
{
    Base_Pair(uLONG left, uLONG right):left(left), right(right){}
    uLONG left;
    uLONG right;
};

using BpArray = vector<Base_Pair>;

class SStructure
{
public:
    SStructure(const string &dot_bracket);
    SStructure(const string &sequence, const string &dot_bracket);
    //SStructure(const BpArray &bp_array, const uLONG length);
    SStructure(const BpArray &bp_array, const uLONG length):length(length), base_pairs(bp_array){ }


    static vector<SStructure> read_ct_file(const string &file_name);

    string to_dot_bracket() const;
    void to_ct(const string &file_name) const;
    void to_ct(ostream &OUT)const ;

    void setTitle(const string &title){ this->title = title; }
    string getTitle()const { return this->title; }

    void setEnergy(const float &energy){ this->energy = energy; }
    float getEnergy()const { return this->energy; }

    void setSequence(const string &sequence){ if(sequence.size() == length) this->sequence = sequence; }
    string getSequence()const { return this->sequence; }

    const BpArray &get_base_pairs()const { return base_pairs; }

private:
    string title;
    float energy;
    uLONG length;
    string sequence;
    BpArray base_pairs;

    void parse_dot_bracket(const string &dot_bracket);
};



inline SStructure::SStructure(const string &dot_bracket)
{
    parse_dot_bracket(dot_bracket);
}

inline SStructure::SStructure(const string &sequence, const string &dot_bracket): sequence(sequence)
{
    if(sequence.size() != dot_bracket.size())
        throw runtime_error("sequence length != dot_bracket length");
    parse_dot_bracket(dot_bracket);
}









}
#endif