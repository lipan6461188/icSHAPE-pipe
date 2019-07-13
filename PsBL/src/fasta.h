#ifndef FASTA_H
#define FASTA_H

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <sstream>
#include <string>
#include <exception>
#include <cstring>
#include <algorithm>
#include <set>
#include "pan_type.h"
#include "string_split.h"
#include "exceptions.h"

namespace pan{

using std::string;

enum NUCTYPE{ DNA, RNA };
enum STRAND{ POSITIVE, NEGATIVE };

// **************************
//  Fasta class
// **************************

class Fasta
{
private:
    // Read fasta sequence and load into the object
    void read_fasta_file(const string &fastaFn, bool strict = false);

public:
    /* 
    Constract a Fasta sequence
    fastaFn             -- Input fasta file name
    strict              -- If true, raise Unexpected_Error when duplicate sequences are found
    */
    Fasta(const string &fastaFn, bool strict = false);
    
    // Default constructor
    Fasta(){};

    // Clear the object
    void clear();

    /*
    Add fasta sequence to object
    fastaFn             -- Input fasta file name
    */
    void add_fasta(const string &fastaFn, bool strict = false);

    /*
    Return the raw chrID list
    */
    StringArray get_chr_ids() const;

    /*
    Return a map of raw chrID => chrLen
    */
    MapStringuLONG get_chr_lens() const;

    /*
    Return the length of a given chrID
    key                 -- The chrID
    */
    uLONG get_chr_len(const string &chrID) const;

    /*
    Return the number of sequence
    */
    uLONG get_chr_num() const;

    /*
    Return the sequence of given chrID
    chrID               -- The chrID
    */
    const string &get_chr_seq(const string &chrID) const;

    /*
    Return the annotation of given chrID
    chrID               -- The chrID
    */
    const string &get_chr_anno(const string &chrID) const;

    /*
    Return a fragment of the sequence of given chrID
    chrID               -- The chrID
    start               -- The start position, 0-based
    len                 -- The sub-sequence length
    strand              -- Strand
    */
    string get_chr_subbseq(const string &chrID, uLONG start, uLONG len=string::npos, STRAND strand=POSITIVE) const;

    /*
    Test if a chrID is in the Fasta
    chrID                 -- The chrID
    */
    bool has_chr(const string &chrID) const { return this->sequence.find(chrID) != this->sequence.cend(); }

    // Iterators
    MapStringString::const_iterator cbegin() const { return this->sequence.cbegin(); }
    MapStringString::const_iterator cend() const { return this->sequence.cend(); }
    MapStringString::iterator begin(){ return this->sequence.begin(); }
    MapStringString::iterator end(){ return this->sequence.end(); }

private:

    // The chrID list, keep the raw order
    StringArray chr_ids;
    // A map of chrID => sequence 
    MapStringString sequence;
    // A map of chrID => annotation
    MapStringString chr_annotation;
};

// **************************
//  Quick fasta class
// **************************

typedef struct {
    uLONG chr_len;
    uLONGLONG chr_start;
    uLONG chr_line_len;
} fai_item;

// Quick fasta: a class to fetch sequence from disk rather than memory
class qFasta
{

public:
    // Read fasta sequence and load into the object
    void load_fasta_file(const string &fastaFn);

    qFasta(){};
    qFasta(const string &fastaFn);
    ~qFasta(){ if(FASTA) FASTA.close(); }

    /*
    Return the raw chrID list
    */
    StringArray get_chr_ids() const;

    /*
    Return a map of raw chrID => chrLen
    */
    MapStringuLONG get_chr_lens() const;

    /*
    Return the length of a given chrID
    key                 -- The chrID
    */
    uLONG get_chr_len(const string &chrID) const;

    /*
    Return the number of sequence
    */
    uLONG get_chr_num() const;

    /*
    Return the sequence of given chrID
    chrID               -- The chrID
    */
    string get_chr_seq(const string &chrID);

    /*
    Return a fragment of the sequence of given chrID
    chrID               -- The chrID
    start               -- The start position, 0-based
    len                 -- The sub-sequence length
    strand              -- Strand
    */
    string get_chr_subbseq(const string &chrID, uLONG start, uLONG len=string::npos/2, STRAND strand=POSITIVE);

    /*
    Test if a chrID is in the Fasta
    chrID                 -- The chrID
    */
    bool has_chr(const string &chrID) const { return this->fa_index.find(chrID) != this->fa_index.cend(); }

private:
    ifstream FASTA;
    MapStringT<fai_item> fa_index;

public:
    /*
        Build fasta index (produce a .fai file)
        fastaFn             -- Input fasta file
    */
    static void build_fasta_index(const string &fastaFn);
};


// **************************
//  Other common functions
// **************************

/*
    Return the reverse-complemenry sequence
    raw_seq             -- Raw sequence
*/
string reverse_comp(const string &raw_seq, NUCTYPE ntype=DNA);

/*
    Flat the sequence into multi-lines
    raw_seq             -- Raw sequence
    line_width          -- The width of each line
*/
string flat_seq(const string &raw_seq, size_t line_width=60);

/*
    Align two sequences
    a                   -- Input sequence 1
    b                   -- Input sequence 2
    a_aligned           -- Aligned sequence 1
    b_aligned           -- Aligned sequence 2
    gap_penalty         -- Penalty for gap
    mismatch_penalty    -- Penalty for mismatch
*/
int global_align_2_seq(const string &a,
        const string &b,
        string &a_aligned, 
        string &b_aligned,
        const int gap_penalty=5,
        const int mismatch_penalty=1);

}
#endif