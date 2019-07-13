#ifndef FASTA_H
#define FASTA_H

#include "pan_type.h"
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <sstream>
#include <string>
#include <exception>
#include "string_split.h"
#include <cstring>
#include <algorithm>
#include "exceptions.h"

namespace pan{

using std::vector;
using std::unordered_map;
using std::string;

struct GTF_Line
{
    enum FRAME{ NO, ZERO, ONE, TWO };

    string chr_id;
    string source;
    string feature;
    
    uLONG start = 0;
    uLONG end = 0;

    double score = 0;
    char strand;
    FRAME frame;

    MapStringString attributes;
};

class Base_GTF_Gene
{
public:
    Base_GTF_Gene() = default;
    Base_GTF_Gene(sp<GTF_Line> sp_gtf_line):sp_gtf_line(sp_gtf_line) { }
    Base_GTF_Gene(const Base_GTF_Gene &other): sp_gtf_line(other.sp_gtf_line) {  }

    virtual const string &gene_id() const = 0;
    virtual const string &gene_name() const = 0;

    const string &chr_id() const{ return sp_gtf_line->chr_id; }
    uLONG genome_start() const{ return sp_gtf_line->start; }
    uLONG genome_end() const{ return sp_gtf_line->end; }
    
    char strand() const{ return sp_gtf_line->strand; }

protected:
    sp<GTF_Line> sp_gtf_line;
};

class Base_GTF_Transcript
{
public:
    Base_GTF_Transcript() = default;
    Base_GTF_Transcript(sp<GTF_Line> sp_gtf_line):sp_gtf_line(sp_gtf_line) { }
    Base_GTF_Transcript(const Base_GTF_Transcript &other): sp_gtf_line(other.sp_gtf_line) {  }

    virtual const string &trans_id() const = 0;
    virtual const string &trans_type() const = 0;

    const string &chr_id() const{ return sp_gtf_line->chr_id; }
    uLONG genome_start() const{ return sp_gtf_line->start; }
    uLONG genome_end() const{ return sp_gtf_line->end; }

    char strand() const{ return sp_gtf_line->strand; }

    virtual const string &uniq_gene_id() const = 0; 

protected:
    sp<GTF_Line> sp_gtf_line;
};

class Base_GTF_Exon
{
public:
    Base_GTF_Exon() = default;
    Base_GTF_Exon(sp<GTF_Line> sp_gtf_line):sp_gtf_line(sp_gtf_line) { }
    Base_GTF_Exon(const Base_GTF_Exon &other): sp_gtf_line(other.sp_gtf_line) {  }

    virtual const string &trans_id() const = 0;

    const string &chr_id() const{ return sp_gtf_line->chr_id; }
    uLONG genome_start() const{ return sp_gtf_line->start; }
    uLONG genome_end() const{ return sp_gtf_line->end; }

    char strand() const{ return sp_gtf_line->strand; }

protected:
    sp<GTF_Line> sp_gtf_line;
};

class Base_GTF_CDS
{
public:
    Base_GTF_CDS() = default;
    Base_GTF_CDS(sp<GTF_Line> sp_gtf_line):sp_gtf_line(sp_gtf_line) { }
    Base_GTF_CDS(const Base_GTF_CDS &other): sp_gtf_line(other.sp_gtf_line) {  }

    virtual const string &trans_id() const = 0;

    const string &chr_id() const{ return sp_gtf_line->chr_id; }
    uLONG genome_start() const{ return sp_gtf_line->start; }
    uLONG genome_end() const{ return sp_gtf_line->end; }

    char strand() const{ return sp_gtf_line->strand; }

protected:
    sp<GTF_Line> sp_gtf_line;
};

bool operator< (const Base_GTF_Exon &exon_1, const Base_GTF_Exon &exon_2);
bool operator< (const Base_GTF_CDS &cds_1, const Base_GTF_CDS &cds_2);


struct GTF2_Parse_Key
{
    string gene_id_key_word = "gene_id";
    string transcript_id_key_word = "transcript_id";
    string transcript_type_key_word = "transcript_type";
    string gene_name_key_word = "gene_name";

    string gene_feature_key = "gene";
    string transcript_feature_key = "transcript";
    string exon_feature_key = "exon";
    string cds_feature_key = "CDS";
};

class GTF2_Exon;
class GTF2_CDS;
class GTF2_Transcript;

class GTF2_Gene: public Base_GTF_Gene
{
public:
    GTF2_Gene() = default;
    GTF2_Gene(sp<GTF_Line> sp_gtf_line, const GTF2_Parse_Key * const gtf_parser):Base_GTF_Gene(sp_gtf_line), gtf_parser(gtf_parser){ }
    GTF2_Gene(const GTF2_Gene& other):Base_GTF_Gene(other), gtf_parser(other.gtf_parser){  }

    virtual const string &gene_id() const override { return sp_gtf_line->attributes.at(gtf_parser->gene_id_key_word); }
    virtual const string &gene_name() const override { return sp_gtf_line->attributes.at(gtf_parser->gene_name_key_word); }

    vector< sp<GTF2_Transcript> > trans_list;

private:
    const GTF2_Parse_Key * gtf_parser;
};


class GTF2_Transcript: public Base_GTF_Transcript
{
public:
    GTF2_Transcript() = default;
    GTF2_Transcript(sp<GTF_Line> sp_gtf_line, const GTF2_Parse_Key * const gtf_parser):Base_GTF_Transcript(sp_gtf_line), gtf_parser(gtf_parser){ }
    GTF2_Transcript(const GTF2_Transcript& other):Base_GTF_Transcript(other), gtf_parser(other.gtf_parser){  }

    virtual const string &trans_id() const override { return sp_gtf_line->attributes.at(gtf_parser->transcript_id_key_word); }
    virtual const string &trans_type() const override { return sp_gtf_line->attributes.at(gtf_parser->transcript_type_key_word); }
    virtual const string &uniq_gene_id() const override { return sp_gtf_line->attributes.at(gtf_parser->gene_id_key_word); }

    sp<GTF2_Gene> gene;

    vector< sp<GTF2_CDS> > cds_list;
    vector< sp<GTF2_Exon> > exon_list;

private:
    const GTF2_Parse_Key * gtf_parser;
};


class GTF2_Exon: public Base_GTF_Exon
{
public:
    GTF2_Exon() = default;
    GTF2_Exon(sp<GTF_Line> sp_gtf_line, const GTF2_Parse_Key * const gtf_parser):Base_GTF_Exon(sp_gtf_line), gtf_parser(gtf_parser){ }
    GTF2_Exon(const GTF2_Exon& other):Base_GTF_Exon(other), gtf_parser(other.gtf_parser){  }

    virtual const string &trans_id() const override { return sp_gtf_line->attributes.at(gtf_parser->transcript_id_key_word); }

private:
    const GTF2_Parse_Key * gtf_parser;
};

class GTF2_CDS: public Base_GTF_CDS
{
public:
    GTF2_CDS() = default;
    GTF2_CDS(sp<GTF_Line> sp_gtf_line, const GTF2_Parse_Key * const gtf_parser):Base_GTF_CDS(sp_gtf_line), gtf_parser(gtf_parser){ }
    GTF2_CDS(const GTF2_CDS& other):Base_GTF_CDS(other), gtf_parser(other.gtf_parser){  }

    virtual const string &trans_id() const override { return sp_gtf_line->attributes.at(gtf_parser->transcript_id_key_word); }

private:
    const GTF2_Parse_Key * gtf_parser;
};



class GTF2_Gencode
{
public:
    GTF2_Gencode(GTF2_Parse_Key * const gtf_parser, const string &gtf2_file_name): gtf_parser(gtf_parser)
    {
        read_file(gtf2_file_name, gtf_parser);
    }

    ~GTF2_Gencode(){ delete gtf_parser; }

    sp<StringArray> get_gene_list() const;
    sp<StringArray> get_trans_list() const;
    sp<MapStringT< vector<sp<GTF2_Transcript> >>> get_gene_trans_map() const;

    bool has_gene(const string &gene_id) const;
    bool has_trans(const string &trans_id) const;
    const sp<GTF2_Gene> get_gene_handle(const string &gene_id) const;
    const sp<GTF2_Transcript> get_trans_handle(const string &trans_id) const;

    void write_genome_coor_annotation(ostream &OUT) const;

private:
    void read_file(const string &gtf2_file_name, const GTF2_Parse_Key * const gtf_parser);

    MapStringT< sp<GTF2_Gene> > genes;
    MapStringT< sp<GTF2_Transcript> > transcripts;
    
    //MapStringT< vector<GTF2_Exon> > cds;

    GTF2_Parse_Key * const gtf_parser;
};
















































}
#endif