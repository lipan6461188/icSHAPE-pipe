
#ifndef STO_H
#define STO_H

#include "pan_type.h"
#include <cstring>
#include <exception>

namespace pan{

class Sto
{

public:
    //using MapSmallCodeString = unordered_map<char*>
    Sto(const string &file_name);


private:
    const string sto_file_name;

    const string version = "STOCKHOLM 1.0";
    
    MapStringString GF_Items; // Description of fields
    MapStringString GS_Items; // 
    MapStringString GR_Items; //
    MapStringString GC_Items;
    MapStringString Sequence;

    void read_sto_file();

    void read_head(istream &IN);        // stop when occur GS
    void read_GS(istream &IN);          // stop when occur blank lines
    void read_a_block(istream &IN); // stop when occur not #GR head not id head not //
}

Sto::Sto(const string &file_name): sto_file_name(file_name)
{
    read_sto_file();
}

void Sto::read_head(istream &IN)
{
    string cur_line;
    getline(IN, cur_line);
    if(cur_line != "# STOCKHOLM 1.0")
        throw runtime_error(sto_file_name+" is not a STOCKHOLM 1.0 file");

    auto last_pos = IN.tellg();
    while(getline(IN, cur_line))
    {
        while(cur_line.empty() and getline(IN, cur_line));
        
        if(cur_line.substr(0, 4) == "#=GF")
        {
            StringArray items;
            string head, tag, annotation;

            istringstream Stream(cur_line);
            Stream >> head >> tag;
            getline(Stream, annotation);

            GF_Items[tag] += annotation;

            last_pos = IN.tellg();
        }else{
            IN.seekg(last_pos);
            return;
        }
    }
}

void Sto::read_GS(istream &IN)
{
    string cur_line;

    auto last_pos = IN.tellg();
    while(getline(IN, cur_line))
    {
        while(cur_line.empty() and getline(IN, cur_line));
        
        if(cur_line.substr(0, 4) == "#=GS")
        {
            StringArray items;
            string head, chr_id, feature,annotation;

            istringstream Stream(cur_line);
            Stream >> head >> chr_id >> feature;
            getline(Stream, annotation);

            if(feature == "DE")
            {
                auto chr_id_items = split(chr_id, '/');
                GS_Items[ chr_id_items[0] ] = annotation;
                Sequence[ chr_id_items[0] ];
            }else{
                throw runtime_error("Undefined GS Feature: "+feature);
            }

            last_pos = IN.tellg();
        }else{
            IN.seekg(last_pos);
            return;
        }
    }
}

void read_a_block(istream &IN)
{
    string cur_line;

    auto last_pos = IN.tellg();
    while(getline(IN, cur_line))
    {
        while(cur_line.empty() and getline(IN, cur_line));
        
        if(cur_line[0] != "#")
        {
            
        }else if(cur_line.substr(0, 4) == "#=GS")
        {
            StringArray items;
            string head, chr_id, feature,annotation;

            istringstream Stream(cur_line);
            Stream >> head >> chr_id >> feature;
            getline(Stream, annotation);

            if(feature == "DE")
            {
                auto chr_id_items = split(chr_id, '/');
                GS_Items[ chr_id_items[0] ] = annotation;
                Sequence[ chr_id_items[0] ];
            }else{
                throw runtime_error("Undefined GS Feature: "+feature);
            }

            last_pos = IN.tellg();
        }else{
            IN.seekg(last_pos);
            return;
        }
    }
}

void Sto::read_sto_file()
{
    ifstream IN(sto_file_name, ifstream::in);
    if(not IN)
    {
        throw runtime_error(sto_file_name+" cannot readable");
    }

    // read head
    string cur_line;
    getline(IN, cur_line);
    if(cur_line != "# STOCKHOLM 1.0")
    {
        throw runtime_error(sto_file_name+" is not a STOCKHOLM 1.0 file");
    }

    while(getline(IN, cur_line))
    {
        if(cur_line.empty())
    }
}













































}
#endif