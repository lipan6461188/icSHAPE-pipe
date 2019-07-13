#include "align.h"
#include <cstdlib>
#include "string_split.h"
#include <sstream>

using namespace std;

namespace pan{


/*
istream &operator>>(istream &IN, Sto_Record &record)
{
    string this_line;
    bool success = static_cast<bool>(getline(IN, this_line));
    while( success and (this_line.empty() or this_line.at(0) == '#') )
        success = static_cast<bool>(getline(IN, this_line));
    if(not success)
    {
        record.chr_id = "";
    }else{
        istringstream align_stream(this_line);
        if(not (align_stream >> record.chr_id >> record.align_seq) )
        {
            cerr << "Bad align line: " << this_line << endl;
            record.chr_id = "";
        }else{
            record.align_length = record.align_seq.size();
            for(char &seq_base: record.align_seq)
                if(seq_base != '-' and seq_base != '.' and seq_base != '_' and seq_base != '*' and seq_base != ' ')
                    record.seq.push_back(seq_base);
                else
                    seq_base = '-';
            record.seq_length = record.seq.size();
        }
    }
    return IN;
}
*/



bool read_an_sto_record(istream &IN, Sto_Record &sto_record)
{
    string this_line;
    bool success = static_cast<bool>(getline(IN, this_line));
    trim(this_line);
    while( success and (this_line.empty() or this_line.at(0) == '#') )
    {
        success = static_cast<bool>(getline(IN, this_line));
        trim(this_line);

        // #=GC SS_cons lines
        if(this_line.substr(0,4) == "#=GC")
        {
            istringstream align_stream(this_line);
            string tag, mode, ss_string;
            align_stream >> tag >> mode >> ss_string;
            if(mode == "SS_cons")
            {
                sto_record.chr_id = "SS_cons";
                sto_record.align_seq = ss_string;

                return true;
            }
        }
    }
    if(not success)
    {
        return false;
    }else if(this_line == "//")
    {
        return false;
    }else{
        istringstream align_stream(this_line);

        align_stream >> sto_record.chr_id >> sto_record.align_seq;
        /*
        string::size_type pos;
        if( (pos = sto_record.chr_id.find_first_of('/')) != string::npos )
        {
            sto_record.chr_id = sto_record.chr_id.substr(0, pos);
        }
        */
        //clog << sto_record.chr_id << endl;
        sto_record.align_length = sto_record.align_seq.size();

        for(char &seq_base: sto_record.align_seq)
            if(seq_base != '-' and seq_base != '.' and seq_base != '_' and seq_base != '*' and seq_base != ' ')
                sto_record.seq.push_back(seq_base);
            else
                seq_base = '-';

        sto_record.seq_length = sto_record.seq.size();
    }
    return true;
}

ostream &operator<<(ostream &OUT, const Sto_Record &sto_record)
{
    OUT << sto_record.chr_id << "\t" << sto_record.align_seq << "\n";
    return OUT;
}

Sto_Record &Sto_Record::operator +=(const Sto_Record& sto_record)
{
    this->align_length += sto_record.align_length;
    this->seq_length += sto_record.seq_length;
    this->align_seq += sto_record.align_seq;
    this->seq += sto_record.seq;

    return *this;
}

void Multi_Align::read_in_align()
{
    ifstream IN(sto_file, ifstream::in);
    if(not IN)
    {
        throw runtime_error(sto_file+" cannot be readable");
    }

    read_sto_head(IN, sequence_annotation);
    while(1)
    {
        shared_ptr<Sto_Record> p_sto_record(new Sto_Record);
        //IN >> *new_align;
        if(read_an_sto_record(IN, *p_sto_record))
        {
            const string chr_id = p_sto_record->chr_id;

            if(find(chr_ids.cbegin(), chr_ids.cend(), chr_id) != chr_ids.end())
            {
                *alignments.at(chr_id) += *p_sto_record;
            }else{
                chr_ids.push_back(chr_id);
                alignments[chr_id] = p_sto_record;
                ++capacity;
            }

            /*
            auto align_iter = alignments.find(p_sto_record->species_name);
            if(align_iter != alignments.cend())
            {
                *(align_iter->second) += *p_sto_record;
            }else{
                species.push_back(p_sto_record->species_name);
                alignments[p_sto_record->species_name] = p_sto_record;
                ++capacity;
            }
            */

        }else{
            break;
        }
    }
    IN.close();
    sort(chr_ids.begin(), chr_ids.end());
    
    auto ss_cons_iter = find(chr_ids.cbegin(), chr_ids.cend(), "SS_cons");
    if(ss_cons_iter != chr_ids.end())
    {
        // secondary structure in file
        chr_ids.erase(ss_cons_iter, ss_cons_iter+1);
        p_ss = new SStructure(alignments.at("SS_cons")->align_seq);
        alignments.erase("SS_cons");
    }

    if(alignments.empty())
    {
        throw runtime_error("Empty multi-alignment file: " + sto_file);
        //cerr << "FATAL Error: Null Multialignment" << endl;
        //exit(EXIT_FAILURE);
    }
   // for(int i=0;i<chr_ids.size();i++)
     //   cout << chr_ids[i] << "\t" << alignments.at(chr_ids[i])->align_length << endl;
    check_align_length();
    build_raw_to_align_coor();
    build_align_to_raw_coor();
}

void Multi_Align::load_in_align(const MapStringString &genome)
{
    // check
    if(genome.size() <= 1)
        throw runtime_error("too less sequences");

    uLONG first_len = genome.cbegin()->second.size();
    bool invalid_length = std::any_of(genome.cbegin(), genome.cend(), [&](const pair<string, string> &seq)->bool{ return first_len != seq.second.size(); });
    if(invalid_length)
        throw runtime_error("different aligned sequences");

    for(auto iter=genome.cbegin(); iter!=genome.cend(); iter++)
    {
        stringstream STREAM(iter->first+"\t"+iter->second);

        shared_ptr<Sto_Record> p_sto_record(new Sto_Record);
        if(read_an_sto_record(STREAM, *p_sto_record))
        {
            const string chr_id = p_sto_record->chr_id;
            chr_ids.push_back(chr_id);
            alignments[chr_id] = p_sto_record;
            ++capacity;
        }else{
            throw runtime_error("Impossible Error in load_in_align/align.cpp");
        }
    }

    sort(chr_ids.begin(), chr_ids.end());
    check_align_length();
    build_raw_to_align_coor();
    build_align_to_raw_coor();
}

void Multi_Align::read_sto_head(istream &IN, MapStringString &annotation)
{
    annotation.clear();

    string cur_line;

    auto last_pos = IN.tellg();
    while(getline(IN, cur_line))
    {
        if(cur_line.empty())
            continue;

        if(cur_line.at(0) == '/' and cur_line.at(1) == '/')
            break;

        if(cur_line.at(0) != '#')
        {
            IN.seekg(last_pos);
            break;
        }
        last_pos = IN.tellg();

        if( cur_line.substr(0, 4) == "#=GS" )
        {
            istringstream anno_stream(cur_line.substr(5));
            string chr_id, tag;
            anno_stream >> chr_id >> tag;
            if(tag == "DE")
            {
                getline(anno_stream, annotation[chr_id]);
                //annotation[chr_id] = anno_stream.str();
                //cout << chr_id << "\t" << annotation[chr_id] << endl;
            }
        }
    }
}

void Multi_Align::check_align_length()
{
    auto iter = alignments.cbegin();
    align_length = iter->second->align_length;
    for(; iter!=alignments.cend(); iter++)
    {
        if(align_length != iter->second->align_length)
        {
            char buffer[2000];
            sprintf(buffer, "Different Alignment Length -- %s(%lu) <==> %s(%lu)", alignments.cbegin()->first.c_str(), align_length, iter->first.c_str(), iter->second->align_length);
            throw runtime_error(buffer);
        }
    }

    /*
    for(const pair<string, shared_ptr<species_align>> &single_align: alignments )
    {
        if(align_length==0)
        {
            align_length = single_align.second->align_length;
            UNDEFINED_COOR = align_length;
        }else if(align_length != single_align.second->align_length)
        {
            cerr << "FATAL Error: Bad Alignment length: (" <<  align_length << ") and (" << single_align.second->align_length << ")" << endl;
            exit(EXIT_FAILURE);
        }
    }
    */
}






void Multi_Align::build_raw_to_align_coor()
{
    raw_to_align_coor.clear();
    for(const string &chr_id: chr_ids)
    {
        const Sto_Record &record = *(alignments.at(chr_id));
        //raw_to_align_coor[chr_id].resize(sa.seq_length);
        uLONG base_idx = 1;
        for(char base: record.align_seq)
        {
            if(base != '-')
            {
                raw_to_align_coor[chr_id].push_back(base_idx);
            }
            base_idx++;
        }
        if(raw_to_align_coor[chr_id].size() != record.seq_length)
        {
            throw runtime_error("Different Length in build_raw_to_align_coor");
            //cerr << "FATAL Error: build_raw_to_align_coor error: " << raw_to_align_coor[chr_id].size() << "\t" << sa.seq_length << endl;
            //exit(EXIT_FAILURE);
        }
    }
}

void Multi_Align::build_align_to_raw_coor()
{
    align_to_raw_coor.clear();
    for(const string &chr_id: chr_ids)
    {
        const Sto_Record &record = *(alignments.at(chr_id));
        //align_to_raw_coor[chr_id].resize(sa.align_length);
        uLONG base_idx = 1;
        for(char base: record.align_seq)
        {
            if(base != '-')
            {
                align_to_raw_coor[chr_id].push_back(base_idx);
                base_idx++;
            }else{
                align_to_raw_coor[chr_id].push_back(UNDEFINED_COOR);
            }
        }
        if(align_to_raw_coor[chr_id].size() != record.align_length)
        {
            throw runtime_error("Different Length in build_align_to_raw_coor");
            //cerr << "FATAL Error: build_align_to_raw_coor error: " << align_to_raw_coor[chr_id].size() << "\t" << sa.align_length << endl;
            //exit(EXIT_FAILURE);
        }
    }
}

uLONG Multi_Align::raw_coor_to_align_coor(const string &chr_id, uLONG coor) const
{
    /*
    if(not has(chr_id))
        return UNDEFINED_COOR;
    if(coor > raw_to_align_coor.at(chr_id).size())
        return UNDEFINED_COOR;
    */
    return raw_to_align_coor.at(chr_id).at(coor);
}

uLONG Multi_Align::align_coor_to_raw_coor(const string &chr_id, uLONG coor) const
{
    /*
    if(not has(chr_id))
        return UNDEFINED_COOR;
    if(coor >= align_length)
        return UNDEFINED_COOR;
    */
    return align_to_raw_coor.at(chr_id).at(coor);
}








}