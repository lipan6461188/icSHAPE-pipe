#include <paris.h>
#include <param.h>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <sstream>

#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#include "version.h"

using namespace std;
using namespace pan;

#define WARNING "The input sam file should be aligned with end-to-end mode"

Color::Modifier RED(Color::FG_RED);
Color::Modifier DEF(Color::FG_DEFAULT);
Color::Modifier YELLOW(Color::FG_YELLOW);

void print_usage()
{
    char buff[2000];
    const char *help_info = 
            "sam2tab - covert sam file to a tab-seperated file\n"
            "=============================================================\n"
            "\e[1mUSAGE:\e[0m\n"
            "\tsam2tab -in input_sam -out output_tab -sort yes -ruj yes\n"
            "\e[1mHELP:\e[0m\n"
            "\t-in: input a sam or bam file, the file extension should be .sam or .bam\n"
            "\t-out: output a tab-seperated file \n"
            "\t-sort: <yes/no> whether to sort the output file (default: yes)\n"
            "\t-ruj: <yes/no> whether to remove reads with unannotated junctions, only useful for .sam file (default: yes)\n\n"

            "\e[1mWARNING:\e[0m\n\t%s\n"
            "\e[1mVERSION:\e[0m\n\t%s\n"
            "\e[1mLIB VERSION:\e[0m\n\t%s\n"
            "\e[1mCOMPILE DATE:\e[0m\n\t%s\n"
            "\e[1mAUTHOR:\e[0m\n\t%s\n";

    ostringstream warning;
    warning << YELLOW << WARNING << DEF;

    sprintf(buff, help_info, warning.str().c_str(), BINVERSION, LIBVERSION, DATE, "Li Pan");
    cerr << buff << endl;
}

struct Param
{
    string input_file;
    string output_file;

    bool will_sort = true;
    bool rem_anno_junc = true;

    operator bool()
    { 
        return input_file.empty() or output_file.empty() ? false : true; 
    }

};


void has_next(int argc, int current)
{
    if(current + 1 >= argc)
    {
        cerr << RED << "FATAL ERROR: Parameter Error" << DEF << endl;
        print_usage();
        exit(-1);
    }
}

Param read_param(int argc, char *argv[])
{
    Param param;
    for(size_t i=1; i<argc; i++)
    {
        if( argv[i][0] == '-' )
        {
            if(not strcmp(argv[i]+1, "in"))
            {
                has_next(argc, i);
                param.input_file = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "out"))
            {
                has_next(argc, i);
                param.output_file = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "sort"))
            {
                has_next(argc, i);
                if(not strcmp(argv[i+1], "yes"))
                    param.will_sort = true;
                else
                    param.will_sort = false;
                i++;
            }else if(not strcmp(argv[i]+1, "ruj"))
            {
                has_next(argc, i);
                if(not strcmp(argv[i+1], "yes"))
                    param.rem_anno_junc = true;
                else
                    param.rem_anno_junc = false;
                i++;
            }else{
                cerr << RED << "FATAL ERROR: unknown option: " << argv[i] << DEF << endl;
                print_usage();
                exit(-1);
            }
        }else{
            cerr << RED << "FATAL ERROR: unknown option: " << argv[i] << DEF << endl;
            print_usage();
            exit(-1);
        }
    }
    return param;
}

enum STRAND{ NEG=0, POS=1 };

struct Map_Record
{

    STRAND strand = POS;

    vector<Region> regions;

    Map_Record(const vector<Region> &regions, const STRAND &strand): strand(strand), regions(regions){ }

    Map_Record() = default;
};

struct RecordArray
{
    vector<Map_Record *> content;

    void clear()
    {
        for(Map_Record *p: content)
            delete p;

        content.clear();
    }

    ~RecordArray()
    {
        clear();
    }
};

struct MapChrRecord
{
    MapStringT<RecordArray *> content;

    void clear()
    {
        for(auto it=content.begin(); it!=content.end(); it++)
            delete it->second;

        content.clear();
    }

    ~MapChrRecord()
    {
        clear();
    }
};

/*  sort duplex group by read coordination for PARIS analysis  */
bool Sort_Map_Record(const Map_Record* const mr_1, const Map_Record* const mr_2)
{
    auto &st1 = mr_1->strand;
    auto &st2 = mr_2->strand;
    
    auto &s1 = mr_1->regions.front().first;
    auto &e1 = mr_1->regions.back().second;

    auto &s2 = mr_2->regions.front().first;
    auto &e2 = mr_2->regions.back().second;

    if(st1<st2)
        return true;
    else if(st1>st2)
        return false;
    else{
        if(st1 == POS)
        {
            if(s1<s2)
                return true;
            else if(s1>s2)
                return false;
            else{
                if(e1<e2)
                    return true;
                else
                    return false;
            }
        }else{
            if(e1<e2)
                return false;
            else if(e1>e2)
                return true;
            else{
                if(s1<=s2)
                    return false;
                else
                    return true;
            }
        }
    }
}


void sam2tab(const Param &param)
{
    srand (time(NULL));
    const string randID = to_string(rand());

    BGZF* bam_hd = nullptr;
    bam_hdr_t *hdr = nullptr;
    Sam_Record *p = nullptr;
    ifstream IN;

    if( endswith(param.input_file, ".bam") )
    {
        bam_hd = bgzf_open(param.input_file.c_str(), "r");
        if(not bam_hd)
        {
            cerr << RED << "FATAL Error: open file " << param.input_file << " failed" << DEF << endl;
            exit(-1);
        }
        hdr = bam_hdr_read(bam_hd);
    }else if( endswith(param.input_file, ".sam") )
    {
        IN.open(param.input_file, ifstream::in);
        if(not IN)
        {
            cerr << RED << "FATAL Error: cannot read " << param.input_file << DEF << endl;
            exit(-1);
        }
    }else{
        cerr << RED << "FATAL Error: input file must be sam or bam file -- " << param.input_file << DEF << endl;
        exit(-1);
    }

    cerr << "Start to init " << param.input_file << "\n\t" << currentDateTime() << endl;

    Sam_Head sam_head;
    MapChrRecord mapChrRecords;
    bool success;

    if(bam_hd)
        success = read_sam_head(hdr, sam_head);
    else
        success = read_sam_head(IN, sam_head);

    if(success and sam_head.trans_len.size() >= 1)
    {
        cerr << "\tTotal chr numbers: " << sam_head.trans_len.size() << endl;
        for(auto it=sam_head.trans_len.cbegin(); it!=sam_head.trans_len.cend(); it++)
        {
            const string &chr_id = it->first;
            mapChrRecords.content[chr_id+"+"] = new RecordArray;
            mapChrRecords.content[chr_id+"-"] = new RecordArray;
        }
    }else{
        cerr << RED << "FATAL Error: Read sam head fail..." << DEF << endl;
        exit(-1);
    }

    cerr << "Start to read records... " << endl;

    uLONG removed_unanno = 0;
    uLONG lineCount = 0;
    Sam_Record read_record;

    while(1) 
    {
        if(bam_hd)
            success = read_a_sam_record(bam_hd, hdr, read_record);
        else
            success = read_a_sam_record(IN, read_record);

        if(not success)
            break;

        ++lineCount;
        if(lineCount % 1000000 == 0)
            cerr << "\t lines " << lineCount << endl;

        if(not read_is_mapped(read_record))
            continue;

        if(param.rem_anno_junc)
        {
            bool remove = false;
            for(const string &attr: read_record.attributes)
            {
                const string &title = attr.substr(0, 6);
                if(title == "jM:B:c")
                {
                    const string code_list = attr.substr(7);
                    StringArray codes;
                    split(code_list, ',', codes);
                    long code = stol(codes[0]);
                    if( code >= 0 and code <= 6 )
                    {
                        remove = true;
                        break;
                    }
                }
            }
            if(remove)
            {
                ++removed_unanno;
                continue;
            }
        }

        RegionArray matchRegion;
        get_global_match_region(read_record.cigar, read_record.pos, matchRegion);
        mapChrRecords.content[read_record.chr_id+read_record.strand()]->content.push_back( new Map_Record(matchRegion, read_record.strand() == '+' ? POS : NEG) );
    }

    if(bam_hd)
        bgzf_close(bam_hd);
    else    
        IN.close();

    if(param.rem_anno_junc)
    {
        cerr << YELLOW << "total reads number: " << lineCount << endl;
        cerr << "remove " << removed_unanno << " reads with unannotated junctions" << endl;
    }

    if(param.will_sort)
    {
        cerr << "Start to sort..." << endl;
        for(auto it=mapChrRecords.content.begin(); it!=mapChrRecords.content.end(); it++)
        {
            const string &chr_id = it->first;
            sort(it->second->content.begin(), it->second->content.end(), Sort_Map_Record);
        }
    }

    cerr << "Start to output..." << endl;
    ofstream OUT(param.output_file, ofstream::out);
    if(not OUT)
    {
        cerr << RED << "FATAL Error: cannot write " << param.output_file << DEF << endl;
        exit(-1);
    }

    StringArray chr_list;
    for(auto it=mapChrRecords.content.begin(); it!=mapChrRecords.content.end(); it++)
        chr_list.push_back(it->first);
    std::sort(chr_list.begin(), chr_list.end());
    
    for(auto it=chr_list.begin(); it!=chr_list.end(); it++)
    {
        const string &chr_id_strand = *it;

        const string chr_id = chr_id_strand.substr(0, chr_id_strand.size()-1);
        const char strand = chr_id_strand[chr_id_strand.size()-1];
        
        const auto &records = mapChrRecords.content.at(chr_id_strand);
        for(auto it2=records->content.cbegin(); it2!=records->content.cend(); it2++)
        {
            OUT << chr_id << '\t' << strand;
            for(const Region &r: (*it2)->regions)
                OUT << '\t' << r.first << '\t' << r.second;
            OUT << '\n';
        }
    }
    OUT.close();
}


int main(int argc, char *argv[])
{
    Param param = read_param(argc, argv);
    if(not param)
    {
        print_usage();
        exit(-1);
    }

    sam2tab(param);

    return 0;
}



