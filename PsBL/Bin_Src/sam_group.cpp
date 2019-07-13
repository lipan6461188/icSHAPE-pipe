/*

    A Programe to group read_id sorted sam file
    group muti-mapped ungapped reads into single gapped reads

*/


#include "paris.h"
#include "param.h"
#include "version.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdexcept>

using namespace std;
using namespace pan;

#define SAM2GROUP_VERSION "1.000"
#define DATE __DATE__

Color::Modifier RED(Color::FG_RED);
Color::Modifier DEF(Color::FG_DEFAULT);

void print_usage()
{
    char buff[2000];
    const char *help_info = 
            "sam_group - group muti-mapped ungapped reads into single gapped reads\n"
            "=====================================================================\n"
            "\e[1mUSAGE:\e[0m\n"
            "\tsam_group -in input_sam -out output_sam [-min_overhang 5 -min_armlen 10 -max_ambiguous_base 10]\n"
            "\e[1mHELP:\e[0m\n"
            "\t-min_overhang: mininum overhang of duplex group(default: 5)\n"
            "\t-min_armlen: mininum arm length of each(left/right) arm(default: 10)\n"
            "\t-max_ambiguous_base: maximum ambiguous base numbers(default: 10) \n"
            "\e[1mVERSION:\e[0m\n\t%s\n"
            "\e[1mLIB VERSION:\e[0m\n\t%s\n"
            "\e[1mDATE:\e[0m\n\t%s\n"
            "\e[1mAUTHOR:\e[0m\n\t%s\n";
    sprintf(buff, help_info, SAM2GROUP_VERSION, VERSION, DATE, "Li Pan");
    cout << buff << endl;
}

struct Param
{
    string input_sam;
    string output_sam;

    uINT min_overhang = 5;
    uINT min_armlen = 10;
    uINT max_ambiguous_base = 10;

    operator bool(){ return input_sam.empty() or output_sam.empty() ? false : true; }
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
            if(not strcmp(argv[i]+1, "min_overhang"))
            {
                has_next(argc, i);
                param.min_overhang = stoul(string(argv[i+1]));
                i++;
            }else if(not strcmp(argv[i]+1, "min_armlen"))
            {
                has_next(argc, i);
                param.min_armlen = stoul(string(argv[i+1]));
                i++;
            }else if(not strcmp(argv[i]+1, "in"))
            {
                has_next(argc, i);
                param.input_sam = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "out"))
            {
                has_next(argc, i);
                param.output_sam = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "h"))
            {
                print_usage();
                exit(-1);
            }else if(not strcmp(argv[i]+1, "max_ambiguous_base"))
            {
                has_next(argc, i);
                param.max_ambiguous_base = stoul(string(argv[i+1]));
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

void group_reads(const Param &param)
{
    using size_type = vector<Sam_Record>::size_type;

    ifstream IN(param.input_sam, ifstream::in);
    ofstream OUT(param.output_sam, ofstream::out);

    uLONGLONG line_count = 0;

    Sam_Head sam_head;
    read_sam_head(IN, sam_head);
    OUT << sam_head;

    vector<Sam_Record> read_records;
    while(not IN.eof())
    {
        bool success = read_a_read_record(IN, read_records);
        if(not success)
        {
            break;
         //   cerr << "Fail..." << endl;
           // exit(-1);
        }
        if(read_records.size() >= 2)
        {
            //bool first = true;
            int index = 1;
            vector<Sam_Record> read_array;
            for(size_type idx=0; idx<read_records.size(); idx++)
                for(size_type idy=idx+1; idy<read_records.size(); idy++)
                {
                    ReadPair read_pair(read_records.at(idx), read_records.at(idy));
                    Sam_Record groupped_read;
                    if(group_read_pair(read_pair, groupped_read) == GROUP_SUCCESS)
                    {
                        RegionArray matchRegion;
                        get_global_match_region(groupped_read.cigar, groupped_read.pos, matchRegion);

                        uLONG gap = matchRegion.at(1).first - matchRegion.at(0).second - 1;
                        uLONG left_hang = matchRegion.at(0).second - matchRegion.at(0).first + 1;
                        uLONG right_hang = matchRegion.at(1).second - matchRegion.at(1).first + 1;

                        if( gap < param.min_overhang or left_hang < param.min_armlen or right_hang < param.min_armlen )
                            continue;
                        if( stoul(groupped_read.get_attr("AB")) > param.max_ambiguous_base )
                            continue;

                        groupped_read.flag &= ~256UL;

                        /*
                        if(first)
                        {
                            first = false;
                            groupped_read.flag &= ~256UL;
                        }else{
                            groupped_read.flag |= 256UL;
                        }
                        */

                        //groupped_read.attributes.push_back("MI:i:"+to_string(index));
                        groupped_read.read_id += "_" + to_string(index);
                        read_array.push_back(groupped_read);
                        index++;
                    }
                }
            for(auto iter=read_array.begin(); iter!=read_array.end(); iter++)
            {
                iter->attributes.push_back("NH:i:"+to_string(read_array.size()));
                OUT << *iter;
            }
        }
        ++line_count;
        if(line_count % 1000000 == 0)
            clog << "Read " << line_count << " lines...\n";
    }
    IN.close();
    OUT.close();
}

int main(int argc, char *argv[])
{
    Param param = read_param(argc, argv);
    if(not param)
    {
        cerr << RED << "Parameter is not valid" << DEF << endl;
        print_usage();
        exit(-1);
    }

    group_reads(param);

    return 0;
}




