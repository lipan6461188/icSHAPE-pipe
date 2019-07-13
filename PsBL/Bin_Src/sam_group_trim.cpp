/*

    A Programe to filter sam file for paris pipeline

*/

#include "sam.h"
#include "param.h"
#include "version.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <sstream>

using namespace std;
using namespace pan;

#define SAM_GROUP_TRIM_VERSION "1.000"
#define DATE __DATE__

Color::Modifier RED(Color::FG_RED);
Color::Modifier DEF(Color::FG_DEFAULT);
Color::Modifier YELLOW(Color::FG_YELLOW);
Color::Modifier CYAN(Color::FG_LIGHT_CYAN);

void print_usage()
{
    char buff[2000];
    const char *help_info = 
            "sam_group_trim - trim or group sam file\n"
            "=======================================\n"
            "\e[1mUSAGE:\e[0m\n"
            "\tsam_group_trim -in input_sam1,input_sam2... -out output_sam [ -out_ungapped output_ungapped_sam  \n"
            "\t                     -max_ungapped_sc 5 -min_overhang 5 -min_armlen 15 -max_ambiguous_base 10 \n"
            "\t                     -remove_antisense no -only_primary no ]    \n"
            "\e[1mHELP:\e[0m\n"
            "\t-in: input sam files, the output sam file head is the head of the first sam file\n"
            "\t-out: output sam file\n"
            "\t-out_ungapped: output all ungapped reads to a sam file\n\n"

            "\t-max_ungapped_sc: max soft clipped bases are allowed in a ungapped read (default: 5)\n"
            "\t-min_overhang: mininum overhang of duplex group (default: 5)\n"
            "\t-min_armlen: mininum arm length of each(left/right) arm (default: 15)\n"
            "\t-max_ambiguous_base: maximum ambiguous base numbers (default: 10) \n"
            "\t-remove_antisense: yes/no. remove any reads map to antisense (default: no)\n"
            "\t-only_primary: yes/no. only preserve primary map for multiple-mapped reads (default: no)\n\n"

            "\e[1mVERSION:\e[0m\n\t%s\n"
            "\e[1mLIB VERSION:\e[0m\n\t%s\n"
            "\e[1mDATE:\e[0m\n\t%s\n"
            "\e[1mAUTHOR:\e[0m\n\t%s\n";

    //ostringstream warning;
    //warning << YELLOW << WARNING << DEF;

    sprintf(buff, help_info, SAM_GROUP_TRIM_VERSION, VERSION, DATE, "Li Pan");
    cout << buff << endl;
}

struct Param
{
    StringArray input_sam_list;
    string output_sam;
    string output_ungapped_sam;

    uINT max_ungapped_sc = 5;
    uINT min_overhang = 5;
    uINT min_armlen = 15;
    uINT max_ambiguous_base = 10;
    bool remove_antisense = false;
    bool only_primary = false;

    operator bool(){ return input_sam_list.empty() or output_sam.empty() ? false : true; }
    bool save_unggapped_reads() const { return not output_ungapped_sam.empty(); }
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
            }else if(not strcmp(argv[i]+1, "max_ungapped_sc"))
            {
                has_next(argc, i);
                param.max_ungapped_sc = stoul(string(argv[i+1]));
                i++;
            }else if(not strcmp(argv[i]+1, "in"))
            {
                has_next(argc, i);
                param.input_sam_list = split(argv[i+1], ',');
                i++;
            }else if(not strcmp(argv[i]+1, "out"))
            {
                has_next(argc, i);
                param.output_sam = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "out_ungapped"))
            {
                has_next(argc, i);
                param.output_ungapped_sam = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "remove_antisense"))
            {
                has_next(argc, i);
                if( not strcmp(argv[i+1], "yes") )
                {
                    param.remove_antisense = true;
                }else if( not strcmp(argv[i+1], "no") )
                {
                    param.remove_antisense = false;
                }else{
                    cerr << RED << "FATAL Error: unknown option of --remove_antisense " << DEF << endl;
                    exit(-1);
                }
                i++;
            }else if(not strcmp(argv[i]+1, "only_primary"))
            {
                has_next(argc, i);
                if( not strcmp(argv[i+1], "yes") )
                {
                    param.only_primary = true;
                }else if( not strcmp(argv[i+1], "no") )
                {
                    param.only_primary = false;
                }else{
                    cerr << RED << "FATAL Error: unknown option of --only_primary " << DEF << endl;
                    exit(-1);
                }
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

sp<Sam_Record> group_reads( const ReadPair &read_pair, 
                            const uLONG min_overhang, 
                            const uLONG min_armlen, 
                            const uLONG max_ambiguous_base)
{
    sp<Sam_Record> sp_grouped_record = make_shared<Sam_Record>();
    if(group_read_pair(read_pair, *sp_grouped_record) == GROUP_SUCCESS)
    {
        RegionArray matchRegion;
        get_global_match_region((*sp_grouped_record).cigar, (*sp_grouped_record).pos, matchRegion);

        uLONG gap = matchRegion.at(1).first - matchRegion.at(0).second - 1;
        uLONG left_hang = matchRegion.at(0).second - matchRegion.at(0).first + 1;
        uLONG right_hang = matchRegion.at(1).second - matchRegion.at(1).first + 1;

        if( gap < min_overhang or left_hang < min_armlen or right_hang < min_armlen )
            return nullptr;
        if( stoul((*sp_grouped_record).get_attr("AB")) > max_ambiguous_base )
            return nullptr;

        return sp_grouped_record;
        // process read flag, read id and other....
    }
    return nullptr;
}

sp<ReadPair> trim_reads(const ReadPair &read_pair,
                        const uLONG min_armlen, 
                        const uLONG max_ambiguous_base)
{
    sp<ReadPair> sp_trimmed_record = make_shared<ReadPair>();
    if(trim_read_pair(read_pair, *sp_trimmed_record) == TRIM_SUCCESS)
    {
        RegionArray matchRegion_1, matchRegion_2;
        get_global_match_region((*sp_trimmed_record).first.cigar, (*sp_trimmed_record).first.pos, matchRegion_1);
        get_global_match_region((*sp_trimmed_record).second.cigar, (*sp_trimmed_record).second.pos, matchRegion_2);

        uLONG hang_1 = matchRegion_1.at(0).second - matchRegion_1.at(0).first + 1;
        uLONG hang_2 = matchRegion_2.at(0).second - matchRegion_2.at(0).first + 1;

        if( hang_1 < min_armlen or hang_2 < min_armlen )
            return nullptr;
        if( stoul((*sp_trimmed_record).first.get_attr("AB")) > max_ambiguous_base )
            return nullptr;

        return sp_trimmed_record;
        // process read flag, read id and other....
    }
    return nullptr;
}

struct Return_Type
{
    bool single_record = false;
    bool read_pair = false;

    sp<Sam_Record> sp_record;
    sp<ReadPair> sp_pair;
};

void trim_or_group( vector<Sam_Record> &records,
                    vector<Return_Type> &return_list,
                    const uLONG min_overhang, 
                    const uLONG min_armlen, 
                    const uLONG max_ambiguous_base,
                    const bool remove_antisense,
                    const bool only_one )
{
    //vector<Return_Type> return_list;
    return_list.clear();

    // filter reads
    for(auto iter=records.begin(); iter != records.end(); )
    {
        RegionArray matchRegion;
        get_global_match_region(iter->cigar, iter->pos, matchRegion);
        if(matchRegion.size() == 2)
        {
            uLONG gap = matchRegion.at(1).first - matchRegion.at(0).second - 1;
            uLONG hang_1 = matchRegion.at(0).second - matchRegion.at(0).first + 1;
            uLONG hang_2 = matchRegion.at(1).second - matchRegion.at(1).first + 1;

            if(remove_antisense)
            {
                bool read_is_antisense = read_is_reverse( *iter );
                if(read_is_antisense)
                {
                    iter = records.erase(iter);
                    continue;
                }
            }

            if( hang_1 < min_armlen or hang_2 < min_armlen or gap < min_overhang )
            {
                iter++;
                continue;
            }
            else{
                Return_Type rt;
                rt.single_record = true;
                rt.sp_record = make_shared<Sam_Record>(*iter);
                rt.sp_record->attributes.clear();
                rt.sp_record->attributes.push_back("AB:i:0");
                return_list.push_back(rt);
                iter = records.erase(iter);

                if(only_one)
                    return;
            }
        }else{
            iter++;
            continue;
        }
    }

    // search valid reads
    for(size_t idx=0; idx<records.size(); idx++)
        for(size_t idy=idx+1; idy<records.size(); idy++)
        {
            bool same_chr = (records[idx].chr_id == records[idy].chr_id);
            bool same_strand = ( records[idx].strand() == records[idy].strand() );
            
            if(remove_antisense)
            {
                bool read_is_antisense = ( read_is_reverse(records[idx]) or read_is_reverse(records[idy]) );
                if(read_is_antisense)
                    continue;
            }

            ReadPair read_pair(records[idx], records[idy]);

            if(same_chr and same_strand)
            {
                // group
                sp<Sam_Record> sp_record = group_reads( read_pair, min_overhang, min_armlen, max_ambiguous_base);
                if(sp_record)
                {
                    Return_Type rt;
                    rt.single_record = true;
                    rt.sp_record = sp_record;
                    return_list.push_back(rt);

                    if(only_one)
                        break;
                }
            }else{
                // trim
                sp<ReadPair> sp_pair = trim_reads(read_pair, min_armlen, max_ambiguous_base);
                if(sp_pair)
                {
                    Return_Type rt;
                    rt.read_pair = true;
                    rt.sp_pair = sp_pair;
                    return_list.push_back(rt);

                    if(only_one)
                        break;
                }
            }
        }
}

void sam_group_trim(const Param &param)
{
    using size_type = vector<Sam_Record>::size_type;

    ifstream IN;
    ofstream OUT(param.output_sam, ofstream::out);
    ofstream OUT_UNGAPPED;
    if(param.save_unggapped_reads())
    {
        OUT_UNGAPPED.open(param.output_ungapped_sam, ofstream::out);
        if(not OUT_UNGAPPED)
        {
            cerr << "FATAL Error: " << param.output_ungapped_sam << " is unwritable" << endl;
            exit(-1);
        }
    }

    if(not OUT)
    {
        cerr << "FATAL Error: " << param.output_sam << " is unwritable" << endl;
        exit(-1);
    }

    Sam_Head sam_head;
    vector<Sam_Record> read_records;
    vector<Return_Type> return_list;

    for(size_type idx=0; idx<param.input_sam_list.size(); idx++)
    {
        IN.open(param.input_sam_list[idx]);
        if(not IN)
            throw runtime_error("File "+param.input_sam_list[idx]+" cannot be readable");

        if(read_sam_head(IN, sam_head) and idx == 0)
        {
            write_sam_head(OUT, sam_head);
            if(param.save_unggapped_reads())
                write_sam_head(OUT_UNGAPPED, sam_head);
        }

        clog << "Start to process " << param.input_sam_list[idx] << "...\n";
        
        uLONGLONG count_index = 0;
        uLONGLONG valid_reads = 0;

        while(not IN.eof())
        {
            count_index++;
            if(count_index % 1000000 == 0)
                clog << "\t processed " << count_index << " reads (" << param.input_sam_list[idx] << ")..." << endl;
            //bool success = read_a_read_record(IN, read_records);
            if(not read_a_read_record(IN, read_records))
                break;

            if(param.save_unggapped_reads())
            {
                for(const Sam_Record &record: read_records)
                {
                    if( not read_is_gapped(record) and read_max_softclip(record.cigar) <= param.max_ungapped_sc )
                    {
                        if(param.remove_antisense and read_is_reverse(record))
                            continue;
                        OUT_UNGAPPED << record;
                    }
                }
            }

            trim_or_group( read_records, return_list, param.min_overhang, param.min_armlen, param.max_ambiguous_base, param.remove_antisense, param.only_primary );
            
            if(return_list.size() != 0)
                ++valid_reads;

            uLONG mapped_index = 0;
            for(const Return_Type &return_item: return_list)
            {
                if(return_item.single_record)
                {
                    if(return_list.size() != 1)
                    {
                        return_item.sp_record->read_id += "_" + to_string(mapped_index);
                        mapped_index++;
                    }
                    remove_read_mutimap( *return_item.sp_record );
                    return_item.sp_record->map_quanlity = 1;
                    OUT << *return_item.sp_record;
                }

                else if(return_item.read_pair)
                {
                    if(return_list.size() != 1)
                    {
                        return_item.sp_pair->first.read_id += "_" + to_string(mapped_index);
                        return_item.sp_pair->second.read_id += "_" + to_string(mapped_index);
                        mapped_index++;
                    }
                    remove_read_mutimap( return_item.sp_pair->first );
                    add_read_mutimap( return_item.sp_pair->second );
                    return_item.sp_pair->first.map_quanlity = 1;
                    return_item.sp_pair->first.map_quanlity = 1;
                    OUT << *return_item.sp_pair;
                }else{
                    cerr << RED << "Unexpected Error" << DEF << endl;
                }
            }
        }
        IN.close();
        if(param.save_unggapped_reads())
            OUT_UNGAPPED.close();

        clog << CYAN << "Summary: " << param.input_sam_list[idx] << ": \n\t" 
                     << "Read Number: " << count_index << "\n\t" 
                     << "Valid Read Number: " << valid_reads << "\n\t"
                     << "Valid Read Ratio: " << 1.0*valid_reads/count_index * 100 << "%" << DEF << endl;
    }
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

    sam_group_trim(param);

    return 0;
}


    



