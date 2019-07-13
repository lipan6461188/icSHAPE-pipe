
/*

    A Programe to covert sam file into tab-seperated dg file

*/

#include "paris.h"
#include "param.h"
#include "version.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <sstream>
#include "version.h"

using namespace std;
using namespace pan;

#define SAM2DG_VERSION "1.000"
#define DATE __DATE__
#define WARNING "The input of sam2dg must be the output of sam_group or sam_trim"

Color::Modifier RED(Color::FG_RED);
Color::Modifier DEF(Color::FG_DEFAULT);
Color::Modifier YELLOW(Color::FG_YELLOW);


void print_usage()
{
    char buff[2000];
    const char *help_info = 
            "sam2dg - covert sam file to a duplex group tab-seperated file\n"
            "=============================================================\n"
            "\e[1mUSAGE:\e[0m\n"
            "\tsam2dg -in input_sam1,input_sam2... -out output_dg [ -min_overhang 5 \n"
            "\t       -min_armlen 10 -s balance|left -out_sam sam_of_dg -mode pair ]\n"
            "\e[1mHELP:\e[0m\n"
            "\t-min_overhang: mininum overhang of duplex group(default: 5)\n"
            "\t-min_armlen: mininum arm length of each(left/right) arm(default: 10)\n"
            "\t-s: sort balance or left(for paris pipeline) priority (default: no sort)\n"
            "\t-out_sam: output a sam file of dg\n"
            "\t-mode: pair|single. single means each sam record will be a line in dg file (default: pair)\n\n"

            "\e[1mWARNING:\e[0m\n\t%s\n"
            "\e[1mVERSION:\e[0m\n\t%s\n"
            "\e[1mLIB VERSION:\e[0m\n\t%s\n"
            "\e[1mDATE:\e[0m\n\t%s\n"
            "\e[1mAUTHOR:\e[0m\n\t%s\n";
    
    ostringstream warning;
    warning << YELLOW << WARNING << DEF;

    sprintf(buff, help_info, warning.str().c_str(), SAM2DG_VERSION, VERSION, DATE, "Li Pan");
    cout << buff << endl;
}

struct Param
{
    enum SORT_TYPE{ no_sort, balance, left_priority };
    enum MODE{ PAIR_MODE, SINGLE_MODE };

    StringArray input_sam_list;
    string output_dg;
    //string output_unggaped_dg;

    string out_sam_of_dg;

    uINT min_overhang = 5;
    uINT min_armlen = 10;

    SORT_TYPE sort = no_sort;
    MODE mode = PAIR_MODE;

    operator bool(){ return input_sam_list.empty() or output_dg.empty() ? false : true; }
    //bool save_unggapped_dg()const { return not output_unggaped_dg.empty(); }
    bool write_single()const{ return mode == SINGLE_MODE; }
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
                param.input_sam_list = split(argv[i+1], ',');
                i++;
            }else if(not strcmp(argv[i]+1, "out"))
            {
                has_next(argc, i);
                param.output_dg = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "h"))
            {
                print_usage();
                exit(-1);
            }else if(not strcmp(argv[i]+1, "s"))
            {
                has_next(argc, i);
                if(not strcmp(argv[i+1], "balance"))
                    param.sort = Param::balance;
                else if(not strcmp(argv[i+1], "left"))
                    param.sort = Param::left_priority;
                else if(not strcmp(argv[i+1], "no_sort"))
                    param.sort = Param::no_sort;
                else{
                    cerr << RED << "FATAL ERROR: unknown -s option: " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                i++;
            }else if(not strcmp(argv[i]+1, "mode"))
            {
                has_next(argc, i);
                if(not strcmp(argv[i+1], "pair"))
                    param.mode = Param::PAIR_MODE;
                else if(not strcmp(argv[i+1], "single"))
                    param.mode = Param::SINGLE_MODE;
                else{
                    cerr << RED << "FATAL ERROR: unknown -mode option: " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                i++;
            }else if(not strcmp(argv[i]+1, "out_sam"))
            {
                has_next(argc, i);
                param.out_sam_of_dg = argv[i+1];
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


/*  sort duplex group by read coordination for PARIS analysis  */
bool SORT_DG_BY_LEFT_COOR(const Duplex_Hang &dh_1, const Duplex_Hang &dh_2)
{
    // chr 1
    if(dh_1.chr_id_1 < dh_2.chr_id_1)
        return true;
    else if(dh_1.chr_id_1 > dh_2.chr_id_1)
        return false;
    else{
        // strand 1
        if(dh_1.strand_1 < dh_2.strand_1)
            return true;
        else if(dh_1.strand_1 > dh_2.strand_1)
            return false;
        else{
            // chr 2
            if(dh_1.chr_id_2 < dh_2.chr_id_2)
                return true;
            else if(dh_1.chr_id_2 > dh_2.chr_id_2)
                return false;
            else{
                // strand 2
                if(dh_1.strand_2 < dh_2.strand_2)
                    return true;
                else if(dh_1.strand_2 > dh_2.strand_2)
                    return false;
                else{
                    // start 1
                    if(dh_1.start_1 < dh_2.start_1)
                        return true;
                    else if(dh_1.start_1 > dh_2.start_1)
                        return false;
                    else{
                        // end 1
                        if(dh_1.end_1 < dh_2.end_1)
                            return true;
                        else if(dh_1.end_1 > dh_2.end_1)
                            return false;
                        else
                        {
                            // start 2
                            if(dh_1.start_2 < dh_2.start_2)
                                return true;
                            else if(dh_1.start_2 > dh_2.start_2)
                                return false;
                            else
                            {
                                // end 2
                                if(dh_1.end_2 < dh_2.end_2)
                                    return true;
                                else if(dh_1.end_2 > dh_2.end_2)
                                    return false;
                                else{
                                    return false;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

bool valid_dg(vector<Sam_Record> records, const Param &param, Duplex_Hang &dh)//   RegionArray &matchRegion)
{
    if(records.size() == 1)
    {
        RegionArray matchRegion;
        const Sam_Record &record = records.front();
        get_global_match_region(record.cigar, record.pos, matchRegion);
        if(matchRegion.size() != 2)
            return false;
        uLONG gap = matchRegion.at(1).first - matchRegion.at(0).second - 1;
        uLONG left_hang = matchRegion.at(0).second - matchRegion.at(0).first + 1;
        uLONG right_hang = matchRegion.at(1).second - matchRegion.at(1).first + 1;
        if( gap >= param.min_overhang and left_hang >= param.min_armlen and right_hang >= param.min_armlen )
        {
            dh.read_id = records[0].read_id;

            dh.chr_id_1 = records[0].chr_id;
            dh.strand_1 = records[0].strand();
            dh.flag_1 = records[0].flag;
            dh.cigar_1 = records[0].cigar;

            dh.chr_id_2 = records[0].chr_id;
            dh.strand_2 = records[0].strand();
            dh.flag_2 = records[0].flag;
            dh.cigar_2 = records[0].cigar;

            // left - left
            if(matchRegion[0].first < matchRegion[1].first)
            {
                dh.start_1 = matchRegion[0].first;
                dh.end_1 = matchRegion[0].second;
                dh.start_2 = matchRegion[1].first;
                dh.end_2 = matchRegion[1].second;
            }else{
                dh.start_1 = matchRegion[1].first;
                dh.end_1 = matchRegion[1].second;
                dh.start_2 = matchRegion[0].first;
                dh.end_2 = matchRegion[0].second;
            }

            return true;
        }
        else
            return false;
    }else if(records.size() == 2)
    {
        RegionArray matchRegion_1, matchRegion_2;
        get_global_match_region(records[0].cigar, records[0].pos, matchRegion_1);
        get_global_match_region(records[1].cigar, records[1].pos, matchRegion_2);
        if(matchRegion_1.size() != 1 or matchRegion_2.size() != 1)
            return false;
        auto r_1 = matchRegion_1[0];
        auto r_2 = matchRegion_2[0];

        // same chromosome and same strand
        if(records[0].chr_id == records[1].chr_id and records[0].strand() == records[1].strand())
        {
            if(r_1.first > r_2.first)
            {
                swap(records[0], records[1]);
                swap(r_1, r_2);
            }
            uLONG gap = r_2.first - r_1.second - 1;
            uLONG left_hang = r_1.second - r_1.first + 1;
            uLONG right_hang = r_2.second - r_2.first + 1;
            if( gap >= param.min_overhang and left_hang >= param.min_armlen and right_hang >= param.min_armlen )
            {
                dh.read_id = records[0].read_id;

                dh.chr_id_1 = records[0].chr_id;
                dh.strand_1 = records[0].strand();
                dh.flag_1 = records[0].flag;
                dh.cigar_1 = records[0].cigar;
                dh.start_1 = r_1.first;
                dh.end_1 = r_1.second;

                dh.chr_id_2 = records[1].chr_id;
                dh.strand_2 = records[1].strand();
                dh.flag_2 = records[1].flag;
                dh.cigar_2 = records[1].cigar;
                dh.start_2 = r_2.first;
                dh.end_2 = r_2.second;

                return true;
            }
            else
                return false;
        }else{
            // keep left side chr < right side chr
            // or left side strand < right side strand

            dh.read_id = records[0].read_id;

            if(records[0].chr_id < records[1].chr_id or (records[0].chr_id == records[1].chr_id and records[0].strand() < records[1].strand()) )
            {
                dh.chr_id_1 = records[0].chr_id;
                dh.strand_1 = records[0].strand();
                dh.flag_1 = records[0].flag;
                dh.cigar_1 = records[0].cigar;
                dh.start_1 = r_1.first;
                dh.end_1 = r_1.second;

                dh.chr_id_2 = records[1].chr_id;
                dh.strand_2 = records[1].strand();
                dh.flag_2 = records[1].flag;
                dh.cigar_2 = records[1].cigar;
                dh.start_2 = r_2.first;
                dh.end_2 = r_2.second;
            }else
            {
                dh.chr_id_2 = records[0].chr_id;
                dh.strand_2 = records[0].strand();
                dh.flag_2 = records[0].flag;
                dh.cigar_2 = records[0].cigar;
                dh.start_2 = r_1.first;
                dh.end_2 = r_1.second;

                dh.chr_id_1 = records[1].chr_id;
                dh.strand_1 = records[1].strand();
                dh.flag_1 = records[1].flag;
                dh.cigar_1 = records[1].cigar;
                dh.start_1 = r_2.first;
                dh.end_1 = r_2.second;
            }
            return true;
        }
    }else{
        return false;
    }
}

void to_dh(const vector<Sam_Record> &records, const Param &param, vector<Duplex_Hang> &dh_array)//   RegionArray &matchRegion)
{
    dh_array.clear();

    //cout << "to_dh" << endl;

    for(const Sam_Record &record: records)
    {
        RegionArray matchRegion;
        get_global_match_region(record.cigar, record.pos, matchRegion);
        if(matchRegion.size() == 1)
        {
            Duplex_Hang dh;

            dh.read_id = record.read_id;

            dh.chr_id_1 = record.chr_id;
            dh.strand_1 = record.strand();
            dh.flag_1 = record.flag;
            dh.cigar_1 = record.cigar;
            dh.start_1 = matchRegion[0].first;
            dh.end_1 = matchRegion[0].second;

            dh.chr_id_2 = record.chr_id;
            dh.strand_2 = record.strand();
            dh.flag_2 = record.flag;
            dh.cigar_2 = record.cigar;
            dh.start_2 = matchRegion[0].first;
            dh.end_2 = matchRegion[0].second;

            dh_array.push_back(dh);
        
            //cout << "to_dh 1" << endl;

        }else if(matchRegion.size() == 2){
            
            //uLONG gap = matchRegion.at(1).first - matchRegion.at(0).second - 1;
            //uLONG left_hang = matchRegion.at(0).second - matchRegion.at(0).first + 1;
            //uLONG right_hang = matchRegion.at(1).second - matchRegion.at(1).first + 1;

            Duplex_Hang dh;

            dh.read_id = record.read_id;

            dh.chr_id_1 = record.chr_id;
            dh.strand_1 = record.strand();
            dh.flag_1 = record.flag;
            dh.cigar_1 = record.cigar;

            dh.chr_id_2 = record.chr_id;
            dh.strand_2 = record.strand();
            dh.flag_2 = record.flag;
            dh.cigar_2 = record.cigar;

            // left - left
            if(matchRegion[0].first < matchRegion[1].first)
            {
                dh.start_1 = matchRegion[0].first;
                dh.end_1 = matchRegion[0].second;
                dh.start_2 = matchRegion[1].first;
                dh.end_2 = matchRegion[1].second;
            }else{
                dh.start_1 = matchRegion[1].first;
                dh.end_1 = matchRegion[1].second;
                dh.start_2 = matchRegion[0].first;
                dh.end_2 = matchRegion[0].second;
            }

            dh_array.push_back(dh);

            //cout << "to_dh 2" << endl;

        }else{
            // no thing

            //cout << "to_dh 3 " << record.cigar << endl;
        }
    }
}

void sam2dg(const Param &param)
{
    bool out_sam = (not param.out_sam_of_dg.empty());

    ifstream IN;
    ofstream DG(param.output_dg, ofstream::out);
    ofstream SAM;

    /*
    if(param.save_unggapped_dg())
    {
        DG_UNGAPED.open(param.output_unggaped_dg, ofstream::out);
        if(not DG_UNGAPED)
        {
            cerr << RED << "FATAL Error: " << param.output_unggaped_dg << " is unreadable" << DEF << endl;
            exit(-1);
        }
    }
    */

    if(not DG)
        throw runtime_error("File "+param.output_dg+" cannot be writeable");
    if(out_sam)
    {
        SAM.open(param.out_sam_of_dg, ofstream::out);
        if(not DG)
            throw runtime_error("File "+param.out_sam_of_dg+" cannot be writeable");
    }

    bool success;
    vector<Sam_Record> read_records;
    Sam_Head sam_head;
    unordered_multimap<string, Sam_Record> records;
    vector<Duplex_Hang> duplex_hangs;

    //MapStringT<uLONGLONG> read_counter;

    for(size_t idx=0; idx<param.input_sam_list.size(); idx++)
    {
        IN.open(param.input_sam_list[idx]);
        if(not IN)
            throw runtime_error("File "+param.input_sam_list[idx]+" cannot be readable");

        if(read_sam_head(IN, sam_head) and out_sam and idx == 0)
            write_sam_head(SAM, sam_head);

        while(not IN.eof())
        {
            success = read_a_read_record(IN, read_records);
            if(not success)
                break;

            Duplex_Hang dh;
            vector<Duplex_Hang> dh_array;

            bool is_valid(true);
            if(param.write_single())
            {
                to_dh(read_records, param, dh_array);
                //cout << dh_array.size() << endl;
            }
            else
                is_valid = valid_dg(read_records, param, dh);

            if(is_valid)
            {
                if(param.sort != Param::no_sort)
                    // sort
                {
                    if(out_sam)
                    {
                        for_each(read_records.cbegin(), read_records.cend(), [&records](const Sam_Record& record){ records.insert({record.read_id, record}); });
                    }
                    if(param.write_single())
                        duplex_hangs.insert(duplex_hangs.end(), dh_array.cbegin(), dh_array.cend());
                        //duplex_hangs += dh_array;
                    else
                        duplex_hangs.push_back(dh);
                }else{
                    // do not sort
                    if(out_sam)
                        SAM << read_records;

                    if(param.write_single())
                    {
                        DG << dh_array;
                        //cout << dh_array << endl;
                    }
                    else
                        DG << dh;
                }
            }
        }

        IN.close();
    }

    if(param.sort != Param::no_sort)
    {
        if(param.sort == Param::balance)
            sort(duplex_hangs.begin(), duplex_hangs.end());
        else if(param.sort == Param::left_priority)
            sort(duplex_hangs.begin(), duplex_hangs.end(), SORT_DG_BY_LEFT_COOR);

        DG << duplex_hangs;
        if(out_sam)
        {
            for(auto iter=duplex_hangs.cbegin(); iter!=duplex_hangs.cend(); iter++)
            {
                auto range = records.equal_range( iter->read_id );
                for(auto begin=range.first; begin!=range.second; begin++)
                    SAM << begin->second;
            }
        }
    }
    DG.close();
    SAM.close();
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

    sam2dg(param);

    return 0;
}



