
/*

    A Programe to trim read_id sorted sam file
    trim muti-mapped ungapped reads to read pairs

*/


#include "paris.h"
#include "param.h"
#include "version.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <thread>
#include <future>


using namespace std;
using namespace pan;

#define SAMTRIM_VERSION "1.000"
#define DATE __DATE__

Color::Modifier RED(Color::FG_RED);
Color::Modifier DEF(Color::FG_DEFAULT);

void print_usage()
{
    char buff[2000];
    const char *help_info = 
            "sam_trim - trim muti-mapped ungapped reads to read pairs\n"
            "=====================================================================\n"
            "\e[1mUSAGE:\e[0m\n"
            "\tsam_trim -in input_sam -out output_sam [-min_armlen 10 -max_ambiguous_base 10\n"
            "\t         -max_duplex 1000 -remove_antisense -remove_groupable_reads -threads 1\n"
            "\t         -threads_load 5000                                                  ]\n"
            "\e[1mHELP:\e[0m\n"
            "\t-min_armlen: mininum arm length of each(left/right) arm(default: 10)\n"
            "\t-max_ambiguous_base: maximum ambiguous base numbers(default: 10) \n"
            "\t-max_duplex: maximum duplex number(default: 1000) \n"
            "\t-remove_antisense: remove any reads map to antisense(default: no) \n"
            "\t-remove_groupable_reads: remove any read pairs which can be groupped(default: no) \n"
            "\t-threads: threads number(default: 1) \n"
            "\t-threads_load: how many reads to process for each thread(default: 5000) \n"
            "\e[1mVERSION:\e[0m\n\t%s\n"
            "\e[1mLIB VERSION:\e[0m\n\t%s\n"
            "\e[1mDATE:\e[0m\n\t%s\n"
            "\e[1mAUTHOR:\e[0m\n\t%s\n";
    sprintf(buff, help_info, SAMTRIM_VERSION, VERSION, DATE, "Li Pan");
    cout << buff << endl;
}

struct Param
{
    string input_sam;
    string output_sam;

    uINT min_armlen = 10;
    uINT max_ambiguous_base = 10;
    uINT max_duplex = 1000;

    uINT threads = 1;
    uINT threads_load = 5000;

    bool remove_antisense = false;
    bool remove_groupable_reads = false;

    operator bool(){ return input_sam.empty() or output_sam.empty() or threads == 0 or threads_load == 0 ? false : true; }
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
            if(not strcmp(argv[i]+1, "min_armlen"))
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
            }else if(not strcmp(argv[i]+1, "max_duplex"))
            {
                has_next(argc, i);
                param.max_duplex = stoul(string(argv[i+1]));
                i++;
            }else if(not strcmp(argv[i]+1, "threads"))
            {
                has_next(argc, i);
                param.threads = stoul(string(argv[i+1]));
                i++;
            }else if(not strcmp(argv[i]+1, "threads_load"))
            {
                has_next(argc, i);
                param.threads_load = stoul(string(argv[i+1]));
                i++;
            }else if(not strcmp(argv[i]+1, "remove_antisense"))
            {
                param.remove_antisense = true;
            }else if(not strcmp(argv[i]+1, "remove_groupable_reads"))
            {
                param.remove_groupable_reads = true;
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


void trim_reads_single(const Param &param)
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
        ++line_count;
        if(line_count % 1000000 == 0)
            clog << "Read " << line_count << " lines...\n";

        bool success = read_a_read_record(IN, read_records);
        if(not success)
        {
            break;
        }

        if(read_records.size() >= 2)
        {
            int index = 0;
            vector<ReadPair> readp_array;
            
            //#pragma omp parallel for
            for(size_type idx=0; idx<read_records.size(); idx++)
            {
                bool idx_reverse_read = read_is_reverse(read_records.at(idx));
                if(param.remove_antisense and idx_reverse_read)
                    continue;
                for(size_type idy=idx+1; idy<read_records.size(); idy++)
                {
                    bool idy_reverse_read = read_is_reverse(read_records.at(idy));
                    if(param.remove_antisense and idy_reverse_read)
                        continue;

                    ReadPair read_pair(read_records.at(idx), read_records.at(idy));
                    ReadPair trimmed_read_pair;

                    bool groupable = (read_pair.first.chr_id == read_pair.second.chr_id) and ( not (idx_reverse_read ^ idy_reverse_read) );
                    if(param.remove_groupable_reads and groupable)
                        continue;

                    TRIM_STATUS success;
                    try{
                        //clog << "call trim_read_pair..." << endl;
                        success = trim_read_pair(read_pair, trimmed_read_pair);
                    }catch(std::out_of_range e)
                    {
                        cerr << e.what() << endl;
                        cerr << read_pair << endl;
                        continue;
                    }

                    if(success == TRIM_SUCCESS)
                    {
                        RegionArray matchRegion_1, matchRegion_2;
                        get_global_match_region(trimmed_read_pair.first.cigar, trimmed_read_pair.first.pos, matchRegion_1);
                        get_global_match_region(trimmed_read_pair.second.cigar, trimmed_read_pair.second.pos, matchRegion_2);

                        uLONG hang_1 = matchRegion_1.at(0).second - matchRegion_1.at(0).first + 1;
                        uLONG hang_2 = matchRegion_2.at(0).second - matchRegion_2.at(0).first + 1;

                        if( hang_1 < param.min_armlen or hang_2 < param.min_armlen )
                            continue;
                        if( stoul(trimmed_read_pair.first.get_attr("AB")) > param.max_ambiguous_base )
                            continue;

                        trimmed_read_pair.first.read_id += "_" + to_string(index);
                        trimmed_read_pair.second.read_id += "_" + to_string(index);
                        readp_array.push_back( trimmed_read_pair );
                        index++;
                    }
                }
            }
            if(readp_array.size() > param.max_duplex)
                continue;
            for(auto iter=readp_array.begin(); iter!=readp_array.end(); iter++)
            {
                iter->first.attributes.push_back("NH:i:"+to_string(readp_array.size()));
                iter->second.attributes.push_back("NH:i:"+to_string(readp_array.size()));
                iter->first.flag &= ~256UL;
                iter->second.flag |= 256UL;
                OUT << *iter;

            }
        }
    }
    IN.close();
    OUT.close();
}



bool trim_single_read(vector<Sam_Record> * const p_read_records, const Param * const param, vector<ReadPair> *readp_array)
{
    int index = 0;
    //readp_array = new vector<ReadPair>();

    //clog << "1..." << endl;

    for(uINT idx=0; idx<p_read_records->size(); idx++)
    {
        bool idx_reverse_read = read_is_reverse(p_read_records->at(idx));
        if(param->remove_antisense and idx_reverse_read)
            continue;
        for(uINT idy=idx+1; idy<p_read_records->size(); idy++)
        {
            bool idy_reverse_read = read_is_reverse(p_read_records->at(idy));
            if(param->remove_antisense and idy_reverse_read)
                continue;

            ReadPair read_pair(p_read_records->at(idx), p_read_records->at(idy));
            ReadPair trimmed_read_pair;

            bool groupable = (read_pair.first.chr_id == read_pair.second.chr_id) and ( not (idx_reverse_read ^ idy_reverse_read) );
            if(param->remove_groupable_reads and groupable)
                continue;

            TRIM_STATUS success;
            try{
                //clog << "call trim_read_pair..." << endl;
                success = trim_read_pair(read_pair, trimmed_read_pair);
            }catch(std::out_of_range e)
            {
                cerr << e.what() << endl;
                cerr << read_pair << endl;
                continue;
            }

            if(success == TRIM_SUCCESS)
            {
                RegionArray matchRegion_1, matchRegion_2;
                get_global_match_region(trimmed_read_pair.first.cigar, trimmed_read_pair.first.pos, matchRegion_1);
                get_global_match_region(trimmed_read_pair.second.cigar, trimmed_read_pair.second.pos, matchRegion_2);

                uLONG hang_1 = matchRegion_1.at(0).second - matchRegion_1.at(0).first + 1;
                uLONG hang_2 = matchRegion_2.at(0).second - matchRegion_2.at(0).first + 1;

                if( hang_1 < param->min_armlen or hang_2 < param->min_armlen )
                    continue;
                if( stoul(trimmed_read_pair.first.get_attr("AB")) > param->max_ambiguous_base )
                    continue;

                trimmed_read_pair.first.read_id += "_" + to_string(index);
                trimmed_read_pair.second.read_id += "_" + to_string(index);
                readp_array->push_back( trimmed_read_pair );
                index++;
            }
        }
    }
    delete p_read_records;

    //clog << "2..." << endl;

    if(readp_array->size() > param->max_duplex)
    {
        //delete readp_array;
        return false;
    }
    for(auto iter=readp_array->begin(); iter!=readp_array->end(); iter++)
    {
        iter->first.attributes.push_back("NH:i:"+to_string(readp_array->size()));
        iter->second.attributes.push_back("NH:i:"+to_string(readp_array->size()));
        iter->first.flag &= ~256UL;
        iter->second.flag |= 256UL;
        //OUT << *iter;
    }
    //clog << "3..." << endl;
    return true;
}

using MULTI_TRIM_TYPE = vector< pair<vector<ReadPair> *, bool> >;

MULTI_TRIM_TYPE trim_multiple_reads(vector< vector<Sam_Record>* >* p_reads_factory, const Param * const param)
{
    MULTI_TRIM_TYPE my_list;

    for(vector<Sam_Record>* pointer: *p_reads_factory)
    {
        vector<ReadPair> *readp_array = new vector<ReadPair>();
        bool success(true);
        success = trim_single_read(pointer, param, readp_array);
        my_list.push_back( make_pair(readp_array,  success) );
    }
    delete p_reads_factory;
    return my_list;
}

void trim_reads_parallel(const Param &param)
{
    //using size_type = vector<Sam_Record>::size_type;

    ifstream IN(param.input_sam, ifstream::in);
    ofstream OUT(param.output_sam, ofstream::out);

    uLONGLONG line_count = 0;

    Sam_Head sam_head;
    read_sam_head(IN, sam_head);
    OUT << sam_head;

    bool finish(false);
    while(not IN.eof())
    {
        uINT threads_number = 0;
        vector< future<MULTI_TRIM_TYPE> > worker;
        //vector< vector<ReadPair>* > readp_arrays(param.threads);

        while(threads_number < param.threads)
        {
            vector< vector<Sam_Record>* >* p_reads_factory = new vector< vector<Sam_Record>* >();
            uINT reads_number = 0;
            while(reads_number < param.threads_load)
            {
                ++line_count;
                if(line_count % 1000000 == 0)
                    clog << "Read " << line_count << " reads...\n";

                vector<Sam_Record> *sp_read_records = new vector<Sam_Record>();
                bool success = read_a_read_record(IN, *sp_read_records);
                if(not success)
                {
                    finish = true;
                    break;
                }
                //cout << sp_read_records->size() << endl;
                if(sp_read_records->size() < 2)
                {
                    delete sp_read_records;
                    continue;
                }
                p_reads_factory->push_back(sp_read_records);
                ++reads_number;
            }
            if(reads_number > 0)
            {
                worker.emplace_back(  std::async(trim_multiple_reads, p_reads_factory, &param) );
                ++threads_number;
            }
            if(finish) 
                break;
        }

        for(uINT i=0; i<threads_number; i++)
        {
            MULTI_TRIM_TYPE work_reward = worker[i].get();
            for(auto iter=work_reward.cbegin(); iter!=work_reward.cend(); iter++ )
            {
                if(iter->second)
                    for(auto read=iter->first->cbegin(); read!=iter->first->cend(); read++)
                        OUT << read->first << read->second;
                delete iter->first;
            }
        }
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

    if(param.threads == 1)
        trim_reads_single(param);
    else
        trim_reads_parallel(param);

    return 0;
}

