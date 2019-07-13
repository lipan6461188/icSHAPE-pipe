/*

    A Programe to covert sam file into tab-seperated dg file

*/

#include "paris.h"
#include "param.h"
#include "fasta.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <sstream>
#include "version.h"

using namespace std;
using namespace pan;

#define SAM2FQ_VERSION "1.000"
#define DATE __DATE__
#define WARNING "The input of sam2dg must be sorted by read id"

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
            "\tsam2fq -in input_sam -out output_fq [ -quick ] \n"
            "\e[1mHELP:\e[0m\n"
            "\t-in: input sam file\n"
            "\t-out: output fastq file\n"
            "\t-quick: quick mode - don't remove duplicated reads to speed up \n\n"

            "\e[1mWARNING:\e[0m\n\t%s\n"
            "\e[1mVERSION:\e[0m\n\t%s\n"
            "\e[1mLIB VERSION:\e[0m\n\t%s\n"
            "\e[1mCOMPILE DATE::\e[0m\n\t%s\n"
            "\e[1mAUTHOR:\e[0m\n\t%s\n";

    ostringstream warning;
    warning << YELLOW << WARNING << DEF;

    sprintf(buff, help_info, warning.str().c_str(), SAM2FQ_VERSION, VERSION, DATE, "Li Pan");
    cout << buff << endl;
}

struct Param
{
    string input_sam;
    string output_fq;

    bool quick_mode = false;

    operator bool(){ return input_sam.empty() or output_fq.empty() ? false : true; }
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
                param.input_sam = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "out"))
            {
                has_next(argc, i);
                param.output_fq = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "h"))
            {
                print_usage();
                exit(-1);
            }else if(not strcmp(argv[i]+1, "quick"))
            {
                param.quick_mode = true;
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


void sam2fq(const Param &param)
{
    ifstream IN(param.input_sam, ifstream::in);
    ofstream OUT(param.output_fq, ofstream::out);

    if(not IN)
    {
        cerr << RED << "Fatal Error: " << param.input_sam << " cannot be readable" << endl;
        exit(-1);
    }

    if(not OUT)
    {
        cerr << RED << "Fatal Error: " << param.output_fq << " cannot be writable" << endl;
        exit(-1);
    }

    Sam_Head sam_head;
    read_sam_head(IN, sam_head);

    if(param.quick_mode)
    {
        Sam_Record read_record;
        while(not IN.eof())
        {
            bool success = read_a_sam_record(IN, read_record);
            if(not success)
                break;
            if(read_is_reverse(read_record))
            {
                OUT << "@" << read_record.read_id << "\n" << reverse_comp(read_record.read_seq) << "\n+\n" << reverse_string(read_record.read_quality) << "\n";
            }else{
                OUT << "@" << read_record.read_id << "\n" << read_record.read_seq << "\n+\n" << read_record.read_quality << "\n";
            }
        }
    }else{
        vector<Sam_Record> read_records;
        while(not IN.eof())
        {
            bool success = read_a_read_record(IN, read_records);
            if(not success)
                break;
            if(read_is_reverse(read_records[0]))
            {
                OUT << "@" << read_records[0].read_id << "\n" << reverse_comp(read_records[0].read_seq) << "\n+\n" << reverse_string(read_records[0].read_quality) << "\n";
            }else{
                OUT << "@" << read_records[0].read_id << "\n" << read_records[0].read_seq << "\n+\n" << read_records[0].read_quality << "\n";
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

    sam2fq(param);

    return 0;
}












