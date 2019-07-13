/*

    A Programe to tag the number of mismatched bases between reads and reference

*/


//#include "paris_plot.h"
#include "param.h"
#include "fasta.h"
#include "sam.h"
#include "version.h"

//#include <QGuiApplication>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <functional>
#include <cstring>
#include <numeric>

using namespace std;
using namespace pan;

#define CALL_SAM_MISMATCH_VERSION "1.000"
#define DATE __DATE__

using std::accumulate;

Color::Modifier RED(Color::FG_RED);
Color::Modifier DEF(Color::FG_DEFAULT);

void print_usage()
{
    char buff[2000];
    const char *help_info = 
            "sam_mismatch - tag the number of mismatched bases between reads and reference\n"
            "=============================================================================\n"
            "\e[1mUSAGE:\e[0m\n"
            "\tsam_mismatch -in input_sam -out output_sam -genome ref_seq.fa [ -tag MM ]\n"
            "\e[1mHELP:\e[0m\n"
            "\t-tag: the tag name of mismatched base(default: MM)\n"
            "\e[1mHELP:\e[0m\n"
            "\e[1mVERSION:\e[0m\n\t%s\n"
            "\e[1mLIB VERSION:\e[0m\n\t%s\n"
            "\e[1mDATE:\e[0m\n\t%s\n"
            "\e[1mAUTHOR:\e[0m\n\t%s\n";
    sprintf(buff, help_info, CALL_SAM_MISMATCH_VERSION, VERSION, DATE, "Li Pan");
    cout << buff << endl;
}

struct Param
{
    string input_sam;
    string output_sam;
    string genome_file;

    string tag = "MM";

    string param_string;

    operator bool(){ return (input_sam.empty() or output_sam.empty() or genome_file.empty()) ? false : true; }
};

void has_next(int argc, int current)
{
    if(current + 1 >= argc)
    {
        cerr << RED << "FATAL ERROR: Parameter Error" << DEF << DEF << endl;
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
                param.output_sam = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "genome"))
            {
                has_next(argc, i);
                param.genome_file = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "tag"))
            {
                has_next(argc, i);
                param.tag = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "h"))
            {
                print_usage();
                exit(0);
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
    param.param_string = param_string(argc, argv);
    return param;
}

/*
    Cigar:

ref     √   x   √   x
read    √   √   x   x
Cigar   M   I   D   P
        X   S   N
            H

*/

uINT MM_number(const Sam_Record &record, const Fasta &fasta)
{
    //string seq_frag(const string &key, uLONG start, uLONG len=string::npos) const;

    auto chr_id = record.chr_id;
    auto start = record.pos - 1;

    /*
    if(not fasta.has(chr_id))
    {
        cerr << RED << chr_id << " not in reference file" << endl;
        return 0;
    }
    */

    const string &ref_seq = fasta.get_chr_seq(chr_id);
    const string &read_seq = record.read_seq;

    vector<int> cigarLen;
    vector<char> cigarAlpha;
    split_cigar(record.cigar, cigarLen, cigarAlpha);

    uINT mm_number(0);

    uINT read_index(0);
    for(size_t i=0; i<cigarLen.size(); i++)
    {
        switch(cigarAlpha[i])
        {
            case 'M': case 'X':
                for(auto iter=read_seq.cbegin()+read_index; iter!=read_seq.cbegin()+read_index+cigarLen[i]; iter++)
                    if(*iter != ref_seq[start++])
                        ++mm_number;
                read_index += cigarLen[i];
                break;
            case 'I': case 'S': case 'H':
                read_index += cigarLen[i];
                break;
            case 'D': case 'N':
                start += cigarLen[i];
                break;
            case 'P':
                break;
            default:
                cerr << "Undefined Cigar Code: " << cigarAlpha[i] << endl;
        }
    }
    return mm_number;
}

void tag_MM_from_sam(const Param &param)
{
    //using size_type = vector<Sam_Record>::size_type;

    Fasta fasta(param.genome_file);

    //clog << fasta.keys() << endl;

    ifstream IN(param.input_sam, ifstream::in);
    ofstream OUT(param.output_sam, ofstream::out);

    uLONGLONG line_count = 0;

    Sam_Head sam_head;
    read_sam_head(IN, sam_head);
    OUT << sam_head;

    Sam_Record read_record;
    while(not IN.eof())
    {
        bool success = read_a_sam_record(IN, read_record);
        
        if(not success)
            break;
        
        if(not read_is_mapped(read_record))
            continue;
        
        uINT mm_number;
        try{
            mm_number = MM_number(read_record, fasta);
        }catch(out_of_range e)
        {
            cerr << read_record.chr_id << " not in reference file" << endl;
            continue;
        }
        
        read_record.attributes.push_back(param.tag+":i:"+to_string(mm_number));
        OUT << read_record;

        ++line_count;
        if(line_count % 100000 == 0)
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

    tag_MM_from_sam(param);

    return 0;
}











