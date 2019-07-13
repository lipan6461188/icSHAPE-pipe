
/*

    A Programe to convert stockholm file to fasta file

*/


//#include "paris_plot.h"
#include "align.h"
#include "param.h"
#include "fasta.h"
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

#define CALL_STO2FA_VERSION "1.000"
#define DATE __DATE__

using std::accumulate;

Color::Modifier RED(Color::FG_RED);
Color::Modifier DEF(Color::FG_DEFAULT);

void print_usage()
{
    char buff[2000];
    const char *help_info = 
            "sto2fa - convert stockholm file to fasta file\n"
            "=============================================================================\n"
            "\e[1mUSAGE:\e[0m\n"
            "\tsto2fa -in input_sto -out output_fasta\n"
            "\e[1mVERSION:\e[0m\n\t%s\n"
            "\e[1mLIB VERSION:\e[0m\n\t%s\n"
            "\e[1mDATE:\e[0m\n\t%s\n"
            "\e[1mAUTHOR:\e[0m\n\t%s\n";
    sprintf(buff, help_info, CALL_STO2FA_VERSION, VERSION, DATE, "Li Pan");
    cout << buff << endl;
}


struct Param
{

    string input_sto;
    string output_fasta;

    string param_string;

    operator bool(){ return (input_sto.empty() or output_fasta.empty()) ? false : true; }
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
                param.input_sto = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "out"))
            {
                has_next(argc, i);
                param.output_fasta = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "h"))
            {
                print_usage();
                exit(-1);
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

void sto2fa(const Param &param)
{
    Multi_Align ma(param.input_sto);
    const StringArray& chr_ids = ma.keys();

    ofstream OUT(param.output_fasta, ofstream::out);
    for(const string &chr_id: chr_ids)
    {
        try{
            string anno = ma.get_annotation(chr_id);
            OUT << ">" << chr_id << "\t" << ma.get_annotation(chr_id) << "\n" << flat_seq(ma.get_align_seq(chr_id)) << "\n";
        }catch(out_of_range e)
        {
            OUT << ">" << chr_id << "\n" << flat_seq(ma.get_align_seq(chr_id)) << "\n";
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

    sto2fa(param);

    return 0;
}





