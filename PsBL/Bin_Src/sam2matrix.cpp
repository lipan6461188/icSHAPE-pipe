/*

    A Programe to convert sam file into matrix

*/


//#include "paris_plot.h"
#include "paris.h"
#include "param.h"
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


#define CALL_SAM2MATRIX_VERSION "1.000"
#define DATE __DATE__

using std::accumulate;

Color::Modifier RED(Color::FG_RED);
Color::Modifier DEF(Color::FG_DEFAULT);

void print_usage()
{
    char buff[2000];
    const char *help_info = 
            "sam2matrix - remove the background in PARIS data\n"
            "===============================================================================================\n"
            "\e[1mUSAGE:\e[0m\n"
            "\tparis_backround -in input_sam -chr chr_id -out output_matrix [-min_overhang 5 -min_armlen 10 -strand +]\n"
            "\e[1mHELP:\e[0m\n"
            "\t-min_overhang: mininum overhang of duplex group(default: 5)\n"
            "\t-min_armlen: mininum arm length of each(left/right) arm(default: 10)\n"
            "\t-strand: strand of reads(+/-) (default: +)\n"
            "\e[1mVERSION:\e[0m\n\t%s\n"
            "\e[1mLIB VERSION:\e[0m\n\t%s\n"
            "\e[1mDATE:\e[0m\n\t%s\n"
            "\e[1mAUTHOR:\e[0m\n\t%s\n";
    sprintf(buff, help_info, CALL_SAM2MATRIX_VERSION, VERSION, DATE, "Li Pan");
    cout << buff << endl;
}

struct Param
{

    string input_file;
    string output_matrix;
    string chr_id;

    uINT min_overhang = 5;
    uINT min_armlen = 10;
    char strand = '+';

    string param_string;

    operator bool(){ return (input_file.empty() or output_matrix.empty() or chr_id.empty()) ? false : true; }
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
            }else if(not strcmp(argv[i]+1, "chr"))
            {
                has_next(argc, i);
                param.chr_id = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "in"))
            {
                has_next(argc, i);
                param.input_file = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "out"))
            {
                has_next(argc, i);
                param.output_matrix = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "strand"))
            {
                has_next(argc, i);
                if( not strcmp(argv[i+1], "+") )
                    param.strand = '+';
                else if( not strcmp(argv[i+1], "-") )
                    param.strand = '-';
                else{
                    cerr << "FATAL Error: unknown option value of -strand " << argv[i+1] << endl;
                    exit(-1);
                }
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


void sam2matrix(const Param &param)
{
    vector<Duplex_Hang> duplex_array;
    uLONG chr_len;
    ofstream OUT;
    
    Matrix<double> matrix;

    try{
        chr_len = get_chromosome_hang( param.input_file, duplex_array, param.min_overhang, param.min_armlen, param.chr_id, param.strand);
        fill_sym_matrix(matrix, duplex_array, chr_len);
    }catch(runtime_error e)
    {
        cerr << RED << "FATAL Error: " << e.what() << DEF << endl;
        print_usage();
        exit(-1);
    }

    OUT.open(param.output_matrix, ofstream::out);
    OUT << matrix;
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

    sam2matrix(param);

    return 0;
}

