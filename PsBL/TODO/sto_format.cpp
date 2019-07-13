/*
	
	MINI Programe -- 
    A Programe to format a sto file with a transcript

*/


//#include "paris_plot.h"
#include "pan_type.h"
#include "param.h"

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


#define CALL_FORMATSTO_VERSION "1.101"
#define DATE "2017-11-2"

using std::accumulate;

Color::Modifier RED(Color::FG_RED);
Color::Modifier DEF(Color::FG_DEFAULT);

void print_usage()
{
    char buff[2000];
    const char *help_info = 
            "sto_format - format a sto file with a transcript\n"
            "=============================================================================\n"
            "\e[1mUSAGE:\e[0m\n"
            "\tsto_format -in input_sto -chr chr_id -out output_sto\n"
            "\e[1mVERSION:\e[0m\n\t%s\n"
            "\e[1mLIB VERSION:\e[0m\n\t%s\n"
            "\e[1mDATE:\e[0m\n\t%s\n"
            "\e[1mAUTHOR:\e[0m\n\t%s\n";
    sprintf(buff, help_info, CALL_FORMATSTO_VERSION, VERSION, DATE, "Li Pan");
    cout << buff << endl;
}


struct Param
{
    enum FILE_TYPE {SAM_FILE, MATRIX_FILE};
    enum METHOD { QUANTILE_METHOD, ESTIMATE_METHOD };

    string input_file;
    string output_matrix;

    string chr_id;

    uINT min_overhang = 5;
    uINT min_armlen = 10;
    double ratio = 0.6;
    uINT surround = 5;

    FILE_TYPE file_type = SAM_FILE;
    METHOD method = ESTIMATE_METHOD;

    string param_string;

    operator bool(){ 
        bool cond1, cond2;
        if(file_type == SAM_FILE)
            cond1 = chr_id.empty() or input_file.empty() or output_matrix.empty() ? false : true;
        else if(file_type == MATRIX_FILE)
            cond1 = input_file.empty() or output_matrix.empty() ? false : true;
        else{
            cond1 = false;
        }
        cond2 = (ratio < 1) and (ratio > 0);
        return cond1 and cond2;
    }
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
            }else if(not strcmp(argv[i]+1, "h"))
            {
                print_usage();
                exit(0);
            }else if(not strcmp(argv[i]+1, "ratio"))
            {
                has_next(argc, i);
                param.ratio = stod(string(argv[i+1]));
                i++;
            }else if(not strcmp(argv[i]+1, "surround"))
            {
                has_next(argc, i);
                param.surround = stoul(string(argv[i+1]));
                i++;
            }else if(not strcmp(argv[i]+1, "method"))
            {
                has_next(argc, i);
                if( not strcmp(argv[i+1], "estimate") )
                    param.method = Param::ESTIMATE_METHOD;
                else if( not strcmp(argv[i+1], "quantile") )
                    param.method = Param::QUANTILE_METHOD;
                else{
                    cerr << RED << "FATAL Error: Unknown Method(-method): " << argv[i+1] << DEF << endl;
                    print_usage();
                    exit(-1);
                }
                i++;
            }else if(not strcmp(argv[i]+1, "file_type"))
            {
                has_next(argc, i);
                if( not strcmp(argv[i+1], "sam") )
                    param.file_type = Param::SAM_FILE;
                else if( not strcmp(argv[i+1], "matrix") )
                    param.file_type = Param::MATRIX_FILE;
                else{
                    cerr << RED << "FATAL Error: Unknown File Type(-file_type): " << argv[i+1] << DEF << endl;
                    print_usage();
                    exit(-1);
                }
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
    param.param_string = param_string(argc, argv);
    return param;
}




