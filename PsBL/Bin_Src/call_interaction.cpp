
/*

    A Programe to call PARIS interaction regions from a sam file or a matrix

*/

#include "paris.h"
#include "param.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include "version.h"

using namespace std;
using namespace pan;

#define CALL_INTTERACTION_VERSION "1.000"
#define VERSION_DATE "2017-11-2"
#define COMPILE_DATE __DATE__

Color::Modifier RED(Color::FG_RED);
Color::Modifier DEF(Color::FG_DEFAULT);

void print_usage()
{
    char buff[2000];
    const char *help_info = 
            "call_interaction - call PARIS interaction regions from a sam file or a matrix\n"
            "=============================================================================\n"
            "\e[1mUSAGE:\e[0m\n"
            "\tcall_interaction -in input_sam/input_matrix -chr chr_id -out output_txt [-min_overhang 5 -min_armlen 10 -min_dist 100\n"
            "\t                 -min_window_size 50 -max_window_size 200 -percep_threshold 100 -extend_threshold 20 -file_type sam]\n"
            "\e[1mHELP:\e[0m\n"
            "\t-min_overhang: mininum overhang of duplex group(default: 5)\n"
            "\t-min_armlen: mininum arm length of each(left/right) arm(default: 10)\n"
            "\t-min_dist: mininum distance of interaction(default: 100)\n"
            "\t-min_window_size: mininum window size of scanning(default: 50)\n"
            "\t-max_window_size: maximun window size of scanning(default: 200)\n"
            "\t-percep_threshold: perception cutoff of scanning(default: 100)\n"
            "\t-extend_threshold: extending cutoff of scanning(default: 20)\n"
            "\t-file_type: input file type -- sam or matrix(default: sam) \n"
            "\e[1mVERSION:\e[0m\n\t%s\n"
            "\e[1mLIB VERSION:\e[0m\n\t%s\n"
            "\e[1mVERSION DATE:\e[0m\n\t%s\n"
            "\e[1mCOMPILE DATE:\e[0m\n\t%s\n"
            "\e[1mAUTHOR:\e[0m\n\t%s\n";
    sprintf(buff, help_info, CALL_INTTERACTION_VERSION, VERSION, VERSION_DATE, COMPILE_DATE, "Li Pan");
    cout << buff << endl;
}


struct Param
{
    enum FILE_TYPE {SAM_FILE, MATRIX_FILE};

    string input_file;
    string output_txt;

    string chr_id;

    uINT min_overhang = 5;
    uINT min_armlen = 10;
    uINT min_dist = 100;
    uINT min_window_size = 50;
    uINT max_window_size = 200;

    double percep_threshold = 100;
    double extend_threshold = 20;

    FILE_TYPE file_type = SAM_FILE;

    string param_string;

    operator bool(){ 
        if(file_type == SAM_FILE)
            return chr_id.empty() or input_file.empty() or output_txt.empty() ? false : true;
        else if(file_type == MATRIX_FILE)
            return input_file.empty() or output_txt.empty() ? false : true;
        else{
            return false;
        }
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
                param.output_txt = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "h"))
            {
                print_usage();
                exit(0);
            }else if(not strcmp(argv[i]+1, "min_dist"))
            {
                has_next(argc, i);
                param.min_dist = stoul(string(argv[i+1]));
                i++;
            }else if(not strcmp(argv[i]+1, "min_window_size"))
            {
                has_next(argc, i);
                param.min_window_size = stoul(string(argv[i+1]));
                i++;
            }else if(not strcmp(argv[i]+1, "max_window_size"))
            {
                has_next(argc, i);
                param.max_window_size = stoul(string(argv[i+1]));
                i++;
            }else if(not strcmp(argv[i]+1, "percep_threshold"))
            {
                has_next(argc, i);
                param.percep_threshold = stod(string(argv[i+1]));
                i++;
            }else if(not strcmp(argv[i]+1, "extend_threshold"))
            {
                has_next(argc, i);
                param.extend_threshold = stod(string(argv[i+1]));
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


void call_interaction(const Param &param)
{
    vector<Duplex_Hang> duplex_array;
    uLONG chr_len;
    
    Matrix<double> matrix;
    vector<InterRegion> interact_regions;

    if(param.file_type == Param::SAM_FILE)
    {
        try{
            chr_len = get_chromosome_hang( param.input_file, duplex_array, param.min_overhang, param.min_armlen, param.chr_id);
            fill_sym_matrix(matrix, duplex_array, chr_len);
        }catch(runtime_error e)
        {
            cerr << RED << "FATAL Error: " << e.what() << DEF << endl;
            print_usage();
            exit(-1);
        }
    }else{
        try{
            chr_len = read_matrix(param.input_file, matrix);
        }catch(runtime_error e)
        {
            cerr << RED << "FATAL Error: " << e.what() << DEF << endl;
            print_usage();
            exit(-1);
        }
    }

    scan_interaction<double>(   matrix, 
                        interact_regions,
                        param.min_dist, 
                        param.min_window_size, 
                        param.max_window_size, 
                        param.percep_threshold, 
                        param.extend_threshold);

    clog << "Report: " << interact_regions.size() << " are found finally" << endl;

    if(not interact_regions.empty())
    {
        ofstream OUT(param.output_txt, ofstream::out);
        if(not OUT)
        {
            cerr << RED << "FATAL Error: " + param.output_txt + " cannot be writebale" << DEF << endl;
            exit(-1);
        }
        OUT << "#" << param.param_string << "\n";
        OUT << interact_regions;
        OUT.close();
    }
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

    call_interaction(param);

    return 0;
}



