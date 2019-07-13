/*

    A Programe to produce a submatris from big matrix

*/


#include "param.h"
#include "paris.h"
#include "fasta.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <functional>
#include <cstring>
#include <numeric>
#include <chrono>       // std::chrono::system_clock
#include <random>
#include <iomanip>
#include <ctime>
#include "version.h"

using namespace std;
using namespace pan;

#define PARIS_PREPARE_VERSION "1.000"
//#define VERSION_DATE "2017-11-25"
#define COMPILE_DATE __DATE__

using std::accumulate;

Color::Modifier RED(Color::FG_RED);
Color::Modifier DEF(Color::FG_DEFAULT);

void print_usage()
{
    char buff[2000];
    const char *help_info = 
            "paris_prepare - prepare paris data to predict RNA structures\n"
            "============================================================\n"
            "\e[1mUSAGE:\e[0m\n"
            "\tparis_prepare -in input_matrix/input_sam -out output_matrix -region x_1,x_2,y_1,y_2 -norm_ratio 0.01,0.50 -min_value 0.0 \n"
            "\t           -chr chr_id [-genome fasta_file -min_overhang 5 -min_armlen 10 -strand + -rem_noise yes -norm full]\n\n"
            "\e[1mHELP:\e[0m\n"
            "\t-in: input a .matrix file or .sam file (default: guess by postfix)\n"
            "\t-out: output normalized PARIS matrix data into a file \n"
            "\t-region: 1-based coordinates of full length RNA (default: full length)\n"
            "\t-genome: a fasta file include the sequence of current matrix to suppress the non AT/GC/GT base pair scores (default: no)\n"
            "\t-chr: the sequence id of current matrix, if input a sam or provide a fasta file, it's neccessary (default: no)\n"
            "\t-rem_noise: yes|no. remove noise or not (default: yes)\n"
            "\t-norm: full|sub|none. normalize the full matrix, sub-matrix or no normalization (default: full)\n"
            "\t-min_value: define the lower bound of normalization (default: 0.0)\n"
            "\t-norm_ratio: lower_ratio,upper_ratio normalization ratio (default: 0.01,0.50)\n\n"

            "\e[1mVERSION:\e[0m\n\t%s\n"
            "\e[1mLIB VERSION:\e[0m\n\t%s\n"
            "\e[1mCOMPILE DATE:\e[0m\n\t%s\n"
            "\e[1mAUTHOR:\e[0m\n\t%s\n";

    sprintf(buff, help_info, PARIS_PREPARE_VERSION, VERSION, COMPILE_DATE, "Li Pan");
    cout << buff << endl;
}

struct Param
{
    struct Matrix_Region{ uLONG x_start, x_end, y_start, y_end; };
    enum Norm_Method{ NORM_FULL, NORM_SUB, NORM_NONE };

    string input_file;
    string output_file;

    string genome_file;
    string chr_id;

    uINT min_overhang = 5;
    uINT min_armlen = 10;
    char strand = '+';

    bool full_region = true;
    Matrix_Region region;
    string param_string;

    Norm_Method norm_method = NORM_FULL;
    bool remove_noise = true;

    double norm_upper_ratio = 0.50;
    double norm_lower_ratio = 0.01;

    double min_value = 0.0;

    operator bool()
    {
        if(input_file.empty() or output_file.empty())
            return false;

        if(not full_region)
        {
            if(region.x_start > region.x_end or region.y_start > region.y_end)
                return false;
            if(region.x_start == 0 or region.x_end == 0 or region.y_start == 0 or region.y_end == 0)
                return false;
        }

        if(not genome_file.empty() and chr_id.empty())
            return false;

        return true;
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
            }else if(not strcmp(argv[i]+1, "strand"))
            {
                has_next(argc, i);
                if( not strcmp(argv[i+1], "+") )
                    param.strand = '+';
                else if( not strcmp(argv[i+1], "-") )
                    param.strand = '-';
                else{
                    cerr << RED << "FATAL Error: unknown option value of -strand " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                i++;
            }else if(not strcmp(argv[i]+1, "in"))
            {
                has_next(argc, i);
                param.input_file = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "out"))
            {
                has_next(argc, i);
                param.output_file = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "genome"))
            {
                has_next(argc, i);
                param.genome_file = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "min_value"))
            {
                has_next(argc, i);
                param.min_value = stod(argv[i+1]);
                i++;
            }else if(not strcmp(argv[i]+1, "chr"))
            {
                has_next(argc, i);
                param.chr_id = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "region"))
            {
                has_next(argc, i);
                auto items = split(argv[i+1], ',');
                if(items.size() == 4)
                {
                    param.region.x_start = stoul(items[0]);
                    param.region.x_end = stoul(items[1]);
                    param.region.y_start = stoul(items[2]);
                    param.region.y_end = stoul(items[3]);
                }else{
                    cerr << RED << "FATAL ERROR: unknown option -region " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                param.full_region = false;
                i++;
            }else if(not strcmp(argv[i]+1, "norm_ratio"))
            {
                has_next(argc, i);
                auto items = split(argv[i+1], ',');
                if(items.size() == 2)
                {
                    param.norm_lower_ratio = stod(items[0]);
                    param.norm_upper_ratio = stod(items[1]);
                }else{
                    cerr << RED << "FATAL ERROR: unknown option -norm_ratio " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                param.full_region = false;
                i++;
            }else if(not strcmp(argv[i]+1, "rem_noise"))
            {
                has_next(argc, i);
                if( not strcmp(argv[i+1], "yes") )
                    param.remove_noise = true;
                else if( not strcmp(argv[i+1], "no") )
                    param.remove_noise = false;
                else{
                    cerr << RED << "FATAL Error: unknown option value of -rem_noise " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                i++;
            }else if(not strcmp(argv[i]+1, "norm"))
            {
                has_next(argc, i);
                if( not strcmp(argv[i+1], "full") )
                    param.norm_method = Param::NORM_FULL;
                else if( not strcmp(argv[i+1], "sub") )
                    param.norm_method = Param::NORM_SUB;
                else if( not strcmp(argv[i+1], "none") )
                    param.norm_method = Param::NORM_NONE;
                else{
                    cerr << RED << "FATAL Error: unknown option value of -norm " << argv[i+1] << DEF << endl;
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

void normalize_matrix(Matrix<double> &matrix, const double upper=0.85, const double lower=0.01, const double min_value=0)
{
    // 1. collect all values greater than min_value
    DoubleArray non_zeros;
    for(auto iter=matrix.cbegin(); iter!=matrix.cend(); iter++)
        for(double value: *iter)
        {
            if(value > min_value)
                non_zeros.push_back(value);
        }

    if(non_zeros.size() < 40)
    {
        cerr << RED << "FATAL Error: The number of values which are greater than min_value are too less: " << non_zeros.size() << DEF << endl;
        exit(-1);
    }
    std::sort(non_zeros.begin(), non_zeros.end());

    Matrix<double>::size_type nums = non_zeros.size();
    //uLONG norm_start_pos = int( nums * 16.0 / 20.0 );
    //uLONG norm_end_pos = int( nums * 18.0 / 20.0 );

    double upper_bound = non_zeros[ nums * upper ];
    double lower_bound = non_zeros[ nums * lower ];
    double divide_factor = upper_bound - lower_bound;

    clog << "\tlower_bound: " << lower_bound << "(" << uLONG(nums * lower) << ")\tupper_bound: " << upper_bound << "(" << uLONG(nums * upper) << ")" << endl;

    //double normalized_factor = std::accumulate( non_zeros.cbegin()+norm_start_pos, non_zeros.cbegin()+norm_end_pos, 0.0) / (norm_end_pos-norm_start_pos+1);
    //double norm_upper_bound = upper_bound / normalized_factor;
    //double norm_lower_bound = 0.0; //lower_bound / normalized_factor;

    for(auto iter=matrix.begin(); iter!=matrix.end(); iter++)
        for(double &value: *iter)
            if(value > upper_bound)
                value = 1.0;
            else if(value <= lower_bound)
                value = -999;
            else
                value = (value - lower_bound)/divide_factor;
}


void matrix_sub(const Param &param)
{
    sp<Fasta> p_fa;
    if(not param.genome_file.empty())
    {
        clog << currentDateTime() << "\tstart to load fasta file: " << param.genome_file << " ......" << endl;
        
        p_fa = make_shared<Fasta>(param.genome_file);
        if(not p_fa->has_chr(param.chr_id))
        {
            cerr << RED << "FATAL Error: " << param.chr_id << " not in " << param.genome_file << DEF << endl;
            exit(-1);
        }
    }

    Matrix<double> matrix;
    Rect<double> rect;

    vector<Duplex_Hang> duplex_array; uLONG chr_len;

    switch(FILE_FORMAT::guess_file_type(param.input_file))
    {
        case FILE_FORMAT::MATRIX_FILE:
            clog << currentDateTime() << "\tstart to load matrix file: " << param.input_file << " ......" << endl;
            read_matrix(param.input_file, matrix);
            break;
        case FILE_FORMAT::SAM_FILE:
            clog << currentDateTime() << "\tstart to load sam file: " << param.input_file << " ......" << endl;
            chr_len = get_chromosome_hang( param.input_file, duplex_array, param.min_overhang, param.min_armlen, param.chr_id, param.strand);
            fill_sym_matrix(matrix, duplex_array, chr_len);
            break;
        default:
            cerr << RED << "FATAL ERROR: unknown file format" << DEF << endl;
            exit(-1);
    }


    if(param.remove_noise)
    {
        clog << currentDateTime() << "\tstart to remove noise ......" << endl;

        Matrix<double> clean_matrix;
        const uLONG window_size = 5;
        remove_paris_background(matrix, clean_matrix, window_size);
        matrix = clean_matrix;
    }

    const double lower=param.norm_lower_ratio; 
    const double upper=param.norm_upper_ratio; 
    const double min_value=param.min_value;

    clog << currentDateTime() << "\tstart to normalize ......" << endl;

    if(param.full_region)
    {
        switch(param.norm_method)
        {
            case Param::NORM_FULL:
                normalize_matrix(matrix, upper, lower, min_value);
                sub_matrix(matrix, rect);
                break;
            case Param::NORM_SUB:
                sub_matrix(matrix, rect);
                normalize_matrix(rect, upper, lower, min_value);
                break;
            case Param::NORM_NONE:
                sub_matrix(matrix, rect);
                break;
        }
    }else{
        switch(param.norm_method)
        {
            case Param::NORM_FULL:
                normalize_matrix(matrix, upper, lower, min_value);
                sub_matrix(matrix, rect, param.region.x_start-1, param.region.x_end, param.region.y_start-1, param.region.y_end);
                break;
            case Param::NORM_SUB:
                sub_matrix(matrix, rect, param.region.x_start-1, param.region.x_end, param.region.y_start-1, param.region.y_end);
                normalize_matrix(rect, upper, lower, min_value);
                break;
            case Param::NORM_NONE:
                sub_matrix(matrix, rect, param.region.x_start-1, param.region.x_end, param.region.y_start-1, param.region.y_end);
                break;
        }
    }

    //rect = matrix;

    ofstream OUT(param.output_file, ofstream::out);
    if(not OUT)
    {
        cerr << RED << "FATAL Error: " << param.output_file << " is unreadable" << DEF << endl;
        exit(-1);
    }

    clog << currentDateTime() << "\tstart to write matrix ......" << endl;

    for(uLONG row_idx=0; row_idx<rect.size(); row_idx++)
    {
        for(uLONG col_idx=0; col_idx<rect[row_idx].size(); col_idx++ )
        {
            // Suppress the non-canonical base pairing
            if(p_fa)
            {
                string base_pair = p_fa->get_chr_subbseq(param.chr_id, row_idx, 1) + p_fa->get_chr_subbseq(param.chr_id, col_idx, 1);
                if( base_pair != "AT" and base_pair != "TA" and base_pair != "AU" 
                    and base_pair != "UA" and base_pair != "GC" and base_pair != "CG" 
                    and base_pair != "GT" and base_pair != "TG" and base_pair != "GU" and base_pair != "UG" )
                    OUT << setprecision(3) << -999;
                else
                    OUT << setprecision(3) << rect[row_idx][col_idx];
            }else{
                OUT << setprecision(3) << rect[row_idx][col_idx];
            }
            if(col_idx != param.region.y_end-1)
                OUT << "\t";
        }
        OUT << "\n";
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

    matrix_sub(param);

    return 0;
}




