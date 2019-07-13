
/*

    A Programe to remove the background in PARIS data

*/


//#include "paris_plot.h"
#include "paris.h"
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
#include "version.h"

using namespace std;
using namespace pan;

#define CALL_BACKGROUND_VERSION "1.000"
#define DATE __DATE__

using std::accumulate;

Color::Modifier RED(Color::FG_RED);
Color::Modifier DEF(Color::FG_DEFAULT);

void print_usage()
{
    char buff[2000];
    const char *help_info = 
            "paris_backround - remove the background in PARIS data\n"
            "=============================================================================\n"
            "\e[1mUSAGE:\e[0m\n"
            "\tparis_backround -in input_sam/input_matrix -chr chr_id -out output_matrix -method estimate \n"
            "\t                [-min_overhang 5 -min_armlen 10 -ratio 0.6 -surround 5 -file_type sam]\n"
            "\e[1mHELP:\e[0m\n"
            "\t-method: estimate or quantile(default: estimate)\n"
            "\t-min_overhang: mininum overhang of duplex group(default: 5)\n"
            "\t-min_armlen: mininum arm length of each(left/right) arm(default: 10)\n"
            "\t-ratio: ratio(0-1) of quantile(default: 0.6)\n"
            "\t-surround: estimate the background from surrounding nucleotide base interactiob(default: 5)\n"
            "\t-file_type: input file type -- sam or matrix(default: sam) \n"
            "\e[1mVERSION:\e[0m\n\t%s\n"
            "\e[1mLIB VERSION:\e[0m\n\t%s\n"
            "\e[1mDATE:\e[0m\n\t%s\n"
            "\e[1mAUTHOR:\e[0m\n\t%s\n";
    sprintf(buff, help_info, CALL_BACKGROUND_VERSION, VERSION, DATE, "Li Pan");
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




template<typename T>
void quantile_method(const Matrix<T> &raw_matrix, 
                        Matrix<T> &matrix,
                        double ratio=0.6)
{
    if(ratio <= 0 or ratio >= 1)
    {
        throw runtime_error("Invalid BackGround Ratio");
    }

    using size_type = typename Matrix<T>::size_type;

    matrix = raw_matrix;
    size_type matrix_size = matrix.size();
    DoubleArray background;

    //vector<T> tmp_vec;
    for(decltype(matrix_size) idx=0; idx<matrix_size; idx++)
    {
        vector<T> tmp_vec = matrix.at(idx);
        sort(tmp_vec.begin(), tmp_vec.end());
        auto lower_bound = tmp_vec.at( matrix_size*ratio );
        background.push_back(lower_bound);
    }

    T zero = T();

    for(size_type idx=0; idx<matrix_size; idx++)
        for(size_type idy=0; idy<matrix_size; idy++)
        {
            matrix.at(idx).at(idy) -= background.at(idx) + background.at(idy);
            if( matrix.at(idx).at(idy) < zero )
                matrix.at(idx).at(idy) = zero;
        }
}

void paris_backround(const Param &param)
{
    vector<Duplex_Hang> duplex_array;
    uLONG chr_len;
    ofstream OUT;
    
    Matrix<double> raw_matrix, norm_matrix;

    if(param.file_type == Param::SAM_FILE)
    {
        try{
            chr_len = get_chromosome_hang( param.input_file, duplex_array, param.min_overhang, param.min_armlen, param.chr_id);
            fill_sym_matrix(raw_matrix, duplex_array, chr_len);
        }catch(runtime_error e)
        {
            cerr << RED << "FATAL Error: " << e.what() << DEF << endl;
            print_usage();
            exit(-1);
        }
    }else{
        try{
            chr_len = read_matrix(param.input_file, raw_matrix);
        }catch(runtime_error e)
        {
            cerr << RED << "FATAL Error: " << e.what() << DEF << endl;
            print_usage();
            exit(-1);
        }
    }

    if(param.method == Param::QUANTILE_METHOD)
        quantile_method(raw_matrix, norm_matrix, param.ratio);
    else if(param.method == Param::ESTIMATE_METHOD)
        remove_paris_background(raw_matrix, norm_matrix, param.surround);

    OUT.open(param.output_matrix, ofstream::out);
    OUT << norm_matrix;
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

    paris_backround(param);

    return 0;
}

