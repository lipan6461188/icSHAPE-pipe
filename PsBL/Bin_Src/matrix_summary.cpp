/*

    A Programe to summary a matrix file

*/


#include <algorithm>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <functional>

#include "param.h"
#include "paris.h"
#include "paris_plot.h"
#include "version.h"

using namespace std;
using namespace pan;

Color::Modifier RED(Color::FG_RED);
Color::Modifier GREEN(Color::FG_GREEN);
Color::Modifier BLUE(Color::FG_BLUE);
Color::Modifier MAGENTA(Color::FG_MAGENTA);
Color::Modifier WHITE(Color::FG_WHITE);
Color::Modifier DEF(Color::FG_DEFAULT);

#define MATRIX_SUMMARY_VERSION "1.000"
#define DATE __DATE__
//"2017-12-18"

void print_usage()
{
    char buff[2000];
    const char *help_info = 
            "matrix_summary - summary and visualize a matrix file\n"
            "======================================================================\n"
            "\e[1mUSAGE:\e[0m\n"
            "\tmatrix_summary -in input_sam/input_matrix [-save file_name -bins 40 -feature max -level auto \n"
            "\t                -transform linear -chr chromosome_id -diagonal white -min_overhang 5 -min_armlen 10 \n"
            "\t                -strand + -dfile domain_file -dchar \"*\" -dfill no -rem_noise no -raw] \n"
            "\e[1mHELP:\e[0m\n"
            
            "\t\e[1mInput/Output Files: \e[0m\n"
            "\t-in: input .matrix file or .sam file  (default: guess by postfix)\n"
            "\t-chr: chr id, which must be consistant with sam file (default: no)\n"
            "\t-save: save matrix into file (default: no)\n\n"

            "\t\e[1mSam Filters: \e[0m\n"
            "\t-min_overhang: mininum overhang of duplex group, ignored when -file_type matrix (default: 5)\n"
            "\t-min_armlen: mininum arm length of each(left/right) arm, ignored when -file_type matrix (default: 10)\n"
            "\t-strand: +/- strand of reads(default: +)\n\n"

            "\t\e[1mData Clean: \e[0m\n"
            "\t-rem_noise: remove noise from raw data (default: no)\n"

            "\t\e[1mPlot Parameters: \e[0m\n"
            "\t-bins: how many bins the raw matrix will reduce to (default: 40)\n"
            "\t-feature: max/mean/sum. how to calculate a value for each bin (default: max)\n"
            "\t-level: v_1,v_2,v_3 or auto. set the color cutoff of the value of each block \n"
            "\t        auto represent 65%%,80%%,90%% of all non-zero values (default: auto)\n"
            "\t-diagonal: white/nomal. white represent white diagonal; normal represent normal color (default: white)\n"
            "\t-raw: show raw matrix (default: no)\n"
            "\t-transform: linear/log. how to transform the quantiles bar plot (default: linear)\n\n"

            "\t\e[1mDomain Parameters: \e[0m\n"
            "\t-dfile: a file to specify domains (default: no file)\n"
            "\t-dchar: domain character (default: *)\n"
            "\t-dfill: fill domain (default: no)\n\n"

            "\e[1mVERSION:\e[0m\n\t%s\n"
            "\e[1mLIB VERSION:\e[0m\n\t%s\n"
            "\e[1mCOMPILE DATE:\e[0m\n\t%s\n"
            "\e[1mAUTHOR:\e[0m\n\t%s\n";

    sprintf(buff, help_info, MATRIX_SUMMARY_VERSION, VERSION, DATE, "Li Pan");
    cout << buff << endl;
}


const string block_unicode = "\u25A0";

struct Statistic_Basic
{
    uLONG dim;

    // 0/20 - 1/20 - .... - 19/20 - 20/20
    DoubleArray quantiles; // all non-zero value quantiles: 21 parts

    uLONG non_zero_num;
};

struct Param
{
    string input_file;

    string domain_file;
    //string lower_domain_file;
    string save_file_name;

    mutable uLONG bins = 40;
    bool show_raw = false;
    FEATURE feature = FEATURE_MAX;

    mutable Color_Control cc;
    bool auto_cc = true;
    bool rem_noise = false;

    uINT min_overhang = 5;
    uINT min_armlen = 10;
    char strand = '+';
    string chr_id;

    char dchar = '*';
    bool dfill = false;

    string transform = "linear";
    bool diagonal = true;

    bool has_domain()const { return not domain_file.empty(); }

    operator bool()
    { 
        return not input_file.empty() and bins > 5 and cc.valid();
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
            }else if(not strcmp(argv[i]+1, "bins"))
            {
                has_next(argc, i);
                param.bins = stoul(string(argv[i+1]));
                i++;
            }else if(not strcmp(argv[i]+1, "feature"))
            {
                has_next(argc, i);
                if( not strcmp(argv[i+1], "max") )
                {
                    param.feature = FEATURE_MAX;
                }else if( not strcmp(argv[i+1], "mean") )
                {
                    param.feature = FEATURE_MEAN;
                }else if( not strcmp(argv[i+1], "sum") )
                {
                    param.feature = FEATURE_SUM;
                }else{
                    cerr << RED << "FATAL ERROR: unknown -feature option: " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                i++;
            }else if(not strcmp(argv[i]+1, "in"))
            {
                has_next(argc, i);
                param.input_file = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "save"))
            {
                has_next(argc, i);
                param.save_file_name = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "level"))
            {
                has_next(argc, i);
                if( not strcmp(argv[i+1], "auto") )
                {
                    param.auto_cc = true;
                }else{
                    auto levels = split(argv[i+1], ',');
                    if(levels.size() != 3)
                    {
                        cerr << RED << "FATAL ERROR: bad -level option: " << argv[i+1] << DEF << endl;
                        exit(-1);
                    }
                    param.cc.level_up_1 = param.cc.level_dn_1 = stod(levels[0]);
                    param.cc.level_up_2 = param.cc.level_dn_2 = stod(levels[1]);
                    param.cc.level_up_3 = param.cc.level_dn_3 = stod(levels[2]);
                    param.auto_cc = false;
                }
                i++;
            }else if(not strcmp(argv[i]+1, "dfill"))
            {
                has_next(argc, i);
                if( not strcmp(argv[i+1], "yes") )
                {
                    param.dfill = true;
                }else if( not strcmp(argv[i+1], "no") ){
                    param.dfill = false;
                }else{
                    cerr << RED << "FATAL ERROR: bad -dfill option: " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                i++;
            }else if(not strcmp(argv[i]+1, "h"))
            {
                print_usage();
                exit(-1);
            }else if(not strcmp(argv[i]+1, "transform"))
            {
                has_next(argc, i);
                param.transform = argv[i+1];
                if(param.transform != "linear" and param.transform != "log")
                {
                    cerr << RED << "FATAL ERROR: bad -transform option: " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                i++;
            }else if(not strcmp(argv[i]+1, "diagonal"))
            {
                has_next(argc, i);
                if( not strcmp( argv[i+1], "white") )
                {
                    param.diagonal = true;
                }else if( not strcmp( argv[i+1], "normal") )
                {
                    param.diagonal = false;
                }else{
                    cerr << RED << "FATAL ERROR: bad -diagonal option: " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                i++;
            }else if(not strcmp(argv[i]+1, "chr"))
            {
                has_next(argc, i);
                param.chr_id = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "dchar"))
            {
                has_next(argc, i);
                param.dchar = argv[i+1][0];
                i++;
            }else if(not strcmp(argv[i]+1, "raw"))
            {
                //has_next(argc, i);
                param.show_raw = true;
                //i++;
            }else if(not strcmp(argv[i]+1, "rem_noise"))
            {
                has_next(argc, i);
                if( not strcmp(argv[i+1], "yes") )
                {
                    param.rem_noise = true;
                }else if( not strcmp(argv[i+1], "no") ){
                    param.rem_noise = false;
                }else{
                    cerr << RED << "FATAL ERROR: bad -rem_noise option: " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                i++;
            }else if(not strcmp(argv[i]+1, "dfile"))
            {
                has_next(argc, i);
                param.domain_file = argv[i+1];
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

namespace TRANSFORM
{
    double linear(double k, double b, double raw)
    {
        return raw * k + b;
    }

    double log(double base, double raw)
    {
        if(base != 2)
            return log2(raw)/log2(base);
        else
            return log2(raw);
    }
}

void print_statistics(const Statistic_Basic &sb, 
    const Color_Control &cc, 
    const uLONG indent=0,
    function<double(double)> transform=[](double x){return x;})
{
    double step = ( transform(sb.quantiles[sb.quantiles.size()-2]) - transform(sb.quantiles.front()) ) / 10;
    uLONGArray blocks_count;
    for_each(sb.quantiles.cbegin(), sb.quantiles.cend(), [&](double value){ blocks_count.push_back( min( transform(value) /step, 10.0) ); });
    //cout << blocks_count << endl;

    //const uLONG width = 20;
    const string indent_string(indent, ' ');

    cout << DEF;
    for(uLONG i=0; i<10; i++)
    {
        cout << DEF << indent_string << "| ";
        for(uLONG j=0; j<21; j++)
        {

            if(j>=19)
                cout << cc.color_3;
            else if(j>=16)
                cout << cc.color_2;
            else if(j>=13)
                cout << cc.color_1;
            else
                cout << DEF;

            if(blocks_count[j] >= 10-i)
                cout << block_unicode+" ";
            else
                cout << "  ";
        }
        cout << "\n";
    }

    cout << DEF << indent_string;
    for(uLONG i=0; i<22; i++)
        cout << "- ";
    cout << "\n";
    
    cout << "  " << indent_string;
    for(uLONG j=0; j<21; j++)
    {
        if(j>=19)
            cout << cc.color_3;
        else if(j>=16)
            cout << cc.color_2;
        else if(j>=13)
            cout << cc.color_1;
        else
            cout << DEF;

        if(j<=8)
            cout << j + 1 << ' ';
        else
            cout << char('A'+j-9) << ' ';
    }
    cout << "\n";
}

void print_legend(const Color_Control &cc)
{
    cc.print_legend();
}

/*
    level_3 95% of sites
    level_2 80% of sites
    level_1 65% of sites
*/
void auto_color_control(const Statistic_Basic &sb, Color_Control &cc)
{
    cc.level_up_1 = cc.level_dn_1 = sb.quantiles.at(13); // 13 means 6.5
    cc.level_up_2 = cc.level_dn_2 = sb.quantiles.at(16); // 16 means 8.0
    cc.level_up_3 = cc.level_dn_3 = sb.quantiles.at(19); // 19 means 9.5
}

void auto_statistic(const Matrix<double> &matrix, Statistic_Basic &sb)
{
    auto matrix_size = matrix.size();

    sb.dim = matrix_size;
    sb.quantiles.clear();

    DoubleArray non_zeros;
    for(uLONG idx=0; idx<matrix_size; idx++)
        for(uLONG idy=0; idy<matrix_size; idy++)
            if(matrix[idx][idy] > 0.0)
                non_zeros.push_back(matrix[idx][idy]);

    std::sort(non_zeros.begin(), non_zeros.end());
    if(non_zeros.size() < 30)
    {
        cerr << RED << "FATAL Error: non-zero values in matrix are too less" << endl;
        exit(-1);
    }

    sb.non_zero_num = non_zeros.size();
    for(uLONG i=0; i<20; i++)
        sb.quantiles.push_back( non_zeros.at(i/20.0 * sb.non_zero_num) );
    sb.quantiles.push_back( non_zeros.back() );
}

void save_matrix(const Matrix<double> &matrix, const string &save_file)
{
    ofstream OUT(save_file, ofstream::out);
    if(not OUT)
    {
        cerr << RED << "FATAL Error: " << save_file << "is unwritable" << DEF << endl;
        exit(-1);
    }
    for(uINT i=0; i<matrix.size(); i++)
    {
        for(uINT j=0; j<matrix.size(); j++)
        {
            OUT << matrix[i][j];
            if(j != matrix.size()-1)
                OUT << "\t";
        }
        OUT << "\n";
    }
    OUT.close();
}

void matrix_summary(const Param &param)
{
    using namespace std::placeholders;
    Matrix<double> matrix, comp_matrix;
    vector<Duplex_Hang> duplex_array;
    uLONG chr_len;
    RegionArray domain_regions, c_domain_regions; //lower_regions, c_upper_regions, c_lower_regions;


    switch(FILE_FORMAT::guess_file_type(param.input_file))
    {
        case FILE_FORMAT::MATRIX_FILE:
            read_matrix(param.input_file, matrix);
            break;
        case FILE_FORMAT::SAM_FILE:
            chr_len = get_chromosome_hang( param.input_file, duplex_array, param.min_overhang, param.min_armlen, param.chr_id, param.strand);
            fill_sym_matrix(matrix, duplex_array, chr_len);
            duplex_array.clear();
            break;
        default:
            cerr << RED << "FATAL ERROR: unknown file format" << DEF << endl;
            exit(-1);
    }

    if(param.rem_noise)
    {
        const uINT window(5);
        Matrix<double> clean_matrix;
        remove_paris_background(matrix, clean_matrix, window);
        matrix = clean_matrix;
    }

    if(param.show_raw)
    {
        param.bins = matrix.size();
        comp_matrix = matrix;
    }else{
        compress_matrix(matrix, comp_matrix, param.bins, param.feature);
    }
    

    Statistic_Basic sb;
    auto_statistic(matrix, sb);

    //Color_Control cc;
    if(param.auto_cc)
        auto_color_control(sb, param.cc);

    if( param.has_domain() )
    {
        read_domain_file(param.domain_file, domain_regions);
        compress_regions(domain_regions, c_domain_regions, matrix.size(), param.bins);
        //cout << c_upper_regions << endl;
    }

    print_heatmap(comp_matrix, param.cc, param.diagonal, 0, "\u25A0", STD_Color::WHITE, c_domain_regions, c_domain_regions, param.dchar, param.dfill);

    if( not param.save_file_name.empty() )
        save_matrix(comp_matrix, param.save_file_name);

/*
void print_heatmap( const Matrix<double> &matrix, 
                    const Color_Control &cc,
                    bool white_diagonal = false,
                    const uLONG indent = 0,
                    const string block_unicode = "\u25A0",
                    const Color::Modifier domain_color = STD_Color::WHITE,
                    const RegionArray &upper_region = RegionArray(),
                    const RegionArray &lower_region = RegionArray());
*/

    print_legend(param.cc);

    // print title
    const uLONG bar_plot_width = 40;
    char buff[2000];
    cout << DEF;
    sprintf(buff, "20 quantiles of all non-zero values(dim: %lu*%lu max value: %.3f nums: %lu(%.2f%%))", 
            sb.dim, sb.dim, sb.quantiles.back(), sb.non_zero_num, 100.0*sb.non_zero_num/(sb.dim*sb.dim));
    uLONG title_len = strlen(buff);
    uLONG indent = 0, title_indent = 0;

    if(title_len > bar_plot_width)
    {
        indent = (title_len - bar_plot_width) / 2;
    }else{
        title_indent = (bar_plot_width - title_len) / 2;
    }

    cout << string(title_indent, ' ') << buff << "\n";
    if(param.transform == "linear")
        print_statistics(sb, param.cc, indent);
    else if(param.transform == "log")
        print_statistics(sb, param.cc, indent, std::bind(TRANSFORM::log, 2, _1));
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

    matrix_summary(param);

    return 0;
}




