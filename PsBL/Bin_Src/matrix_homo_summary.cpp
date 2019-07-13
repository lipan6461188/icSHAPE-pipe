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
#include "fasta.h"
#include "align.h"
#include "version.h"

using namespace std;
using namespace pan;

Color::Modifier RED(Color::FG_RED);
Color::Modifier GREEN(Color::FG_GREEN);
Color::Modifier BLUE(Color::FG_BLUE);
Color::Modifier MAGENTA(Color::FG_MAGENTA);
Color::Modifier WHITE(Color::FG_WHITE);
Color::Modifier DEF(Color::FG_DEFAULT);

#define MATRIX_HOMO_SUMMARY_VERSION "1.000"
#define DATE __DATE__
//"2017-12-18"

void print_usage()
{
    char buff[3000];
    const char *help_info = 
            "matrix_homo_summary - summary and visualize two homologous matrix files\n"
            "======================================================================\n"
            "\e[1mUSAGE:\e[0m\n"
            "\tmatrix_homo_summary -in input_sam1,input_sam2/input_matrix1,input_matrix2 -chr chromosome_id1,chromosome_id2 \n"
            "\t               [-genome fasta | -msa sto_alignment] [-bins 40 -feature max -level auto -diagonal white -min_overhang 5 \n"
            "\t               -min_armlen 10 -strand +,+ -transform linear -dleft domain_file -dright domain_file -dchar \"*\" -dfill no] \n"
            "\e[1mHELP:\e[0m\n"
            
            "\t\e[1mInput Files: \e[0m\n"
            "\t-in: input .matrix file or .sam file  (default: guess by postfix)\n"
            "\t-chr: chr id, which must be consistant with sam file, genome fasta (default: no chr)\n"
            "\t-genome: fasta sequence for homologous alignment, mutual to -msa (default: no file)\n"
            "\t-msa: multiple alignment stockholm file, mutual to -genome (default: no file)\n\n"

            "\t\e[1mSam Filters: \e[0m\n"
            "\t-min_overhang: mininum overhang of duplex group, ignored when -file_type matrix (default: 5)\n"
            "\t-min_armlen: mininum arm length of each(left/right) arm, ignored when -file_type matrix (default: 10)\n"
            "\t-strand: +/- strand of reads(default: +)\n\n"

            "\t\e[1mPlot Parameters: \e[0m\n"
            "\t-bins: how many bins the raw matrix will reduce to (default: 40)\n"
            "\t-feature: max/mean/sum. how to calculate a value for each bin (default: max)\n"
            "\t-level: v1_1,v1_2,v1_3,V2_1,V2_2,V3_3 or auto. set the color cutoff of the value of each block \n"
            "\t        auto represent 65%%,80%%,90%% of all non-zero values (default: auto)\n"
            "\t-diagonal: white/nomal. white represent white diagonal; normal represent normal color (default: white)\n"
            "\t-transform: linear/log. how to transform the quantiles bar plot (default: linear)\n\n"
            
            "\t\e[1mDomain Parameters: \e[0m\n"
            "\t-dleft: a file to specify domains of left part (default: no file)\n"
            "\t-dright: a file to specify domains of right part (default: no file)\n"
            "\t-dchar: domain character (default: *)\n"
            "\t-dfill: fill domain (default: no)\n\n"

            "\e[1mVERSION:\e[0m\n\t%s\n"
            "\e[1mLIB VERSION:\e[0m\n\t%s\n"
            "\e[1mCOMPILE DATE:\e[0m\n\t%s\n"
            "\e[1mAUTHOR:\e[0m\n\t%s\n";

    sprintf(buff, help_info, MATRIX_HOMO_SUMMARY_VERSION, VERSION, DATE, "Li Pan");
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
    string input_file_1;
    string input_file_2;

    string left_domain_file;
    string right_domain_file;

    string genome;
    string sto_alignment;

    uLONG bins = 40;
    FEATURE feature = FEATURE_MAX;

    mutable Color_Control cc;

    bool auto_cc = true;

    uINT min_overhang = 5;
    uINT min_armlen = 10;
    
    char strand_1 = '+';
    char strand_2 = '+';

    string chr_id_1;
    string chr_id_2;

    char dchar = '*';
    bool dfill = false;

    string transform = "linear";
    bool diagonal = true;

    bool has_left_domain()const { return not left_domain_file.empty(); }
    bool has_right_domain()const { return not right_domain_file.empty(); }

    operator bool()
    {
        if(input_file_1.empty() or input_file_2.empty() or bins <= 5 or chr_id_1.empty() or chr_id_2.empty())
            return false;
        if(not cc.valid())
            return false;
        if(genome != "" and sto_alignment != "")
        {
            std::cerr << " -genome and -msa cannot be specified  " << endl;
            return false;
        }
        if(genome == "" and sto_alignment == "")
        {
            std::cerr << " Please specify -genome or -msa " << endl;
            return false;
        }
        return true;
        //return not input_file.empty() and bins > 5 and cc.valid();
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
            }else if(not strcmp(argv[i]+1, "genome"))
            {
                has_next(argc, i);
                param.genome = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "msa"))
            {
                has_next(argc, i);
                param.sto_alignment = argv[i+1];
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
                auto items = split(argv[i+1], ',');
                if(items.size() != 2)
                {
                    cerr << RED << "FATAL ERROR: bad -in option: " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                param.input_file_1 = items[0];
                param.input_file_2 = items[1];
                i++;
            }else if(not strcmp(argv[i]+1, "level"))
            {
                has_next(argc, i);
                if( not strcmp(argv[i+1], "auto") )
                {
                    param.auto_cc = true;
                }else{
                    auto levels = split(argv[i+1], ',');
                    if(levels.size() != 6)
                    {
                        cerr << RED << "FATAL ERROR: bad -level option: " << argv[i+1] << DEF << endl;
                        exit(-1);
                    }
                    param.cc.level_dn_1 = stod(levels[0]);
                    param.cc.level_dn_2 = stod(levels[1]);
                    param.cc.level_dn_3 = stod(levels[2]);

                    param.cc.level_up_1 = stod(levels[3]);
                    param.cc.level_up_2 = stod(levels[4]);
                    param.cc.level_up_3 = stod(levels[5]);

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
                auto items = split(argv[i+1], ',');
                if(items.size() != 2)
                {
                    cerr << RED << "FATAL ERROR: bad -chr option: " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                param.chr_id_1 = items[0];
                param.chr_id_2 = items[1];
                i++;
            }else if(not strcmp(argv[i]+1, "dchar"))
            {
                has_next(argc, i);
                param.dchar = argv[i+1][0];
                i++;
            }else if(not strcmp(argv[i]+1, "dright"))
            {
                has_next(argc, i);
                param.right_domain_file = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "dleft"))
            {
                has_next(argc, i);
                param.left_domain_file = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "strand"))
            {
                has_next(argc, i);
                auto items = split(argv[i+1], ',');
                if(items.size() != 2)
                {
                    cerr << RED << "FATAL ERROR: bad -strand option: " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                param.strand_1 = items[0] == "+" ? '+' : '-';
                param.strand_2 = items[0] == "+" ? '+' : '-';
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
void auto_color_control(const Statistic_Basic &sb_up, const Statistic_Basic &sb_dn, Color_Control &cc)
{
    cc.level_up_1 = sb_up.quantiles.at(13);
    cc.level_up_2 = sb_up.quantiles.at(16);
    cc.level_up_3 = sb_up.quantiles.at(19);

    cc.level_dn_1 = sb_dn.quantiles.at(13);
    cc.level_dn_2 = sb_dn.quantiles.at(16);
    cc.level_dn_3 = sb_dn.quantiles.at(19);
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


void print_statistics_figure(const Statistic_Basic &sb, const Param &param)
{
    using namespace std::placeholders;

    cout << DEF;
    // print title
    const uLONG bar_plot_width = 40;
    char buff[2000];
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

    cout << DEF << string(title_indent, ' ') << buff << "\n";
    if(param.transform == "linear")
        print_statistics(sb, param.cc, indent);
    else if(param.transform == "log")
        print_statistics(sb, param.cc, indent, std::bind(TRANSFORM::log, 2, _1));
}

template<typename T>
void expand_matrix(const Matrix<T> &raw_matrix, Matrix<T> &expanded_matrix, const Multi_Align &alignment, const string &chromosome_name)
{
    using namespace std::placeholders;
    auto coor_covert = std::bind(&Multi_Align::raw_coor_to_align_coor, _1, chromosome_name, _2);
    auto matrix_size = raw_matrix.size();

    uLONG align_length = alignment.length();
    init_matrix(expanded_matrix, align_length);

    uLONGArray raw_to_align_map;
    for(uLONG idx=0; idx<matrix_size; idx++)
        raw_to_align_map.push_back( coor_covert(alignment, idx)-1 );

    for(uLONG idx=0; idx<matrix_size; idx++)
    {
        uLONG x_covert = raw_to_align_map[idx];
        for(uLONG idy=idx; idy<matrix_size; idy++)
        {
            uLONG y_covert = raw_to_align_map[idy];
            expanded_matrix[x_covert][y_covert] = raw_matrix[idx][idy];
            expanded_matrix[y_covert][x_covert] = raw_matrix[idy][idx];
        }
    }
}

void expand_domain_region(RegionArray &domains, const Multi_Align &alignment, const string &chromosome_name)
{
    using namespace std::placeholders;
    auto coor_covert = std::bind(&Multi_Align::raw_coor_to_align_coor, _1, chromosome_name, _2);

    for(uLONG idx=0; idx<domains.size();idx++)
    {
        domains[idx].first = coor_covert(alignment, domains[idx].first-1);
        domains[idx].second = coor_covert(alignment, domains[idx].second-1);
    }
}

template<typename T>
void combine_matrix(   const Matrix<T> &matrix_lower, const Matrix<T> &matrix_upper, Matrix<T> &c_matrix )
{
    if(matrix_lower.size() != matrix_upper.size())
        throw Unexpected_Error("Different Size of 2 input ");

    auto matrix_size = matrix_lower.size();

    c_matrix = matrix_lower;

    for(uLONG idx=0; idx<matrix_size; idx++)
        for(uLONG idy=idx+1; idy<matrix_size; idy++)
            // upper part
            c_matrix[idx][idy] = matrix_upper[idx][idy]; //rm_2.at(idx).at(idy).count;
}


void mask_gap_alignment_region( Matrix<double> &matrix,
                                const Multi_Align &alignment, 
                                const string &chromosome_name,
                                double mask_index=0.25)
{
    using namespace std::placeholders;
    auto matrix_size = matrix.size();

    if(not alignment.has(chromosome_name))
        throw Unexpected_Error("Chromosome "+chromosome_name+" is not in alignment file");

    auto coor_covert = std::bind(&Multi_Align::raw_coor_to_align_coor, _1, chromosome_name, _2);

    uLONG align_length = alignment.length();
    uLONG raw_length = alignment.get_sto_record(chromosome_name).seq_length;
    double step_size = 1.0 * align_length / matrix_size;

    uLONGArray countArray(matrix_size);
    for(uLONG x=0; x<raw_length; x++)
    {
        uLONG x_convert = coor_covert(alignment, x) - 1;
        x_convert /= step_size;
        countArray.at(x_convert)++;
    }

    uLONG count_each_bin = 1.0*align_length / matrix_size;
    for(uLONG x=0; x<countArray.size(); x++)
    {
        if(countArray[x] < mask_index * count_each_bin)
        {
            for(uLONG idx=0; idx<matrix_size; idx++)
                matrix[x][idx] = matrix[idx][x] = -1;
        }
    }
}

void matrix_homo_summary(const Param &param)
{
    using namespace std::placeholders;
    Matrix<double> matrix_1, comp_matrix_1, expand_matrix_1,  matrix_2, comp_matrix_2, expand_matrix_2, combined_matrix;
    vector<Duplex_Hang> duplex_array_1, duplex_array_2;
    uLONG chr_len_1, chr_len_2;
    RegionArray upper_regions, lower_regions, c_upper_regions, c_lower_regions;

    clog << "Start to load " << param.input_file_1 << "..." << endl;
    switch(FILE_FORMAT::guess_file_type(param.input_file_1))
    {
        case FILE_FORMAT::MATRIX_FILE:
            read_matrix(param.input_file_1, matrix_1);
            break;
        case FILE_FORMAT::SAM_FILE:
            chr_len_1 = get_chromosome_hang( param.input_file_1, duplex_array_1, param.min_overhang, param.min_armlen, param.chr_id_1, param.strand_1);
            fill_sym_matrix(matrix_1, duplex_array_1, chr_len_1);
            duplex_array_1.clear();
            break;
        default:
            cerr << RED << "FATAL ERROR: unknown file format" << DEF << endl;
            exit(-1);
    }

    clog << "Start to load " << param.input_file_2 << "..." << endl;
    switch(FILE_FORMAT::guess_file_type(param.input_file_2))
    {
        case FILE_FORMAT::MATRIX_FILE:
            read_matrix(param.input_file_2, matrix_2);
            break;
        case FILE_FORMAT::SAM_FILE:
            chr_len_2 = get_chromosome_hang( param.input_file_2, duplex_array_2, param.min_overhang, param.min_armlen, param.chr_id_2, param.strand_2);
            fill_sym_matrix(matrix_2, duplex_array_2, chr_len_2);
            duplex_array_2.clear();
            break;
        default:
            cerr << RED << "FATAL ERROR: unknown file format" << DEF << endl;
            exit(-1);
    }

    Multi_Align *p_ma = nullptr;

    if(param.genome != "")
    {
        clog << "Start to load " << param.genome << "..." << endl;
        Fasta genome(param.genome);
        if(not genome.has_chr(param.chr_id_1) or not genome.has_chr(param.chr_id_2))
        {
            cerr << RED << "FATAL ERROR: chr not in fasta file" << DEF << endl;
            exit(-1);
        }

        MapStringString align_seqs;
        if(param.chr_id_1 == param.chr_id_2)
        {
            align_seqs[param.chr_id_1] = genome.get_chr_seq(param.chr_id_2);
            align_seqs[param.chr_id_1+"_2"] = genome.get_chr_seq(param.chr_id_2);
        }else{
            clog << "Start to global_align_2_seq" << "..." << endl;
            global_align_2_seq(genome.get_chr_seq(param.chr_id_1), genome.get_chr_seq(param.chr_id_2), align_seqs[param.chr_id_1], align_seqs[param.chr_id_2], 2, 1);
        }
        clog << "Start to Multi_Align" << "..." << endl;
        
        p_ma =new Multi_Align(align_seqs);
    }else{
        p_ma =new Multi_Align(param.sto_alignment);
    }

    clog << "Start to expand_matrix 1" << "..." << endl;
    expand_matrix(matrix_1, expand_matrix_1, *p_ma, param.chr_id_1);
    clog << "Start to expand_matrix 2" << "..." << endl;
    expand_matrix(matrix_2, expand_matrix_2, *p_ma, param.chr_id_2);

    clog << "Start to compress_matrix 1" << "..." << endl;
    compress_matrix(expand_matrix_1, comp_matrix_1, param.bins, param.feature);
    double mask_index = 0.4;
    mask_gap_alignment_region( comp_matrix_1, *p_ma, param.chr_id_1, mask_index );
    clog << "Start to compress_matrix 2" << "..." << endl;
    compress_matrix(expand_matrix_2, comp_matrix_2, param.bins, param.feature);
    mask_gap_alignment_region( comp_matrix_2, *p_ma, param.chr_id_2, mask_index );

    clog << "Start to combine_matrix" << "..." << endl;
    combine_matrix(comp_matrix_1, comp_matrix_2, combined_matrix);
    //clog << combined_matrix << endl;

    Statistic_Basic sb_1, sb_2;
    clog << "Start to auto_statistic" << "..." << endl;
    auto_statistic(matrix_1, sb_1);
    auto_statistic(matrix_2, sb_2);

    //Color_Control cc;
    if(param.auto_cc)
    {
        auto_color_control(sb_2, sb_1, param.cc);
    }

    if( param.has_right_domain() )
    {
        read_domain_file(param.right_domain_file, upper_regions);
        expand_domain_region(upper_regions, *p_ma, param.chr_id_2);
        compress_regions(upper_regions, c_upper_regions, expand_matrix_2.size(), param.bins);
        //cout << c_upper_regions << endl;
    }
    if( param.has_left_domain() )
    {
        read_domain_file(param.left_domain_file, lower_regions);
        expand_domain_region(lower_regions, *p_ma, param.chr_id_1);
        compress_regions(lower_regions, c_lower_regions, expand_matrix_1.size(), param.bins);
    }

    print_heatmap(combined_matrix, param.cc, param.diagonal, 0, "\u25A0", STD_Color::WHITE, c_upper_regions, c_lower_regions, param.dchar, param.dfill);

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

    // upper part
    print_statistics_figure(sb_2, param);
    // lower part
    print_statistics_figure(sb_1, param);
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

    matrix_homo_summary(param);

    return 0;
}




