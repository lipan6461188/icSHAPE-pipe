
/*

    A Programe to plot PARIS heatmap with 2 homologous species

*/

#include "paris.h"
#include "paris_plot.h"
#include "param.h"
#include "fasta.h"

#include <QGuiApplication>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <cstring>
#include <iomanip>      // std::setprecision

using namespace std;
using namespace pan;

#define PARIS_HOMO_HEATMAP_VERSION "1.101"
#define DATE "2017-12-06"

Color::Modifier RED(Color::FG_RED);
Color::Modifier DEF(Color::FG_DEFAULT);
Color::Modifier YELLOW(Color::FG_YELLOW);

/*
bool calc_quantitle(const Matrix<double> &matrix, 
            double &value_1, 
            double &value_2, 
            double &value_3, 
            double &value_4, 
            const double ratio_1, 
            const double ratio_2);
*/

bool calc_quantitle(const Matrix<double> &matrix, 
                double &value_1, 
                double &value_2, 
                const double ratio_1, 
                const double ratio_2);

void print_usage()
{
    char buff[4000];
    const char *help_info = 
            "paris_homo_heatmap - plot PARIS heatmap with 2 homologous species \n"
            "=========================================================================\n"
            "\e[1mUSAGE:\e[0m\n"
            "\tparis_homo_heatmap -in input_sam_1,input_sam_2 -chr chromosome_id_1,chromosome_id_2 -align sto_file        \n"
            "\t              -out heatmap.pdf [ -min_overhang 5 -min_armlen 10 -color #C44E52 -gap_color #D1CFCF -min 5,5 \n"
            "\t              -max 30,30 -auto no -bins 200 -bin_image_size 3 -strand +,+ -h] \n"
            "\e[1mHELP:\e[0m\n"
            "\t-in: sam files of 2 chromosomes (default: no)\n"
            "\t-strand: +/- strand of reads(default: +)\n"
            "\t-chr: chr ids(the first chr data wiil be in the left panel), which must be consistant with the corresponding sam file (default: no)\n"
            "\t-align: sequence alignment sto file (default: no)\n"
            "\t-out: output pdf file (default: no)\n\n"

            "\t-min_overhang: mininum overhang of duplex group, ignored when -file_type matrix (default: 5)\n"
            "\t-min_armlen: mininum arm length of each(left/right) arm, ignored when -file_type matrix (default: 10)\n\n"

            "\t-color: heatmap color (default: #C44E52) \n"
            "\t-gap_color: heatmap color of gap region (default: #D1CFCF) \n"
            "\t-min: min reads per bin, will be ignored if -auto yes (default: 5) \n"
            "\t-max: max reads per bin, will be ignored if -auto yes (default: 30) \n"
            "\t-auto: yes/no. 5%% of all sorted non-zero value as min value, 85%% of all sorted non-zero value as max value (default: no) \n"
            "\t-bin_image_size: bin image size (default: 3) \n\n"

            "\t-h: show this information\n\n"

            "\e[1mVERSION:\e[0m\n\t%s\n"
            "\e[1mLIB VERSION:\e[0m\n\t%s\n"
            "\e[1mDATE:\e[0m\n\t%s\n"
            "\e[1mAUTHOR:\e[0m\n\t%s\n";
    sprintf(buff, help_info, PARIS_HOMO_HEATMAP_VERSION, VERSION, DATE, "Li Pan");
    cout << buff << endl;
}

struct Param
{
    enum FILE_TYPE {SAM_FILE, MATRIX_FILE};

    string input_sam_1;
    string input_sam_2;

    string chr_id_1;
    string chr_id_2;

    string output_pdf;

    string align_file;

    uINT min_overhang = 5;
    uINT min_armlen = 10;
    
    string color = "#C44E52";
    string gap_color = "#D1CFCF";

    double min_rc_1 = 5;
    double min_rc_2 = 5;
    double max_rc_1 = 30;
    double max_rc_2 = 30;
    bool auto_rc = false;
    
    uINT bins = 200;

    char strand_1 = '+';
    char strand_2 = '+';

    double bin_image_size = 3;

    operator bool()const { 
        
        if(input_sam_1.empty() or input_sam_2.empty() or chr_id_1.empty() or chr_id_2.empty() or output_pdf.empty() or align_file.empty())
            return false;

        if(min_rc_1 > max_rc_1 or min_rc_2 > max_rc_2)
            return false;

        return true;
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

    string structure;

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
            }else if(not strcmp(argv[i]+1, "in"))
            {
                has_next(argc, i);
                
                StringArray input_files = split(argv[i+1], ',');
                if(input_files.size() != 2)
                {
                    cerr << RED << "FATAL Error: Bad option value of -in " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                param.input_sam_1 = input_files[0];
                param.input_sam_2 = input_files[1];

                i++;
            }else if(not strcmp(argv[i]+1, "align"))
            {
                has_next(argc, i);
                param.align_file = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "chr"))
            {
                has_next(argc, i);
                
                StringArray chromosomes = split(argv[i+1], ',');
                if(chromosomes.size() != 2)
                {
                    cerr << RED << "FATAL Error: Bad option value of -chr " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                param.chr_id_1 = chromosomes[0];
                param.chr_id_2 = chromosomes[1];

                i++;
            }else if(not strcmp(argv[i]+1, "h"))
            {
                print_usage();
                exit(0);
            }else if(not strcmp(argv[i]+1, "out"))
            {
                has_next(argc, i);
                param.output_pdf = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "color"))
            {
                has_next(argc, i);
                param.color = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "gap_color"))
            {
                has_next(argc, i);
                param.color = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "min"))
            {
                has_next(argc, i);

                StringArray min_cutoffs = split(argv[i+1], ',');
                if(min_cutoffs.size() != 2)
                {
                    cerr << RED << "FATAL Error: Bad option value of -min " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                param.min_rc_1 = stod(min_cutoffs[0]);
                param.min_rc_2 = stod(min_cutoffs[1]);

                i++;
            }else if(not strcmp(argv[i]+1, "max"))
            {
                has_next(argc, i);
                
                StringArray max_cutoffs = split(argv[i+1], ',');
                if(max_cutoffs.size() != 2)
                {
                    cerr << RED << "FATAL Error: Bad option value of -max " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                param.max_rc_1 = stod(max_cutoffs[0]);
                param.max_rc_2 = stod(max_cutoffs[1]);

                i++;
            }else if(not strcmp(argv[i]+1, "auto"))
            {
                has_next(argc, i);
                if( not strcmp(argv[i+1], "yes") )
                    param.auto_rc = true;
                else if( not strcmp(argv[i+1], "yes") )
                    param.auto_rc = true;
                else
                {
                    cerr << RED << "FATAL Error: unknown option value of -auto " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                i++;
            }else if(not strcmp(argv[i]+1, "bins"))
            {
                has_next(argc, i);
                param.bins = stoul(string(argv[i+1]));
                i++;
            }else if(not strcmp(argv[i]+1, "bin_image_size"))
            {
                has_next(argc, i);
                param.bin_image_size = stod(string(argv[i+1]));
                i++;
            }else if(not strcmp(argv[i]+1, "strand"))
            {
                has_next(argc, i);

                StringArray strands = split(argv[i+1], ',');
                if(strands.size() != 2)
                {
                    cerr << RED << "FATAL Error: Bad option value of -max " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                param.strand_1 = ( strands[0] == "+" ? '+' : '-' );
                param.strand_2 = ( strands[1] == "+" ? '+' : '-' );

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


void paris_homo_heatmap(const Param &param)
{
    vector<Duplex_Hang> duplex_array_1, duplex_array_2;
    uLONG chr_len_1, chr_len_2;
    Matrix<double> matrix_1, matrix_2;
    Heatmap_Matrix heatmap_1, heatmap_2, c_heatmap;

    clog << "Start to load alignment..." << endl;
    Multi_Align align(param.align_file);

    clog << "Start to get_chromosome_hang..." << endl;
    chr_len_1 = get_chromosome_hang( param.input_sam_1, duplex_array_1, param.min_overhang, param.min_armlen, param.chr_id_1, param.strand_1);
    chr_len_2 = get_chromosome_hang( param.input_sam_2, duplex_array_2, param.min_overhang, param.min_armlen, param.chr_id_2, param.strand_2);

    clog << "Start to fill_sym_matrix_with_aligned_dh..." << endl;
    fill_sym_matrix_with_aligned_dh(   matrix_1, duplex_array_1,  align,  param.chr_id_1, param.bins );
    fill_sym_matrix_with_aligned_dh(   matrix_2, duplex_array_2,  align,  param.chr_id_2, param.bins );

    if(param.auto_rc)
    {
        clog << "Start to calc_quantitle..." << endl;
        if(not calc_quantitle(matrix_1, const_cast<Param &>(param).min_rc_1, const_cast<Param &>(param).max_rc_1, 0.05, 0.85) or
            not calc_quantitle(matrix_2, const_cast<Param &>(param).min_rc_2, const_cast<Param &>(param).max_rc_2, 0.05, 0.85))
        {
            cerr << RED << "FATAL Error: too less value to quantile, please use -min and -max to replace -auto " << DEF << endl;
            exit(-1);
        }
    }

    clog << "Start to matrix_to_heatmap..." << endl;
    matrix_to_heatmap<double>(matrix_1, heatmap_1, param.min_rc_1, param.max_rc_1, param.color);
    matrix_to_heatmap<double>(matrix_2, heatmap_2, param.min_rc_2, param.max_rc_2, param.color);

    clog << "Start to dark_gap_alignment_region..." << endl;
    dark_gap_alignment_region( heatmap_1, align, param.chr_id_1, param.gap_color);
    dark_gap_alignment_region( heatmap_2, align, param.chr_id_2, param.gap_color);

    clog << "Start to combine_heatmap..." << endl;
    combine_heatmap( heatmap_1, heatmap_2, c_heatmap );

    clog << "Start to plot_intergrated_heatmap..." << endl;
    plot_intergrated_heatmap( c_heatmap, param.output_pdf, param.bin_image_size);
}

int main(int argc, char *argv[])
{

    //push_params(argc, argv, "-platform");
    //push_params(argc, argv, "offscreen");

    QGuiApplication app(argc, argv);

    Param param = read_param(argc, argv);
    if(not param)
    {
        cerr << RED << "Parameter is not valid" << DEF << endl;
        print_usage();
        exit(-1);
    }

    paris_homo_heatmap(param);

    return 0;
}

/*
bool calc_quantitle(const Matrix<double> &matrix, 
            double &value_1, 
            double &value_2, 
            double &value_3, 
            double &value_4, 
            const double ratio_1, 
            const double ratio_2)
{
    DoubleArray non_zero_1, non_zero_2;
    for(uLONG i=0; i<matrix.size(); i++)
        for(uLONG j=0; j<matrix.size(); j++)
            if(matrix[i][j] > 0)
                if(i > j)
                    non_zero_1.push_back(matrix[i][j]);
                else
                    non_zero_2.push_back(matrix[i][j]);
    if( non_zero_1.size() < 20 or non_zero_2.size() < 20 )
        return false;
    sort(non_zero_1.begin(), non_zero_1.end());
    sort(non_zero_2.begin(), non_zero_2.end());
    value_1 = non_zero_1[ non_zero_1.size() * ratio_1 ];
    value_2 = non_zero_1[ non_zero_1.size() * ratio_2 ];
    value_3 = non_zero_2[ non_zero_1.size() * ratio_1 ];
    value_4 = non_zero_2[ non_zero_1.size() * ratio_2 ];
    return true;
}
*/

bool calc_quantitle(const Matrix<double> &matrix, double &value_1, double &value_2, const double ratio_1, const double ratio_2)
{
    DoubleArray non_zero;
    for(uLONG i=0; i<matrix.size(); i++)
        for(uLONG j=0; j<matrix.size(); j++)
            if(matrix[i][j] > 0)
                non_zero.push_back(matrix[i][j]);
    if( non_zero.size() < 20 )
        return false;
    sort(non_zero.begin(), non_zero.end());
    value_1 = non_zero[ non_zero.size() * ratio_1 ];
    value_2 = non_zero[ non_zero.size() * ratio_2 ];
    return true;
}













