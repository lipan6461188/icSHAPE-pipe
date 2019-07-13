
/*

    A Programe to plot a heatmap of PARIS data from SAM file or matrix file

    Examples:
        ./paris_heatmap -in 59.sam -chr KU501215.1 -bin_image_size 1 -max 35 -min 0 -bins 0 -regionFile /tmp/long_range/input.txt
        ./paris_heatmap -in 59.sam -chr KU501215.1 -out_pdf /tmp/59.pdf -bin_image_size 1 -max 35 -min 0 -region 591,667,2881,3002 -bins 0 -regionFile /tmp/long_range/input.txt

*/


#include "paris.h"
#include "paris_plot.h"
#include "param.h"
#include "fasta.h"
#include "sstructure.h"
#include "fold.h"

#include <QGuiApplication>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <cstring>
#include <iomanip>      // std::setprecision

using namespace std;
using namespace pan;

#define PARIS_HEATMAP_VERSION "1.101"
#define DATE "2017-11-23"

Color::Modifier RED(Color::FG_RED);
Color::Modifier DEF(Color::FG_DEFAULT);
Color::Modifier YELLOW(Color::FG_YELLOW);

void build_base_pairs(const string &dot_bracket, 
                    Heatmap_Matrix &heatmap,
                    const uLONG start_1, const uLONG end_1,
                    const uLONG start_2, const uLONG end_2);

void build_cross_linkings(const string &dot_bracket,
                        Heatmap_Matrix &heatmap,
                        const uLONG start_1, const uLONG end_1,
                        const uLONG start_2, const uLONG end_2,
                        const string &chromosome);

void clear_base_pairs(Heatmap_Matrix &heatmap);

bool calc_quantitle(const Matrix<double> &matrix, 
                double &value_1, 
                double &value_2, 
                const double ratio_1, 
                const double ratio_2);

void print_regionFile_format()
{
    const char *regionFile_format = 
        "\e[1mregionFile format:\e[0m\n"
        "\tstart_1,end_1,start_2,end_2 out_pdf_file_name out_matrix_file_name NULL\n"
        "\tstart_1,end_1,start_2,end_2 out_pdf_file_name NULL ..((..))..\n"
        "\tstart_1,end_1,start_2,end_2 out_pdf_file_name NULL ..((..))..\n"
        "\t.............\n";
    cout << regionFile_format << endl;
}

void print_usage()
{
    char buff[4000];
    const char *help_info = 
            "paris_heatmap - plot a heatmap of PARIS data from SAM file or matrix file\n"
            "=========================================================================\n"
            "\e[1mUSAGE:\e[0m\n"
            "\tparis_heatmap -in input_sam/input_matrix -chr chromosome_id [-file_type sam -out_pdf output.pdf -out_matrix output.txt \n"
            "\t              -min_overhang 5 -min_armlen 10 -color #C44E52 -min 10 -max 100 -auto no -bins 200 -bin_image_size 3 \n"
            "\t              -region 0,0,0,0 -regionFile NULL -faFile NULL -shapeFile NULL -structure NULL -strand + -local yes -h] \n"
            "\e[1mHELP:\e[0m\n"
            "\t-in: sam file or matrix file(default: no)\n"
            "\t-file_type: input file type -- sam or matrix (default: sam) \n"
            "\t-strand: +/- strand of reads(default: +)\n"
            "\t-chr: chr id, which must be consistant with sam file, genome fasta and shape file (default: no)\n"
            "\t-out_pdf: output pdf file (default: no)\n"
            "\t-out_matrix: output matrix file of process region (default: no)\n\n"

            "\t-min_overhang: mininum overhang of duplex group, ignored when -file_type matrix (default: 5)\n"
            "\t-min_armlen: mininum arm length of each(left/right) arm, ignored when -file_type matrix (default: 10)\n"

            "\t-color: heatmap color (default: #C44E52) \n"
            "\t-min: min reads per bin, will be ignored if -auto yes (default: 10) \n"
            "\t-max: max reads per bin, will be ignored if -auto yes (default: 100) \n"
            "\t-auto: yes/no. 5%% of all sorted non-zero value as min value, 85%% of all sorted non-zero value as max value (default: no) \n"
            "\t-bin_image_size: bin image size (default: 3) \n\n"

            "\t-bins: bin numbers, should little than 1/2 of raw size, -bins 0 means no compress (default: 200) \n"
            "\t-region: 1-based coordination of chromosome to plot, format: start_1,end_1,start_2,end_2 (default: 0,0,0,0) \n"
            "\t-regionFile: Batch process, a file contain many regions to plot (default: no file) \n\n"
            
            "\t-faFile: a genome fasta file, if bins is not 0, it is ignored (default: no file) \n"
            "\t-shapeFile: a icSHAPE file, if bins is not 0, it is ignored (default: no file) \n"
            "\t-structure: auto or dot-bracket secondary structure of interaction regions(cat them directly), if bins is not 0, it is ignored (default: no) \n"
            "\t-local: yes/no, permit local base pairing (default: yes) \n\n"
            
            "\t-h: show this information and regionFile format\n\n"

            "\e[1mVERSION:\e[0m\n\t%s\n"
            "\e[1mLIB VERSION:\e[0m\n\t%s\n"
            "\e[1mDATE:\e[0m\n\t%s\n"
            "\e[1mAUTHOR:\e[0m\n\t%s\n";
    sprintf(buff, help_info, PARIS_HEATMAP_VERSION, VERSION, DATE, "Li Pan");
    cout << buff << endl;
}

struct Region_Plot
{
    uLONG start_1 = 0;
    uLONG end_1 = 0;
    uLONG start_2 = 0;
    uLONG end_2 = 0;

    string out_pdf;
    string out_matrix;

    string structure;
};

struct Param
{
    enum FILE_TYPE {SAM_FILE, MATRIX_FILE};

    string input_file;
    string out_matrix;
    string out_pdf;
    string chr_id;

    uINT min_overhang = 5;
    uINT min_armlen = 10;
    string color = "#C44E52";
    double min_rc = 10;
    double max_rc = 100;
    bool auto_rc = false;
    uINT bins = 200;

    string faFile;
    string shapeFile;

    char strand = '+';

    bool local = true;

    mutable vector<Region_Plot> plot_regions;

    FILE_TYPE file_type = SAM_FILE;

    double bin_image_size = 3;

    operator bool()const { 
        
        if(file_type == SAM_FILE)
        {
            if( chr_id.empty() or input_file.empty() )
                return false;
        }
        else if(file_type == MATRIX_FILE)
        {
            if ( input_file.empty() )
                return false;
            if ( not faFile.empty() or not shapeFile.empty() )
                return false;
        }
        else
        {
            return false;
        }

        if(plot_regions.size() == 0 and out_pdf.empty())
            return false;

        for(auto iter=plot_regions.begin(); iter!=plot_regions.end(); iter++)
        {
            if(iter->out_pdf.empty() and iter->out_matrix.empty())
                return false;
            if(bins != 0)
            {
                iter->structure.clear();
            }
        }

        if(std::any_of(plot_regions.cbegin(), plot_regions.cend(), [](const Region_Plot &region){ return region.structure == "auto" ? true : false; }))
            if(faFile.empty())
                return false;

        return true;
    }

    bool will_plot_region() const
    {
        return plot_regions.size() != 0;
    }
};

void read_regionFile(const string &file_name, Param &param)
{
    ifstream IN(file_name, ifstream::in);
    if(not IN)
    {
        cerr << "FATAL Error: " << file_name << " is unreadable" << endl;
        exit(-1);
    }

    string cur_line;
    while(getline(IN, cur_line))
    {
        if(cur_line.empty())
            continue;
        StringArray items;
        split(cur_line, items);
        if(items.size() != 4)
        {
            cerr << "FATAL Error: " << file_name << " format error" << endl;
            exit(0);
        }
        Region_Plot plot_region;

        StringArray region = split(items[0], ',');
        
        plot_region.start_1 = stoul(region[0]);
        plot_region.end_1 = stoul(region[1]);
        plot_region.start_2 = stoul(region[2]);
        plot_region.end_2 = stoul(region[3]);

        if(items[1] != "NULL")
            plot_region.out_pdf = items[1];
        if(items[2] != "NULL")
            plot_region.out_matrix = items[2];
        if(items[3] != "NULL")
            plot_region.structure = items[3];
        param.plot_regions.push_back(plot_region);
    }

    IN.close();
}


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
    int region_idx = -1;

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
                param.input_file = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "chr"))
            {
                has_next(argc, i);
                param.chr_id = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "h"))
            {
                print_usage();
                print_regionFile_format();
                exit(0);
            }else if(not strcmp(argv[i]+1, "out_pdf"))
            {
                has_next(argc, i);
                param.out_pdf = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "out_matrix"))
            {
                has_next(argc, i);
                param.out_matrix = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "color"))
            {
                has_next(argc, i);
                param.color = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "min"))
            {
                has_next(argc, i);
                param.min_rc = stod(string(argv[i+1]));
                i++;
            }else if(not strcmp(argv[i]+1, "max"))
            {
                has_next(argc, i);
                param.max_rc = stod(string(argv[i+1]));
                i++;
            }else if(not strcmp(argv[i]+1, "auto"))
            {
                has_next(argc, i);
                if( not strcmp(argv[i+1], "yes") )
                    param.auto_rc = true;
                else if( not strcmp(argv[i+1], "no") )
                    param.auto_rc = false;
                else
                {
                    cerr << RED << "FATAL Error: unknown option value of -auto " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                i++;
            }else if(not strcmp(argv[i]+1, "local"))
            {
                has_next(argc, i);
                if( not strcmp(argv[i+1], "yes") )
                    param.local = true;
                else if( not strcmp(argv[i+1], "no") )
                    param.local = false;
                else
                {
                    cerr << RED << "FATAL Error: unknown option value of -local " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                i++;
            }else if(not strcmp(argv[i]+1, "bins"))
            {
                has_next(argc, i);
                param.bins = stoul(string(argv[i+1]));
                i++;
            }else if(not strcmp(argv[i]+1, "faFile"))
            {
                has_next(argc, i);
                param.faFile = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "shapeFile"))
            {
                has_next(argc, i);
                param.shapeFile = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "structure"))
            {
                has_next(argc, i);
                structure = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "regionFile"))
            {
                has_next(argc, i);
                read_regionFile(argv[i+1], param);
                //param.regionFile = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "bin_image_size"))
            {
                has_next(argc, i);
                param.bin_image_size = stod(string(argv[i+1]));
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
            }else if(not strcmp(argv[i]+1, "region"))
            {
                has_next(argc, i);
                StringArray region = split(argv[i+1], ',');
                if(region.size() != 4)
                {
                    cerr << RED << "FATAL ERROR: -region " << argv[i+1] << " format error" << DEF << endl;
                    print_usage();
                    exit(-1);
                }
                Region_Plot plot_region;

                plot_region.start_1 = stoul(region[0]);
                plot_region.end_1 = stoul(region[1]);
                plot_region.start_2 = stoul(region[2]);
                plot_region.end_2 = stoul(region[3]);

                param.plot_regions.push_back(plot_region);
                region_idx = param.plot_regions.size() - 1;
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


    //string faFile, shapeFile, structure;

    if(region_idx >= 0)
    {
        if( not param.out_pdf.empty() )
            param.plot_regions[region_idx].out_pdf = param.out_pdf;
        if( not param.out_matrix.empty() )
            param.plot_regions[region_idx].out_matrix = param.out_matrix;
        if( not structure.empty() )
            param.plot_regions[region_idx].structure = structure;
    }

    if(param.bins != 0)
    {
        param.faFile.clear();
        param.shapeFile.clear();
    }

    return param;
}

void save_matrix(const string &file_name, 
                const Matrix<double> &matrix, 
                uLONG start_1=1, uLONG end_1=-1UL,
                uLONG start_2=1, uLONG end_2=-1UL)
{
    using size_type = Matrix<double>::size_type;
    uLONG matrix_size = matrix.size();

    /* Check Range */
    if(start_1 > matrix_size or start_1 == 0)
        throw runtime_error("Out of Matrix Range");

    if(end_1 >= matrix_size)
        end_1 = matrix_size;

    if(start_1 > end_1)
        throw runtime_error("Bad Start End Range");


    if(start_2 > matrix_size or start_2 == 0)
        throw runtime_error("Out of Matrix Range");

    if(end_2 >= matrix_size)
        end_2 = matrix_size;

    if(start_2 > end_2)
        throw runtime_error("Bad Start End Range");

    //clog << "Start to write matrix..." << endl;
    ofstream OUT(file_name, ofstream::out);
    if(not OUT)
    {
        cerr << RED << "FATAL Error: " + file_name + " cannot be writebale" << DEF << endl;
        exit(-1);
    }

    OUT << std::setprecision(7) << std::fixed;
    for(size_type idy=start_2-1;idy<end_2;idy++)
    {
        for(size_type idx=start_1-1;idx<end_1;idx++)
        {
            OUT << matrix[idx][idy];
            if(idx != end_1-1)
                OUT << "\t";
        }
        OUT << "\n";
    }

    OUT.close();
}

void paris_heatmap(const Param &param)
{
    vector<Duplex_Hang> duplex_array;
    uLONG chr_len;
    Matrix<double> matrix, comp_matrix;
    Heatmap_Matrix heatmap;

    if(param.file_type == Param::SAM_FILE)
    {
        try{
            chr_len = get_chromosome_hang( param.input_file, duplex_array, param.min_overhang, param.min_armlen, param.chr_id, param.strand);
            fill_sym_matrix(matrix, duplex_array, chr_len);
            duplex_array.clear();
        }catch(runtime_error e)
        {
            cerr << RED << "FATAL Error: " << e.what() << DEF << endl;
            print_usage();
            exit(-1);
        }
        
    }
    else{
        try{
            chr_len = read_matrix(param.input_file, matrix);
        }catch(runtime_error e)
        {
            cerr << RED << "FATAL Error: " << e.what() << DEF << endl;
            print_usage();
            exit(-1);
        }
    }

    if(param.bins != 0)
    {
        compress_matrix(matrix,
                        comp_matrix,
                        param.bins,
                        FEATURE_MEAN);
    }else{
        comp_matrix = matrix;
    }

    if( param.auto_rc  )
    {
        if(not calc_quantitle(comp_matrix, const_cast<Param &>(param).min_rc, const_cast<Param &>(param).max_rc, 0.05, 0.85))
        {
            cerr << RED << "FATAL Error: too less value to quantile, please use -min and -max to replace -auto " << DEF << endl;
            exit(-1);
        }
    }

    matrix_to_heatmap<double>(comp_matrix,
                            heatmap,
                            param.min_rc,     // min value
                            param.max_rc,     // max value
                            param.color);

    string genome_seq;
    Chr_Shape shape_array;
    if(not param.faFile.empty())
    {
        Fasta fasta(param.faFile);
        genome_seq = fasta.seq(param.chr_id);
    }
    if(not param.shapeFile.empty())
    {
        icSHAPE shape(param.shapeFile);
        shape_array = *shape.get_shape(param.chr_id);
    }

    if(param.will_plot_region())
    {
        for(auto iter=param.plot_regions.cbegin(); iter!=param.plot_regions.cend(); iter++)
        {

            clear_base_pairs(heatmap);

            if(not iter->structure.empty())
            {
                string dot_bracket;
                //clog << iter->structure << endl;
                if(iter->structure == "auto")
                {
                    const string seq_1 = genome_seq.substr(iter->start_1-1, iter->end_1-iter->start_1+1);
                    const string seq_2 = genome_seq.substr(iter->start_2-1, iter->end_2-iter->start_2+1);
                    clog << seq_1 << "III" << seq_2 << endl;
                    string structure = fold_two_seq(seq_1, seq_2, not param.local)->at(1);
                    clog << structure << endl;
                    structure = structure.substr(0, seq_1.size()) + structure.substr(seq_1.size()+3);
                    dot_bracket = structure;
                }else{
                    string dot_bracket = iter->structure;
                }
                build_base_pairs(dot_bracket, heatmap, iter->start_1, iter->end_1, iter->start_2, iter->end_2);
                build_cross_linkings(dot_bracket, heatmap, iter->start_1, iter->end_1, iter->start_2, iter->end_2, genome_seq);
            }

            if(not iter->out_pdf.empty())
                plot_intergrated_heatmap( heatmap, iter->out_pdf, param.bin_image_size, iter->start_1, iter->end_1, iter->start_2, iter->end_2, genome_seq, shape_array);
            if(not iter->out_matrix.empty())
                save_matrix(iter->out_matrix, comp_matrix, iter->start_1, iter->end_1, iter->start_2, iter->end_2);
        }
    }else{
        if(not param.out_pdf.empty())
            plot_intergrated_heatmap( heatmap, param.out_pdf, param.bin_image_size, 1, -1UL, 1, -1UL, genome_seq, shape_array);
        if(not param.out_matrix.empty())
            save_matrix(param.out_matrix, comp_matrix);
    }
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

    paris_heatmap(param);

    return 0;
}

void clear_base_pairs(Heatmap_Matrix &heatmap)
{
    for(auto iter=heatmap.matrix.begin(); iter!=heatmap.matrix.end(); iter++)
        for(auto iter_iter=iter->begin(); iter_iter!=iter->end(); iter_iter++)
        {
            iter_iter->border = false;
        }
}

void build_base_pairs(const string &dot_bracket, 
                    Heatmap_Matrix &heatmap,
                    const uLONG start_1, const uLONG end_1,
                    const uLONG start_2, const uLONG end_2)
{
    uLONG plot_matrix_width = end_1 - start_1 + 1;
    uLONG plot_matrix_height = end_2 - start_2 + 1;

    if(dot_bracket.size() != plot_matrix_width+plot_matrix_height )
    {
        cerr << RED << "FATAL Error: dot_bracket length error" << DEF << endl;
        exit(-1);
    }

    SStructure structure(dot_bracket);
    BpArray bps = structure.get_base_pairs();

    for(auto iter=bps.cbegin(); iter!=bps.cend(); iter++)
    {
        uLONG left_index, right_index;
        if(iter->left <= plot_matrix_width)
            left_index = start_1 + iter->left - 1;
        else
            left_index = start_2 + (iter->left - plot_matrix_width) - 1;

        if(iter->right <= plot_matrix_width)
            right_index = start_1 + iter->right - 1;
        else
            right_index = start_2 + (iter->right - plot_matrix_width) - 1;

        heatmap.matrix[left_index-1][right_index-1].border = true;
        heatmap.matrix[left_index-1][right_index-1].border_color = QColor("#64b5cd");
        heatmap.matrix[left_index-1][right_index-1].border_width = 30;

        heatmap.matrix[right_index-1][left_index-1].border = true;
        heatmap.matrix[right_index-1][left_index-1].border_color = QColor("#64b5cd");
        heatmap.matrix[right_index-1][left_index-1].border_width = 30;
    }
}


void build_cross_linkings(const string &dot_bracket,
                        Heatmap_Matrix &heatmap,
                        const uLONG start_1, const uLONG end_1,
                        const uLONG start_2, const uLONG end_2,
                        const string &chromosome)
{
    const string &tethered_seq = chromosome.substr(start_1-1, end_1-start_1+1) + chromosome.substr(start_2-1, end_2-start_2+1);

    uLONG plot_matrix_width = end_1 - start_1 + 1;
    uLONG plot_matrix_height = end_2 - start_2 + 1;

    if(dot_bracket.size() != plot_matrix_width+plot_matrix_height )
    {
        cerr << RED << "FATAL Error: dot_bracket length error" << DEF << endl;
        exit(-1);
    }

    SStructure structure(tethered_seq, dot_bracket);
    const vector<Point> points = search_TT_cross_linking(structure);

    for(auto iter=points.cbegin(); iter!=points.cend(); iter++)
    {
        uLONG left_index, right_index;
        if(iter->first <= plot_matrix_width)
            left_index = start_1 + iter->first - 1;
        else
            left_index = start_2 + (iter->first - plot_matrix_width) - 1;

        if(iter->second <= plot_matrix_width)
            right_index = start_1 + iter->second - 1;
        else
            right_index = start_2 + (iter->second - plot_matrix_width) - 1;

        heatmap.matrix[left_index-1][right_index-1].border = true;
        heatmap.matrix[left_index-1][right_index-1].border_color = QColor("#8172b2");
        heatmap.matrix[left_index-1][right_index-1].border_width = 60;

        heatmap.matrix[right_index-1][left_index-1].border = true;
        heatmap.matrix[right_index-1][left_index-1].border_color = QColor("#8172b2");
        heatmap.matrix[right_index-1][left_index-1].border_width = 60;
    }
}

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





