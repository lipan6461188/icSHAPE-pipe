
//#include "Fold.h"
#include <RNA.h>
#include <src/ErrorChecker.h>
#include <src/ParseCommandLine.h>

#include "fold.h"
#include "pan_type.h"
#include "param.h"
#include "string_split.h"
#include "version.h"

#include <time.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
using namespace pan;

Color::Modifier RED(Color::FG_RED);
Color::Modifier GREEN(Color::FG_GREEN);
Color::Modifier BLUE(Color::FG_BLUE);
Color::Modifier MAGENTA(Color::FG_MAGENTA);
Color::Modifier WHITE(Color::FG_WHITE);
Color::Modifier DEF(Color::FG_DEFAULT);

#define FOLD_SAMPLE_VERSION "1.000"
#define DATE __DATE__
//"2017-12-18"

void print_usage()
{
    char buff[2000];
    const char *help_info = 
            "fold_param_sample - sample parameters of folding a RNA structure prediction\n"
            "===========================================================================\n"
            "\e[1mUSAGE:\e[0m\n"
            "\tfold_param_sample -in input_ct -out output_file -out_matrix out_matrix_prefix [-constraint constraint_file \n"
            "\t               -paris paris_matrix -shape shape_file -pi -2,0.2,21 -ps 0,-0.05,11 -si 0,-0.2,16 -ss 0,0.2,26 \n"
            "\t               -bp_max -1 -loop_max 30 -window -1 -verbose ] \n"
            "\e[1mHELP:\e[0m\n"

            "\t\e[1mInput Files: \e[0m\n"
            "\t-in: input a ct file with standard structure and sequence \n"
            "\t-out: output all information of each prediction\n"
            "\t-out_matrix: output sensitivity and ppv as a matrix\n\n"
            
            "\t\e[1mSam Filters: \e[0m\n"
            "\t-constraint: a constraint filr with same format as Fold\n"
            "\t-paris: a paris matrix\n"
            "\t-pi: start,step,count. specify paris interception params \n"
            "\t-ps: start,step,count. specify paris slop params \n"
            "\t-shape: a shape file\n"
            "\t-si: start,step,count. specify shape interception params \n"
            "\t-ss: start,step,count. specify shape slop params \n\n"

            "\t\e[1mOther Constraints: \e[0m\n"
            "\t-bp_max: max bp of each base pair (default: -1, unlimited)\n"
            "\t-loop_max: the maximun loop size (default: 30)\n"
            "\t-window: the window size for prediction (default: -1, auto) \n\n"

            "\t\e[1mOther: \e[0m\n"
            "\t-global: search best SHAPE/PARIS parameter (default: no) \n\n"
            "\t-verbose: show structure prediction process \n\n"

            "\e[1mVERSION:\e[0m\n\t%s\n"
            "\e[1mLIB VERSION:\e[0m\n\t%s\n"
            "\e[1mCOMPILE DATE:\e[0m\n\t%s\n"
            "\e[1mAUTHOR:\e[0m\n\t%s\n";

    sprintf(buff, help_info, FOLD_SAMPLE_VERSION, VERSION, DATE, "Li Pan");
    cout << buff << endl;
}

struct Param {

    string input_ct;            // The standard ct file.
    string output_matrix;       // The output matrix file
    string output_structure;    // The output structures       

    string constraint_file;     // The optional folding constraints file.
    string paris_file;          // The optional input bonus file.
    string shape_file;          // The optional SHAPE constraints file.

    /////////   paris

    double pi_start = -2.0;
    double pi_step = 0.2;
    uLONG pi_count = 21;

    double ps_start = 0.0;
    double ps_step = -0.05;
    uLONG ps_count = 11;

    /////////   icSHAPE

    double si_start = 0.0;
    double si_step = -0.2;
    double si_count = 16;

    double ss_start = 0.0;
    double ss_step = 0.2;
    double ss_count = 26;

    // mode
    bool search_global = false;

    // The maximum pairing distance.
    int bp_max = -1;

    // The maximum internal bulge loop size.
    int loop_max = 30;

    // The window size for calculation.
    int windowSize = -1;

    bool verbose = false;

    operator bool()
    {
        if(input_ct.empty())
            return false;
        if(shape_file.empty() and paris_file.empty())
            return false;
        return true;
    }

    void show() const;
};

void Param::show() const
{
    clog    << "input_ct: " << input_ct << "\n"
            << "output_matrix: " << output_matrix << "\n"
            << "output_structure: " << output_structure << "\n"
            << "constraint_file: " << constraint_file << "\n"
            << "paris_file: " << paris_file << "\n"
            << "shape_file: " << shape_file << "\n\n"
            
            << "pi_start: " << pi_start << "\n"
            << "pi_step: " << pi_step << "\n"
            << "pi_count: " << pi_count << "\n\n"
            
            << "ps_start: " << ps_start << "\n"
            << "ps_step: " << ps_step << "\n"
            << "ps_count: " << ps_count << "\n\n"
            
            << "si_start: " << si_start << "\n"
            << "si_step: " << si_step << "\n"
            << "si_count: " << si_count << "\n\n"
            
            << "ss_start: " << ss_start << "\n"
            << "ss_step: " << ss_step << "\n"
            << "ss_step: " << ss_step << "\n\n"

            << "bp_max: " << bp_max << "\n"
            << "loop_max: " << loop_max << "\n"
            << "windowSize: " << windowSize << "\n"
            << "search_global: " << search_global << "\n"
            << "verbose: " << verbose << "\n\n";
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
    for(int i=1; i<argc; i++)
    {
        if( argv[i][0] == '-' )
        {
            if(not strcmp(argv[i]+1, "in"))
            {
                has_next(argc, i);
                param.input_ct = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "out"))
            {
                has_next(argc, i);
                param.output_structure = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "out_matrix"))
            {
                has_next(argc, i);
                param.output_matrix = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "constraint"))
            {
                has_next(argc, i);
                param.constraint_file = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "paris"))
            {
                has_next(argc, i);
                param.paris_file = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "shape"))
            {
                has_next(argc, i);
                param.shape_file = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "pi"))
            {
                has_next(argc, i);
                auto items = split(argv[i+1], ',');
                if(items.size() != 3)
                {
                    cerr << RED << "FATAL Error: Bad option -pi " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                param.pi_start = stod(items[0]);
                param.pi_step = stod(items[1]);
                param.pi_count = stoul(items[2]);
                i++;
            }else if(not strcmp(argv[i]+1, "ps"))
            {
                has_next(argc, i);
                auto items = split(argv[i+1], ',');
                if(items.size() != 3)
                {
                    cerr << RED << "FATAL Error: Bad option -ps " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                param.ps_start = stod(items[0]);
                param.ps_step = stod(items[1]);
                param.ps_count = stoul(items[2]);
                i++;

            }else if(not strcmp(argv[i]+1, "si"))
            {
                has_next(argc, i);
                auto items = split(argv[i+1], ',');
                if(items.size() != 3)
                {
                    cerr << RED << "FATAL Error: Bad option -si " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                param.si_start = stod(items[0]);
                param.si_step = stod(items[1]);
                param.si_count = stoul(items[2]);
                i++;
            }else if(not strcmp(argv[i]+1, "ss"))
            {
                has_next(argc, i);
                auto items = split(argv[i+1], ',');
                if(items.size() != 3)
                {
                    cerr << RED << "FATAL Error: Bad option -ss " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                param.ss_start = stod(items[0]);
                param.ss_step = stod(items[1]);
                param.ss_count = stoul(items[2]);
                i++;
            }else if(not strcmp(argv[i]+1, "bp_max"))
            {
                has_next(argc, i);
                param.bp_max = stoi(argv[i+1]);
                i++;
            }else if(not strcmp(argv[i]+1, "loop_max"))
            {
                has_next(argc, i);
                param.loop_max = stoi(argv[i+1]);
                i++;
            }else if(not strcmp(argv[i]+1, "window"))
            {
                has_next(argc, i);
                param.windowSize = stoi(argv[i+1]);
                i++;
            }else if(not strcmp(argv[i]+1, "verbose"))
            {
                //has_next(argc, i);
                param.verbose = true;
                //i++;
            }else if(not strcmp(argv[i]+1, "global"))
            {
                //has_next(argc, i);
                param.search_global = true;
                //i++;
            }else{
                cerr << RED << "FATAL ERROR: unknown option: " << argv[i] << DEF << endl;
                print_usage();
                exit(-1);
            }
        }
        else{
            cerr << RED << "FATAL ERROR: unknown option: " << argv[i] << DEF << endl;
            print_usage();
            exit(-1);
        }
    }
    return param;
}

using Structure_Pairs = vector< pair<uLONG, uLONG> >;

uLONG common_bp_count(const Structure_Pairs &pairs_1, const Structure_Pairs &pairs_2)
{
    int common_pair_count( 0 );
    for(auto iter=pairs_1.cbegin(); iter!=pairs_1.cend(); iter++)
        if( std::binary_search(pairs_2.begin(), pairs_2.end(), *iter) )
            ++common_pair_count;

    return common_pair_count;
}

struct Prediction_Item
{
    double shape_intercept = -999;
    double shape_slop = -999;

    double paris_intercept = -999;
    double paris_slop = -999;

    uLONG test_bp_num = 0;
    uLONG standard_bp_num = 0;
    uLONG common_bp_num = 0;

    double sensitivity() const { return 1.0 * common_bp_num / test_bp_num; };
    double ppv() const { return 1.0 * common_bp_num / standard_bp_num; };
    double weighted_ave() const { return sqrt(this->sensitivity() * this->ppv()); };

    string structure;
};

void save_matrix(ofstream &OUT, const Rect<Prediction_Item> &rect)
{
    OUT << std::setprecision(3);

    for(uLONG idx=0; idx<rect.size(); idx++)
    {
        for(uLONG idy=0; idy<rect[idx].size(); idy++)
        {
            OUT << rect[idx][idy].sensitivity() << ";" << rect[idx][idy].ppv();
            if(idy != rect[idx].size()-1)
                OUT << "\t";
        }
        OUT << "\n";
    }
    OUT.close();
}

void save_index(ofstream &OUT, const Rect<Prediction_Item> &rect)
{
    for(uLONG idx=0; idx<rect.size(); idx++)
    {
        for(uLONG idy=0; idy<rect[idx].size(); idy++)
        {
            const auto &cur_cell = rect[idx][idy];
            if(cur_cell.paris_intercept == -999)
            {
                OUT << "si: " << cur_cell.shape_intercept << "\tss: " << cur_cell.shape_slop 
                    << "\t" << cur_cell.test_bp_num << "\t" << cur_cell.standard_bp_num 
                    << "\t" << cur_cell.common_bp_num << "\t" << cur_cell.sensitivity() 
                    << "\t" << cur_cell.ppv() << "\t" << cur_cell.structure << "\n";
            }else{
                OUT << "si: " << cur_cell.shape_intercept << "\tss: " << cur_cell.shape_slop << "\t"
                    << "pi: " << cur_cell.paris_intercept << "\tps: " << cur_cell.paris_slop
                    << "\t" << cur_cell.test_bp_num << "\t" << cur_cell.standard_bp_num 
                    << "\t" << cur_cell.common_bp_num << "\t" << cur_cell.sensitivity() 
                    << "\t" << cur_cell.ppv() << "\t" << cur_cell.structure << "\n";
            }
        }
    }
}

inline bool warning_on_error(const int error_code)
{
    if(error_code != 0)
    {
        cerr << RED << "FATAL Error: Error Occurred" << DEF << endl;
        //exit(-1);
    }
    return error_code != 0;
}

void exit_on_error(const int error_code)
{
    if(error_code != 0)
    {
        cerr << RED << "FATAL Error: Error Occurred" << DEF << endl;
        exit(-1);
    }
}




void sample_params(const Param &param) {

    int error = 0; const int percent(10); const bool quickfold(true); const string saveFile;
    const bool simple_iloops(true); const bool disablecoax(false); const bool allowUnkownBases(true);
    const bool skipThermo(true); const int maxStructures(20); const double temperature(310.15);

    //clog << "Initializing nucleic acids..." << endl;
    RNA strand( param.input_ct.c_str(), FILE_CT, DT_RNA, allowUnkownBases, skipThermo);
    Structure_Pairs standard_pairs = get_pairs(strand, 1);

    //strand.getStructure()->RemoveAllStructures();

    ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( &strand );
    error = checker->isErrorStatus();
    exit_on_error(error);

    //clog << "done" << endl;

    /* estimate window size */
    int windowSize = param.windowSize;
    if( windowSize == -1 ) {
        int length = strand.GetSequenceLength();
        windowSize =
            ( length > 1200 ) ? 20 :
            ( length > 800 ) ? 15 :
            ( length > 500 ) ? 11 :
            ( length > 300 ) ? 7 :
            ( length > 120 ) ? 5 :
            ( length > 50 ) ? 3 :
            2;
    }

    //ofstream OUT(param.output_structure);
    double best_si(param.si_start), best_ss(param.ss_start), current_max_ppv( 0 ), current_max_sensi( 0 ), best_weighted_ave( 0 );

    Rect<Prediction_Item> shape_rect;
    //clog << "init_rect" << endl;
    init_rect(shape_rect, param.ss_count, param.si_count);
    //clog << "done" << endl;
    /* Step 1: Fit the best SHAPE parameters */

    for( uLONG si_index=0; si_index<param.si_count; si_index++ )
    {

        // if no shape file. break
        if( param.shape_file.empty() ) break;

        double si = param.si_start + si_index * param.si_step;
        //clog << "outer" << endl;

        for( uLONG ss_index=0; ss_index<param.ss_count; ss_index++ )
        {
            //clog << "inner" << endl;
            double ss = param.ss_start + ss_index * param.ss_step;

            if(param.verbose)
                clog << "Select SHAPE Params: si: " << si << "; ss: " << ss << endl; 

            //clog << "Access shape_rect" << endl;
            shape_rect[ss_index][si_index].shape_intercept = si;
            shape_rect[ss_index][si_index].shape_slop = ss;
            //clog << "done" << endl;

            const bool IsRNA=true;
            //clog << "Create a new RNA" << endl;
            RNA *new_strand = new RNA(strand.GetSequence(), IsRNA);
            //clog << "done" << endl;

            int tempError = new_strand->SetTemperature( temperature );
            error = checker->isErrorStatus( tempError );
            if(warning_on_error(error)) continue;

           // clog << "Part 1" << endl;

            if( param.bp_max != -1 ) {
                int distError = new_strand->ForceMaximumPairingDistance( param.bp_max );
                error = checker->isErrorStatus( distError );
            }

            if(not param.constraint_file.empty())
            {
                int constraintError = new_strand->ReadConstraints( param.constraint_file.c_str() );
                error = checker->isErrorStatus( constraintError );
                if(warning_on_error(error)) continue;
            }

            int slopeSingle(0), interceptSingle(0);
            int shapeError = new_strand->ReadSHAPE( param.shape_file.c_str(), ss, si, slopeSingle, interceptSingle, RESTRAINT_SHAPE );
            error = checker->isErrorStatus( shapeError );
            if(warning_on_error(error)) continue;

            /*
            clog    << "\npercent: " << percent << "\n"
                    << "maxStructures: " << maxStructures << "\n"
                    << "windowSize: " << windowSize << "\n"
                    << "saveFile.c_str(): " << saveFile.c_str() << "\n"
                    << "maxLoop: " << param.loop_max << "\n"
                    << "quickfold: " << quickfold << "\n"
                    << "simple_iloops: " << simple_iloops << "\n"
                    << "disablecoax: " << disablecoax << "\n\n";
            */

            int mainCalcError = new_strand->FoldSingleStrand( percent, maxStructures, windowSize, saveFile.c_str(), param.loop_max, quickfold, simple_iloops, disablecoax);
            error = checker->isErrorStatus( mainCalcError );
            if(warning_on_error(error)) continue;

            if(error == 0) // no error
            {
                Structure_Pairs new_pairs = get_pairs(*new_strand, 1);
                shape_rect[ss_index][si_index].test_bp_num = new_pairs.size();
                shape_rect[ss_index][si_index].standard_bp_num = standard_pairs.size();
                shape_rect[ss_index][si_index].common_bp_num = common_bp_count(new_pairs, standard_pairs);
                shape_rect[ss_index][si_index].structure = get_dot_bracket(*new_strand)->at(0);
            }else{
                shape_rect[ss_index][si_index].test_bp_num = 0;
                shape_rect[ss_index][si_index].standard_bp_num = 0;
                shape_rect[ss_index][si_index].common_bp_num = 0;
                shape_rect[ss_index][si_index].structure = "";
            }

            if(shape_rect[ss_index][si_index].weighted_ave() > best_weighted_ave)
            {
                best_si = si;
                best_ss = ss;
                current_max_ppv = shape_rect[ss_index][si_index].ppv();
                current_max_sensi = shape_rect[ss_index][si_index].sensitivity();
                best_weighted_ave = shape_rect[ss_index][si_index].weighted_ave();
            }

            delete new_strand;
        }
    }

    clog << "Find the best SHAPE Fit: si=" << best_si << "; ss=" << best_ss << "; ppv=" << current_max_ppv << "; sensitivity=" << current_max_sensi << "; weighted_ave=" << best_weighted_ave << endl;

    /* Step 2: Fit the best PARIS parameters with the best SHAPE parameter */

    double best_pi(param.si_start), best_ps(param.ss_start);

    Rect<Prediction_Item> paris_rect;
    init_rect(paris_rect, param.ps_count, param.pi_count);

    for( uLONG pi_index=0; pi_index<param.pi_count; pi_index++ )
    {
        // if no paris file. break
        if( param.paris_file.empty() ) break;

        double pi = param.pi_start + pi_index * param.pi_step;

        for( uLONG ps_index=0; ps_index<param.ps_count; ps_index++ )
        {
            double ps = param.ps_start + ps_index * param.ps_step;

            if(param.verbose)
                clog << "Select PARIS Params: pi: " << pi << "; ps: " << ps << endl; 

            paris_rect[ps_index][pi_index].shape_intercept = best_si;
            paris_rect[ps_index][pi_index].shape_slop = best_ss;

            paris_rect[ps_index][pi_index].paris_intercept = pi;
            paris_rect[ps_index][pi_index].paris_slop = ps;

            const bool IsRNA=true;
            RNA *new_strand = new RNA(strand.GetSequence(), IsRNA);

            int tempError = new_strand->SetTemperature( temperature );
            error = checker->isErrorStatus( tempError );
            if(warning_on_error(error)) continue;

            if( param.bp_max != -1 ) {
                int distError = new_strand->ForceMaximumPairingDistance( param.bp_max );
                error = checker->isErrorStatus( distError );
            }

            if(not param.constraint_file.empty())
            {
                int constraintError = new_strand->ReadConstraints( param.constraint_file.c_str() );
                error = checker->isErrorStatus( constraintError );
                if(warning_on_error(error)) continue;
            }

            if(not param.shape_file.empty())
            {
                int slopeSingle(0), interceptSingle(0);
                int shapeError = new_strand->ReadSHAPE( param.shape_file.c_str(), best_ss, best_si, slopeSingle, interceptSingle, RESTRAINT_SHAPE );
                error = checker->isErrorStatus( shapeError );
                if(warning_on_error(error)) continue;
            }

            int parisError = new_strand->ReadExperimentalPairBonus( param.paris_file.c_str(), pi, ps );
            error = checker->isErrorStatus( parisError );
            if(warning_on_error(error)) continue;

            int mainCalcError = new_strand->FoldSingleStrand( percent, maxStructures, windowSize, saveFile.c_str(), param.loop_max, quickfold, simple_iloops, disablecoax);
            error = checker->isErrorStatus( mainCalcError );
            if(warning_on_error(error)) continue;
            //}
            if(error == 0) // no error
            {
                Structure_Pairs new_pairs = get_pairs(*new_strand, 1);
                paris_rect[ps_index][pi_index].test_bp_num = new_pairs.size();
                paris_rect[ps_index][pi_index].standard_bp_num = standard_pairs.size();
                paris_rect[ps_index][pi_index].common_bp_num = common_bp_count(new_pairs, standard_pairs);
                paris_rect[ps_index][pi_index].structure = get_dot_bracket(*new_strand)->at(0);
            }else{
                paris_rect[ps_index][pi_index].test_bp_num = 0;
                paris_rect[ps_index][pi_index].standard_bp_num = 0;
                paris_rect[ps_index][pi_index].common_bp_num = 0;
                paris_rect[ps_index][pi_index].structure = "";
            }

            if(paris_rect[ps_index][pi_index].weighted_ave() > best_weighted_ave)
            //if(paris_rect[ps_index][pi_index].ppv() > current_max_ppv)
            {
                best_pi = pi;
                best_ps = ps;
                current_max_ppv = paris_rect[ps_index][pi_index].ppv();
                current_max_sensi = paris_rect[ps_index][pi_index].sensitivity();
                best_weighted_ave = paris_rect[ps_index][pi_index].weighted_ave();
            }

            delete new_strand;
        }
    }

    delete checker;

    clog << "Find the best PARIS Fit: pi=" << best_pi << "; ps=" << best_ps << "; ppv=" << current_max_ppv << "; sensitivity=" << current_max_sensi << "; weighted_ave=" << best_weighted_ave << endl;


    ofstream OUT_LIST;
    if(not param.output_structure.empty())
        OUT_LIST.open(param.output_structure, ofstream::out);
    else
        OUT_LIST.open("sample_fold_param", ofstream::out);

    string matrix_prefix = param.output_matrix;
    if(matrix_prefix.empty())
        matrix_prefix = "sample_fold_param";

    if(not param.shape_file.empty())
    {
        ofstream OUT(matrix_prefix+"_shape.matrix", ofstream::out);
        save_matrix(OUT, shape_rect);
        save_index(OUT_LIST, shape_rect);
        OUT.close();
    }

    if(not param.paris_file.empty())
    {
        ofstream OUT(matrix_prefix+"_paris.matrix", ofstream::out);
        save_matrix(OUT, paris_rect);
        save_index(OUT_LIST, paris_rect);
        OUT.close();
    }

    OUT_LIST.close();
}


void sample_params_global(const Param &param) {

    int error = 0; const int percent(10); const bool quickfold(true); const string saveFile;
    const bool simple_iloops(true); const bool disablecoax(false); const bool allowUnkownBases(true);
    const bool skipThermo(true); const int maxStructures(20); const double temperature(310.15);

    //clog << "Initializing nucleic acids..." << endl;
    RNA strand( param.input_ct.c_str(), FILE_CT, DT_RNA, allowUnkownBases, skipThermo);
    Structure_Pairs standard_pairs = get_pairs(strand, 1);

    //strand.getStructure()->RemoveAllStructures();

    ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( &strand );
    error = checker->isErrorStatus();
    exit_on_error(error);

    //clog << "done" << endl;

    /* estimate window size */
    int windowSize = param.windowSize;
    if( windowSize == -1 ) {
        int length = strand.GetSequenceLength();
        windowSize =
            ( length > 1200 ) ? 20 :
            ( length > 800 ) ? 15 :
            ( length > 500 ) ? 11 :
            ( length > 300 ) ? 7 :
            ( length > 120 ) ? 5 :
            ( length > 50 ) ? 3 :
            2;
    }

    double best_si(param.si_start), best_ss(param.ss_start), current_max_ppv( 0 ), current_max_sensi( 0 ), best_weighted_ave( 0 );
    double best_pi(param.si_start), best_ps(param.ss_start);

    if( param.shape_file.empty() )
    {
        cerr << RED << "FATAL Error: Please specify -shape parameter in global mode" << DEF << endl;
        exit(-1);
    }
    if( param.paris_file.empty() )
    {
        cerr << RED << "FATAL Error: Please specify -paris parameter in global mode" << DEF << endl;
        exit(-1);
    }

    for( uLONG si_index=0; si_index<param.si_count; si_index++ )
    {
        // if no shape file. break
        
        double si = param.si_start + si_index * param.si_step;

        for( uLONG ss_index=0; ss_index<param.ss_count; ss_index++ )
        {
            double ss = param.ss_start + ss_index * param.ss_step;

            if(param.verbose)
                clog << "Select icSHAPE Params: si: " << si << "; ss: " << ss << endl; 

            for( uLONG pi_index=0; pi_index<param.pi_count; pi_index++ )
            {
                // if no paris file. break
                
                double pi = param.pi_start + pi_index * param.pi_step;

                for( uLONG ps_index=0; ps_index<param.ps_count; ps_index++ )
                {
                    double ps = param.ps_start + ps_index * param.ps_step;

                    if(param.verbose)
                        clog << "\tSelect PARIS Params: pi: " << pi << "; ps: " << ps << endl; 

                    const bool IsRNA=true;
                    RNA *new_strand = new RNA(strand.GetSequence(), IsRNA);

                    int tempError = new_strand->SetTemperature( temperature );
                    error = checker->isErrorStatus( tempError );
                    if(warning_on_error(error)) continue;

                    if( param.bp_max != -1 ) {
                        int distError = new_strand->ForceMaximumPairingDistance( param.bp_max );
                        error = checker->isErrorStatus( distError );
                    }

                    if(not param.constraint_file.empty())
                    {
                        int constraintError = new_strand->ReadConstraints( param.constraint_file.c_str() );
                        error = checker->isErrorStatus( constraintError );
                        if(warning_on_error(error)) continue;
                    }

                    if(not param.shape_file.empty())
                    {
                        int slopeSingle(0), interceptSingle(0);
                        int shapeError = new_strand->ReadSHAPE( param.shape_file.c_str(), ss, si, slopeSingle, interceptSingle, RESTRAINT_SHAPE );
                        error = checker->isErrorStatus( shapeError );
                        if(warning_on_error(error)) continue;
                    }

                    int parisError = new_strand->ReadExperimentalPairBonus( param.paris_file.c_str(), pi, ps );
                    error = checker->isErrorStatus( parisError );
                    if(warning_on_error(error)) continue;

                    int mainCalcError = new_strand->FoldSingleStrand( percent, maxStructures, windowSize, saveFile.c_str(), param.loop_max, quickfold, simple_iloops, disablecoax);
                    error = checker->isErrorStatus( mainCalcError );
                    if(warning_on_error(error)) continue;
                    //}
                    if(error == 0) // no error
                    {
                        Structure_Pairs new_pairs = get_pairs(*new_strand, 1);
                        Prediction_Item item;
                        item.test_bp_num = new_pairs.size();
                        item.standard_bp_num = standard_pairs.size();
                        item.common_bp_num = common_bp_count(new_pairs, standard_pairs);
                        item.structure = get_dot_bracket(*new_strand)->at(0);

                        if(item.weighted_ave() > best_weighted_ave)
                        //if(paris_rect[ps_index][pi_index].ppv() > current_max_ppv)
                        {
                            best_si = si;
                            best_ss = ss;
                            best_pi = pi;
                            best_ps = ps;
                            current_max_ppv = item.ppv();
                            current_max_sensi = item.sensitivity();
                            best_weighted_ave = item.weighted_ave();
                        }
                    }else{
                        // nothing happen
                    }

                    delete new_strand;
                }
            }
        }
    }

    delete checker;

    clog << "Find the best Fit: si=" << best_si << "; ss=" << best_ss << "; pi=" << best_pi << "; ps=" << best_ps << "; ppv=" << current_max_ppv << "; sensitivity=" << current_max_sensi << "; weighted_ave=" << best_weighted_ave << endl;

}


int main( int argc, char* argv[] ) {

    Param param = read_param(argc, argv);
    if(param.verbose)
        param.show();
    if(not param)
    {
        cerr << RED << "Parameter is not valid" << DEF << endl;
        print_usage();
        exit(-1);
    }

    if(param.search_global)
        sample_params_global(param);
    else
        sample_params(param);
    return 0;
}



