
#include <iostream>
#include <param.h>
#include <string_split.h>
#include <exceptions.h>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <iomanip>

using namespace std;
using namespace pan;

#define null -999.0

Color::Modifier RED(Color::FG_RED);
Color::Modifier DEF(Color::FG_DEFAULT);
Color::Modifier YELLOW(Color::FG_YELLOW);

void print_usage()
{
    char buff[4000];
    const char *help_info = 
        "## --------------------------------------\n"
        "Calculate enrichment file using RT stop as foreground and base density as background\n\n"

        "Command:\n"
        " %s -o output_enrichment_file -f foreground_file -b background_file\n"
        "# what it is:\n"
        " -o     output enrichment file\n"
        " -f     target RT stop file (must be normalized)\n"
        " -b     background RT stop file (must be normalized)\n\n"

        "# more options:\n"
        " -e     enrich method, subtraction or dividing or complex (default: complex)\n"
        " -w     winsor method, defined as ( default: factor5:scaling1 )\n"
        " -x     subtraction factor (default: 0.25)\n"
        " -y     dividing factor (default: 10)\n\n"

        " -d     number of leading nucleotides to crop\n"
        " -l     number of tailing nucleotides to crop\n"
        " -g     input signals are in log rather than in normal space, then enrich method will be forced to subtraction \n\n";

    sprintf(buff, help_info, "calcEnrich");
    cout << buff << endl;
}

enum EnrichMethod { SUBSTRACTION, DIVIDING, COMPLEX };

struct Param
{
    string input_fg_file;
    string input_bg_file;
    string output_enrichment_file;

    EnrichMethod enrich_method = COMPLEX;
    double winsor_factor = 0.05;
    uINT winsor_scaling = 1;

    double sub_factor = 0.25;
    double div_factor = 10;

    bool log_op = false;

    operator bool() const {
        if( input_fg_file.empty() or input_bg_file.empty() or output_enrichment_file.empty() )
        {
            cerr << RED << "Error: please specify -o, -f, and -b" << DEF << endl;
            return false;
        }
        
        if(winsor_factor >= 0.5)
        {   
            cerr << RED << "Error: Winsor factor must be less than 0.5." << DEF << endl;
            return false;
        }
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

    for(size_t i=1; i<argc; i++)
    {
        if( argv[i][0] == '-' )
        {
            if(not strcmp(argv[i]+1, "o"))
            {
                has_next(argc, i);
                param.output_enrichment_file = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "f"))
            {
                has_next(argc, i);
                param.input_fg_file = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "b"))
            {
                has_next(argc, i);
                param.input_bg_file = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "e"))
            {
                has_next(argc, i);
                if(not strcmp(argv[i+1], "subtraction"))
                    param.enrich_method = SUBSTRACTION;
                else if(not strcmp(argv[i+1], "dividing"))
                    param.enrich_method = DIVIDING;
                else
                    param.enrich_method = COMPLEX;
                i++;
            }else if(not strcmp(argv[i]+1, "w"))
            {
                has_next(argc, i);

                //factor5:scaling10
                string winsor_method(argv[i+1]);
                StringArray data = split(winsor_method, ':');
                for(string part: data)
                {
                    if(part.substr(0, 6) == "factor")
                        param.winsor_factor = stoul(part.substr(6))/100.0;
                    else if(part.substr(0, 7) == "scaling")
                        param.winsor_scaling = stoul(part.substr(7));
                    else{
                        cerr << RED << "FATAL Error: invalid -w winsor_method " << endl;
                        exit(-1);
                    }
                }
                i++;
            }else if(not strcmp(argv[i]+1, "x"))
            {
                has_next(argc, i);
                param.sub_factor = stod(argv[i+1]);
                i++;
            }else if(not strcmp(argv[i]+1, "y"))
            {
                has_next(argc, i);
                param.div_factor = stod(argv[i+1]);
                i++;
            }else if(not strcmp(argv[i]+1, "g"))
            {
                //has_next(argc, i);
                param.log_op = true;
                //i++;
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

/*
void ave_rpkm_string(const string &rpkm_string, double &rpkm)
{
    if(rpkm_string.find(',') != string::npos)
    {
        StringArray rpkm_list;
        split(rpkm_string, ',', rpkm_list);
        
        double Sum;
        for(const string &it: rpkm_list)
            Sum += stod(it);
        
        rpkm = Sum / rpkm_list.size();
    }else{
        rpkm = stod(rpkm_string);
    }
}
*/

void readSignal(const string &signalFile, MapStringuLONG &trans_len, MapStringString &trans_rpkm, 
    MapStringDouble &trans_scalingFactor_bd, MapStringT<DoubleArray> &trans_baseDensity, 
    MapStringDouble &trans_scalingFactor_rt, MapStringT<DoubleArray> &trans_RTstop)
{
    ifstream SG(signalFile, ifstream::in);
    string line;

    uLONG lineCount = 0;
    while(getline(SG, line))
    {
        if(line[0]=='#') continue;

        ++lineCount;
        if(lineCount%10000==0)
            cerr << "\tlines " << lineCount << endl;

        trim(line);
        StringArray data;
        split(line, '\t', data);
        string id(data[0]);
        uINT len( stoul(data[1]) );
        string type(data[2]);
        string rpkm(data[3]);
        double scalingFactor( stod(data[4]) );
        //DoubleArray baseDensities;
        
        trans_len[id] = len;
        trans_rpkm[id] = rpkm;

        if(type == "baseDensity")
        {
            //ave_rpkm_string(rpkm, trans_rpkm[id]);
            trans_scalingFactor_bd[id] = scalingFactor;
            DoubleArray &baseDensities = trans_baseDensity[id];
            for(auto it=data.cbegin()+5; it!=data.cend(); it++) baseDensities.push_back( stod(*it) );
        }else if(type == "RTstop"){
            //ave_rpkm_string(rpkm, trans_rpkm[id]);
            trans_scalingFactor_rt[id] = scalingFactor;
            DoubleArray &rtStop = trans_RTstop[id];
            for(auto it=data.cbegin()+5; it!=data.cend(); it++) rtStop.push_back( stod(*it) );

        }
    }
}

void calcEnrich(const MapStringuLONG &trans_bg_len, const MapStringT<DoubleArray> &trans_bg_baseDensity, 
                const MapStringT<DoubleArray> &trans_bg_RTstop, const MapStringuLONG &trans_fg_len, 
                const MapStringT<DoubleArray> &trans_fg_baseDensity, const MapStringT<DoubleArray> &trans_fg_RTstop, 
                MapStringT<DoubleArray> &trans_enrichment, const Param &param )
{
    uLONG skipHead = 0;
    uLONG skipTail = 0;
    double sub_fac = param.sub_factor;
    double div_fac = param.div_factor;

    uLONG transCount = 0;
    for(auto it=trans_fg_RTstop.cbegin(); it!=trans_fg_RTstop.cend(); it++)
    {
        string trans(it->first);
        
        bool in_bg_len = (trans_bg_len.find(trans) != trans_bg_len.end());
        bool in_bg_bd = (trans_bg_baseDensity.find(trans) != trans_bg_baseDensity.end());
        bool in_bg_rt = (trans_bg_RTstop.find(trans) != trans_bg_RTstop.end());
        if(not in_bg_len or not in_bg_bd or not in_bg_rt)
        {
            cerr << RED << "Warning! transcript " << trans << " not defined in background file.\n" << DEF;
            continue;
        }
        bool in_fg_len = (trans_fg_len.find(trans) != trans_fg_len.end());
        if(not in_fg_len)
        {
            cerr << RED << "Warning! transcript " << trans << " not defined in foreground file.\n" << DEF;
            continue;
        }
        if(trans_fg_len.at(trans) != trans_bg_len.at(trans))
        {
            cerr << RED << "Warning! transcript " << trans << " is of different length in background and foreground files.\n" << DEF;
            continue;
        }

        ++transCount;
        if(transCount%10000==0)
            cerr << "\tprocess " << transCount << "..." << endl;

        DoubleArray &siganlEnrichment = trans_enrichment[trans];

        const DoubleArray &fg_rt_array = trans_fg_RTstop.at(trans);
        const DoubleArray &bg_rt_array = trans_bg_RTstop.at(trans);
        const DoubleArray &fg_bd_array = trans_fg_baseDensity.at(trans);
        const DoubleArray &bg_bd_array = trans_bg_baseDensity.at(trans);

        for(int i=0; i<skipHead; i++) siganlEnrichment.push_back(null);

        if(param.log_op or param.enrich_method==SUBSTRACTION)
        {
            for(int i=skipHead; i<trans_fg_len.at(trans)-skipTail; i++)
            {
                //char buff[100];
                //sprintf(buff, "%.3f", fg_rt_array[i]-sub_fac*bg_rt_array[i]);
                siganlEnrichment.push_back( fg_rt_array[i]-sub_fac*bg_rt_array[i] );
            }
        }else if(param.enrich_method==DIVIDING)
        {
            for(int i=skipHead; i<trans_fg_len.at(trans)-skipTail; i++)
                if(bg_bd_array[i] > 0)
                {
                    //char buff[100];
                    //sprintf(buff, "%.3f", div_fac * fg_bd_array[i] / bg_bd_array[i]);
                    siganlEnrichment.push_back( div_fac * fg_bd_array[i] / bg_bd_array[i] );
                }else{
                    siganlEnrichment.push_back(null);
                }

        }else{
            for(int i=skipHead; i<trans_fg_len.at(trans)-skipTail; i++)
                if(bg_bd_array[i] > 0)
                {
                    //char buff[100];
                    //sprintf(buff, "%.3f", div_fac * (fg_rt_array[i]-sub_fac*bg_rt_array[i]) / bg_bd_array[i]);
                    siganlEnrichment.push_back( div_fac * (fg_rt_array[i]-sub_fac*bg_rt_array[i]) / bg_bd_array[i] );
                }else{
                    siganlEnrichment.push_back(null);
                }
        }
        for(int i=trans_fg_len.at(trans)-skipTail; i<trans_fg_len.at(trans); i++) siganlEnrichment.push_back(null);
    }
}

void winsorWindow(  const DoubleArray &enrichment, const uINT &headTopSkip, const uINT &trimed_last, 
                    const double &winsor_factor, double &winsorLower, double &winsorUpper)
{
    /*
    winsor_factor /= 100;
    if(winsorFactor >= 0.5)
    {
        cerr << RED << "FATAL Error! Winsor factor must be less than 0.5.\n";
        exit(-1);
    }
    */
    
    //clog << "winsorWindow..." << endl;
    //clog << "enrichment length: " << enrichment.size() << endl;
    DoubleArray sortedTrimed_array;
    for(uINT i=headTopSkip; i<=trimed_last; i++)
    {
       // cout << i << ": " << enrichment[i] << endl;
        if(enrichment[i] != null)
            sortedTrimed_array.push_back(enrichment[i]);
    }

    //clog << "start to sort" << endl;
    sort(sortedTrimed_array.begin(), sortedTrimed_array.end());

    //for(double v: sortedTrimed_array)
    //    cout << v << " ";
    //cout << "\n\n";

    uLONG len = sortedTrimed_array.size();
    uLONG winsorLen = winsor_factor * len;

    winsorLower = sortedTrimed_array[winsorLen];
    winsorUpper = sortedTrimed_array[len-winsorLen-1];
}
 

void winsorization(const MapStringuLONG &trans_len, MapStringT<DoubleArray> &trans_signal, 
                    const uINT &headTopSkip, const uINT &tailToSkip,
                    const double &winsor_factor, const uINT &winsor_scaling)
{
    cerr << "Winsorization...\n\t" << currentDateTime() << endl;

    uLONG transCount = 0;
    for(auto it=trans_signal.begin(); it!=trans_signal.end(); it++)
    {

        ++transCount;
        if(transCount%10000==0)
            cerr << "\tprocess " << transCount << "...\n";

        string trans(it->first);
        uINT trimed_last = trans_len.at(trans) - tailToSkip - 1;

        //cerr << headTopSkip+1 << "\t" << trimed_last << endl;

        double winsorUpper, winsorLower;
        winsorWindow(it->second, headTopSkip+1, trimed_last, winsor_factor, winsorLower, winsorUpper);
        //clog << "winsorLower: " << winsorLower << "\t" << "winsorUpper: " << winsorUpper << endl;

        if(winsorLower < 0) winsorLower = 0;
        if(winsorLower >= winsorUpper)
        {
            cerr << "Not enough resolution! Skip transcript " << trans << endl;
            continue;
        }
        if(winsor_scaling != 0)
        {
            uLONG len = trans_len.at(it->first);
            for(uLONG i=0; i<len; i++)
            {
                if(it->second[i] != null)
                {
                    if(it->second[i] > winsorUpper)
                        it->second[i] = winsorUpper;
                    else if(it->second[i] < winsorLower)
                        it->second[i] = winsorLower;
                    it->second[i] = (it->second[i] - winsorLower)*winsor_scaling/(winsorUpper-winsorLower);
                }
            }
        }else{
            uLONG len = trans_len.at(it->first);
            for(uLONG i=0; i<len; i++)
            {
                if(it->second[i] != null)
                {
                    if(it->second[i] > winsorUpper)
                        it->second[i] = winsorUpper;
                    else if(it->second[i] < winsorLower)
                        it->second[i] = winsorLower;
                }
            }
        }
    }
}


void outputEnrichment(const string &outputFile, const MapStringuLONG &trans_len, const MapStringT<DoubleArray> &trans_enrichment, const MapStringString &trans_fg_rpkm, 
       const MapStringDouble &trans_fg_scalingFactor_rt, const MapStringT<DoubleArray> &trans_fg_RTstop, const MapStringString &trans_bg_rpkm, 
       const MapStringDouble &trans_bg_scalingFactor_bd, const MapStringDouble &trans_bg_scalingFactor_rt, const MapStringT<DoubleArray> &trans_bg_baseDensity, 
       const MapStringT<DoubleArray> &trans_bg_RTstop)
{
    cerr << "Output enrichment scores to file $outputFile...\n\t" << currentDateTime() << endl;

    ofstream OUT(outputFile, ofstream::out);
    if(not OUT)
    {
        cerr << RED << "FATAL Error: cannot write to " << outputFile << endl;
        exit(-1);
    }

    OUT << "#transcript\tlength\tenrichment score, start from position 1.\n";
    for(auto it=trans_enrichment.cbegin(); it!=trans_enrichment.cend(); it++)
    {
        string trans(it->first);

        OUT << trans << "\t" << trans_len.at(trans) << "\t" << trans_fg_rpkm.at(trans) << "," << trans_bg_rpkm.at(trans) << "\t" << 
        trans_fg_scalingFactor_rt.at(trans) << "," << trans_bg_scalingFactor_bd.at(trans) << "," << trans_bg_scalingFactor_rt.at(trans);

        const auto &enrichment = trans_enrichment.at(trans);
        const auto &fg_rt = trans_fg_RTstop.at(trans);
        const auto &bg_rt = trans_bg_RTstop.at(trans);
        const auto &bg_bd = trans_bg_baseDensity.at(trans);
        const double fg_sf_rt = trans_fg_scalingFactor_rt.at(trans);
        const double bg_sf_rt = trans_bg_scalingFactor_rt.at(trans);
        const double bg_sf_bd = trans_bg_scalingFactor_bd.at(trans);
        for(uLONG i=0; i<trans_len.at(trans); i++)
        {
            std::streamsize ss = OUT.precision();
            OUT << std::fixed << std::setprecision(3);

            if(enrichment[i] == null)
                OUT << "\tNULL";
            else
                OUT << "\t" << enrichment[i];
            OUT << "," << fg_sf_rt*fg_rt[i];
            OUT << "," << bg_sf_rt*bg_rt[i];
            OUT << "," << bg_sf_bd*bg_bd[i];

            OUT.precision(ss);
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
        print_usage();
        return -1;
    }

   // double pp = null;
   // if(null == pp)
   //     cerr << "Yean" << endl;

    MapStringuLONG trans_bg_len, trans_fg_len;
    MapStringString trans_bg_rpkm, trans_fg_rpkm;
    MapStringDouble trans_bg_scalingFactor_bd, trans_fg_scalingFactor_bd;
    MapStringDouble trans_bg_scalingFactor_rt, trans_fg_scalingFactor_rt;
    MapStringT<DoubleArray> trans_bg_baseDensity, trans_fg_baseDensity;
    MapStringT<DoubleArray> trans_bg_RTstop, trans_fg_RTstop;
    MapStringT<DoubleArray> trans_enrichment;
    //MapStringDouble trans_bg_avgRTstop, trans_fg_avgRTstop;

    cerr << "Start to read " << param.input_bg_file << endl;
    readSignal(param.input_bg_file, trans_bg_len, trans_bg_rpkm, trans_bg_scalingFactor_bd, 
            trans_bg_baseDensity, trans_bg_scalingFactor_rt, trans_bg_RTstop);
    
    cerr << "Start to read " << param.input_fg_file << endl;
    readSignal(param.input_fg_file, trans_fg_len, trans_fg_rpkm, trans_fg_scalingFactor_bd, 
            trans_fg_baseDensity, trans_fg_scalingFactor_rt, trans_fg_RTstop);

    cerr << "Start to calculate enrichment " << endl;
    calcEnrich(trans_bg_len, trans_bg_baseDensity, trans_bg_RTstop, 
               trans_fg_len, trans_fg_baseDensity, trans_fg_RTstop, 
               trans_enrichment, param );

    cerr << "Start to calculate winsorization " << endl;
    uLONG headToSkip = 5;  uLONG tailToSkip = 32;
    winsorization(trans_fg_len, trans_enrichment, headToSkip, tailToSkip,
                    param.winsor_factor, param.winsor_scaling);

    cerr << "Start to calculate outputEnrichment " << endl;
    outputEnrichment(param.output_enrichment_file, trans_fg_len, trans_enrichment, trans_fg_rpkm, 
       trans_fg_scalingFactor_rt, trans_fg_RTstop, trans_bg_rpkm, 
       trans_bg_scalingFactor_bd, trans_bg_scalingFactor_rt, trans_bg_baseDensity, 
       trans_bg_RTstop);

    return 0;
}

