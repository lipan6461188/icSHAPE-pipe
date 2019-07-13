
#include <iostream>
#include <param.h>
#include <string_split.h>
#include <exceptions.h>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <iomanip>
#include <numeric>

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
        "Calculate average signals\n\n"

        "Command:\n"
        " %s -i input_signal_file -o output_shape_file \n"
        "# what it is:\n"
        " -i     input signal file \n"
        " -o     input shape file \n\n"

        "# more options:\n"
        " -t     threshold of minimun coverage (default: 200)\n"
        " -T     threshold of average RT stop (default: 2)\n"
        " -s     head to skip (default: 5)\n"
        " -e     end to skip (default: 30)\n\n";

    sprintf(buff, help_info, "filterEnrich");
    cout << buff << endl;
}

//enum EnrichMethod { SUBSTRACTION, DIVIDING, COMPLEX };

struct Param
{
    string input_file;
    string output_file;

    uINT bd_cutoff = 200;
    double rt_cutoff = 2;
    uINT head_skip = 5;
    uINT tail_skip = 30;

    operator bool() const {
        if( input_file.empty() or output_file.empty() )
        {
            cerr << RED << "Error: please specify -i, -o" << DEF << endl;
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
            if(not strcmp(argv[i]+1, "i"))
            {
                has_next(argc, i);
                param.input_file = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "o"))
            {
                has_next(argc, i);
                param.output_file = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "t"))
            {
                has_next(argc, i);
                param.bd_cutoff = stoul(argv[i+1]);
                i++;
            }else if(not strcmp(argv[i]+1, "T"))
            {
                has_next(argc, i);
                param.rt_cutoff = stod(argv[i+1]);
                i++;
            }else if(not strcmp(argv[i]+1, "s"))
            {
                has_next(argc, i);
                param.head_skip = stoul(argv[i+1]);
                i++;
            }else if(not strcmp(argv[i]+1, "e"))
            {
                has_next(argc, i);
                param.tail_skip = stoul(argv[i+1]);
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

double ave_rpkm_string(const string &rpkm_string)
{
    double rpkm;
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

    return rpkm;
}

void readSignals(const Param &param)
{
    uINT bd_cutoff = param.bd_cutoff;
    double rt_cutoff = param.rt_cutoff;
    uINT head_skip = param.head_skip;
    uINT tail_skip = param.tail_skip;

    ofstream OUT(param.output_file, ofstream::out);
    if(not OUT)
    {
        cerr << RED << "FATAL Error: cannot write to " << param.output_file << DEF << endl;
        exit(-1);
    }

    ifstream IN(param.input_file, ifstream::in);
    if(not IN)
    {
        cerr << RED << "FATAL Error: cannot read " << param.input_file << DEF << endl;
        exit(-1);
    }

    cerr << "Start to read " << param.input_file << endl;

    string line;
    uLONG lineCount = 0;
    while(getline(IN, line))
    {
        if(line[0] == '#') continue;

        lineCount += 1;
        if(lineCount % 1000 == 0)
            cerr << "\t line " << lineCount << endl;

        StringArray data;
        split(line, '\t', data);

        string id(data.at(0));
        uLONG len = stoul(data.at(1));
        double rpkm = ave_rpkm_string(data.at(2));
        string scalingFactors(data.at(3));

        StringArray scores;
        DoubleArray fgRT, bgDB;
        for(uINT i=0; i<len; i++)
        {
            StringArray mini_data;
            split(data.at(i+4), ',', mini_data);
            scores.push_back(mini_data.at(0));
            fgRT.push_back( stod(mini_data.at(1)) );
            bgDB.push_back( stod(mini_data.at(3)) );
        }

        double coverage = std::accumulate(fgRT.cbegin(), fgRT.cend(), 0.0) / fgRT.size();
        if(coverage > rt_cutoff)
        {
            OUT << id << "\t" << len << "\t" << rpkm;

            for(uINT j=0; j<len; j++)
            {
                if(j>head_skip and j<len-tail_skip)
                    if(bgDB.at(j) >= bd_cutoff)
                        OUT << "\t" << scores.at(j);
                    else
                        OUT << "\tNULL";
                else
                    OUT << "\tNULL";
            }
            OUT << "\n";
        }
    }
    OUT.close();
}

int main(int argc, char *argv[])
{
    Param param = read_param(argc, argv);
    if(not param)
    {
        print_usage();
        exit(-1);
    }

    readSignals(param);
}




