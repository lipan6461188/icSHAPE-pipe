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
        " %s -i input_normedRT -r output_rt_bedGraph -b output_bd_bedGraph \n\n";

    sprintf(buff, help_info, "normedRT2bedGraph");
    cout << buff << endl;
}


struct Param
{
    string input_file;
    string output_rt_file;
    string output_bd_file;

    operator bool() const {
        if( input_file.empty() or output_rt_file.empty() or output_bd_file.empty()  )
        {
            cerr << RED << "Error: please specify input and output file" << DEF << endl;
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
                param.input_file= argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "r"))
            {
                has_next(argc, i);
                param.output_rt_file = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "b"))
            {
                has_next(argc, i);
                param.output_bd_file = argv[i+1];
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

void writeBedGraph(const string &outRTFile, const string &outBDFile, MapStringDouble &trans_scalingFactor_bd, 
    MapStringT<DoubleArray> &trans_baseDensity, MapStringT<DoubleArray> &trans_RTstop)
{
    ofstream RT(outRTFile, ofstream::out);
    if(not RT)
    {
        cerr << RED << "FATAL Error: cannot read " << outRTFile << DEF << endl;
        exit(-1);
    }

    ofstream BD(outBDFile, ofstream::out);
    if(not BD)
    {
        cerr << RED << "FATAL Error: cannot read " << outBDFile << DEF << endl;
        exit(-1);
    }

    for(auto it=trans_RTstop.cbegin(); it!=trans_RTstop.cend(); it++)
    {
        string trans(it->first);

        const DoubleArray &RT_array = it->second;
        const DoubleArray &BD_array = trans_baseDensity.at(trans);
        const double BD_scaling = trans_scalingFactor_bd.at(trans);
        for(int i=0;i<RT_array.size();i++)
            RT << trans << "\t" << i << "\t" << i+1 << "\t" << RT_array.at(i) << "\n";
        for(int i=0;i<BD_array.size();i++)
            BD << trans << "\t" << i << "\t" << i+1 << "\t" << round(BD_scaling*BD_array.at(i)) << "\n";
    }

    RT.close();
    BD.close();
}

int main(int argc, char *argv[])
{
    Param param = read_param(argc, argv);
    if(not param)
    {
        print_usage();
        return -1;
    }

    MapStringuLONG trans_len;
    MapStringString trans_rpkm;
    MapStringDouble trans_scalingFactor_bd, trans_scalingFactor_rt;
    MapStringT<DoubleArray> trans_baseDensity, trans_RTstop;

    readSignal(param.input_file, trans_len, trans_rpkm,  trans_scalingFactor_bd, trans_baseDensity, 
        trans_scalingFactor_rt, trans_RTstop);

    writeBedGraph(param.output_rt_file, param.output_bd_file, trans_scalingFactor_bd, 
    trans_baseDensity, trans_RTstop);

    return 0;
}


