
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
        "calculate RT stops from sam file\n\n"

        "Command:\n"
        "%s -i input_sam_file -o output_RTstop_file -r transcript_rpkm_file\n\n"

        "# what it is:\n"
        " -i     input sam file\n"
        " -o     output RTstop file (with base density information)\n"
        " -r     rpkm file\n\n"

        "# more options:\n"
        " -c     cutoff of RPKM\n\n";

    sprintf(buff, help_info, "calcRT");
    cout << buff << endl;
}

struct Param
{
    string input_file;
    string output_file;
    string rpkm_file;

    double rpkm_cutoff = 5.0;

    operator bool() const 
    {
        if( input_file.empty() or output_file.empty() or rpkm_file.empty() )
        {
            cerr << RED << "Error: please specify -i, -o and -r" << DEF << endl;
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
            }else if(not strcmp(argv[i]+1, "r"))
            {
                has_next(argc, i);
                param.rpkm_file = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "c"))
            {
                has_next(argc, i);
                param.rpkm_cutoff = stod(argv[i+1]);
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

void readRPKM(const string &rpkmFile, const double &minLoad, MapStringuLONG &trans_len, MapStringDouble &trans_rpkm)
{
    cerr << "Read transcript abundance information from file " << rpkmFile << "...\n\t" << currentDateTime() << endl;
    ifstream IN(rpkmFile, ifstream::in);
    if(not IN)
    {
        cerr << RED << "FATAL Error: cannot read " << rpkmFile << DEF << endl;
        exit(-1);
    }

    string line;
    while(getline(IN, line))
    {
        if(line[0] == '#') continue;

        StringArray data;
        split(line, '\t', data);

        string trans(data[0]);
        uLONG len = stoul(data[1]);
        double uniqRead = stod(data[2]);
        double multiRead = stod(data[3]);
        double rpkm = stod(data[4]);

        if(rpkm < minLoad) continue;

        trans_len[trans] = len;
        trans_rpkm[trans] = rpkm;
    }

    IN.close();
}

void calcBaseDensity(const string &inputSamFile, MapStringT<uLONGArray> &baseDensity, MapStringT<uLONGArray> &RTstop,
                    const MapStringuLONG &trans_len, const MapStringDouble &trans_rpkm)
{
    cerr << "Calculate base density from file " << inputSamFile <<"...\n\t" << currentDateTime() << endl;

    ifstream IN(inputSamFile, ifstream::in);
    if(not IN)
    {
        cerr << RED << "FATAL Error: cannot read " << inputSamFile << "..." << DEF << endl;
        exit(-1);
    }

    string line;
    string readID;
    uLONG lineCount = 0;
    MapStringuLONG hitID_start, hitID_end;
    uLONG hitCount = 0;

    while(getline(IN, line))
    {
        if(line[0] == '@') continue;

        ++lineCount;
        //cout << lineCount << endl;
        if(lineCount % 1000000 == 0)
            cerr << "\tlines " << lineCount << endl;

        StringArray data;
        split(line, '\t', data);
        string read(data[0]);
        string tag(data[1]);
        string hit(data[2]);
        uLONG pos = stoul(data[3]);
        uLONG tlen = stoul(data[8]);
        string seq(data[9]);

        if(tag == "99" or tag == "355")
        {
            if(data[6] == "=")
                continue;
        }else if(tag == "0" or tag == "256")
        {
            tlen = seq.size();
            if(not tlen)
                continue;
        }else{
            continue;
        }

        if(read != readID)
        {
            if(not readID.empty())
                for(auto it=hitID_start.cbegin(); it!=hitID_start.cend(); it++)
                {
                    string trans(it->first);
                    //if(trans_len.find(trans) == trans_len.end()) continue;
                    for(uLONG i=it->second; i<hitID_end.at(trans); i++)
                        baseDensity[trans][i] += 1.0/hitCount;
                    RTstop[trans][it->second-1] += 1.0/hitCount;
                }
            hitCount = 0;
            hitID_start.clear();
            hitID_end.clear();
        }

        readID = read;
        hitCount++;
        if(trans_len.find(hit) != trans_len.end())
        {
            hitID_start[hit] = pos;
            hitID_end[hit] = pos+tlen;
            if(trans_len.at(hit) < hitID_end[hit])
                hitID_end[hit] = trans_len.at(hit) + 1;
        }
    }

    if(not readID.empty())
        for(auto it=hitID_start.cbegin(); it!=hitID_start.cend(); it++)
        {
            string trans(it->first);
            //if(trans_len.find(trans) == trans_len.end()) continue;
            for(uLONG i=it->second; i<hitID_end.at(trans); i++)
                baseDensity[trans][i] += 1.0/hitCount;
            RTstop[trans][it->second-1] += 1.0/hitCount;
        }

    IN.close();
}

void output_baseDensity(const string &outputFile, const MapStringuLONG &trans_len, const MapStringDouble &trans_rpkm,
        const MapStringT<uLONGArray> &baseDensity, const MapStringT<uLONGArray> &RTstop)
{
    cerr << RED << "Output base density to file " << outputFile << "...\n\t" << currentDateTime() << endl;

    ofstream OUT(outputFile, ofstream::out);
    if(not OUT)
    {
        cerr << RED << "FATAL Error: cannot write to " << outputFile << DEF << endl;
        exit(-1);
    }
    OUT << "#transcript\tbase frequency, start from position 0.\n";
    for(auto it=baseDensity.cbegin(); it!=baseDensity.cend(); it++)
    {
        string trans(it->first);

        // BaseDensity
        OUT << trans << "\t" << trans_len.at(trans) << "\t" << trans_rpkm.at(trans);
        std::streamsize ss = OUT.precision();
        OUT << std::fixed << std::setprecision(3);
        for(uLONG i=0; i<=trans_len.at(trans); i++)
            OUT << "\t" << double(baseDensity.at(trans).at(i));
        OUT.precision(ss);

        OUT << "\n";

        // RTstop
        OUT << trans << "\t" << trans_len.at(trans) << "\t" << trans_rpkm.at(trans);
        ss = OUT.precision();
        OUT << std::fixed << std::setprecision(3);
        for(uLONG i=0; i<=trans_len.at(trans); i++)
            OUT << "\t" << double(RTstop.at(trans).at(i));
        OUT.precision(ss);

        OUT << "\n";
    }

    OUT.close();
}

void init_bd_rt(const MapStringuLONG &trans_len, MapStringT<uLONGArray> &baseDensity, MapStringT<uLONGArray> &RTstop)
{
    for(auto it=trans_len.cbegin(); it!=trans_len.cend(); it++)
    {
        baseDensity[it->first];
        baseDensity[it->first].resize(it->second+1);
        RTstop[it->first];
        RTstop[it->first].resize(it->second+1);
    }
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
    MapStringDouble trans_rpkm;
    MapStringT<uLONGArray> baseDensity;
    MapStringT<uLONGArray> RTstop;
    readRPKM(param.rpkm_file, param.rpkm_cutoff, trans_len, trans_rpkm);
    cerr << "Total number: " << trans_len.size() << endl;

    init_bd_rt(trans_len, baseDensity, RTstop);

    calcBaseDensity(param.input_file, baseDensity, RTstop, trans_len, trans_rpkm);
    output_baseDensity(param.output_file, trans_len, trans_rpkm, baseDensity, RTstop);
}












