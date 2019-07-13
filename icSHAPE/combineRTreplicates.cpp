
#include <iostream>
#include <param.h>
#include <string_split.h>
#include <exceptions.h>
#include <fstream>

using namespace std;
using namespace pan;

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
            "%s -i input_rt_1:input_rt_2 -o output_combined_signal_file \n"
            " # what it is:\n"
            " -i     input signal files\n"
            " -o     output combined signal files\n\n"
            
            "# more options:\n"
            " -f     input format (normalized or count, default: count)\n";
    sprintf(buff, help_info, "combineRTreplicates");
    cout << buff << endl;
}

struct Param
{
    
    StringArray input_rt_list;
    string output_rt;

    string format = "count";

    operator bool() const { 
        if( input_rt_list.size() <= 1 or output_rt.empty() ) return false;
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
                param.input_rt_list = split(argv[i+1], ':');
                i++;
            }else if(not strcmp(argv[i]+1, "o"))
            {
                has_next(argc, i);
                param.output_rt = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "f"))
            {
                has_next(argc, i);
                param.format = argv[i+1];
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


void combineSignals(const string &rt_file, MapStringuLONG &trans_len, MapStringString &trans_rpkm, 
                    MapStringString &trans_scalingFactor, MapStringT<DoubleArray> &trans_rtstop, 
                    MapStringT<DoubleArray> &trans_bd, const string &format)
{
    cerr << "read signal from " << rt_file << "\n\t" << currentDateTime() << endl;

    string id, idOld, type, lineOld;
    uLONG len;
    string rpkm, scalingFactor;

    ifstream IN(rt_file, ifstream::in);
    if(not IN)
    {
        cerr << "FATAL Error: " << rt_file << " is not exists" << endl;
        exit(-1);
    }

    string line; uLONG lineCount(0);
    while(getline(IN, line))
    {
        if(line[0] == '#') continue;

        ++lineCount;
        if(lineCount % 10000 == 0)
            cerr << "  line: " << lineCount << endl;

        trim(line);
        StringArray data;
        split(line, '\t', data);

        DoubleArray rtstops, baseDensities;

        if(data[0] == idOld)
        {
            if(format == "normalized")
            {
                if(line.find("baseDensity") != string::npos)
                {
                    if(lineOld.find("RTstop") != string::npos)
                    {
                        //( $id, $len, $type, $rpkm, $scalingFactor, @baseDensities ) = split ( /\t/, $line );
                        //( $id, $len, $type, $rpkm, $scalingFactor, @rtstops ) = split ( /\t/, $lineOld );

                        id = data[0];
                        len = stoul(data[1]);
                        type = data[2];
                        rpkm = data[3];
                        scalingFactor = data[4];
                        for(auto it=data.cbegin()+5; it!=data.cend(); it++) baseDensities.push_back( stod(*it) );
                        // baseDensities.assign(data.cbegin()+5, data.cend());

                        StringArray oldData;
                        split(lineOld, '\t', oldData);
                        for(auto it=oldData.cbegin()+5; it!=oldData.cend(); it++) rtstops.push_back( stod(*it) );
                        //rtstops.assign(oldData.cbegin()+5, oldData.cend());
                    }
                }else if(line.find("RTstop") != string::npos){
                    if(lineOld.find("baseDensity") != string::npos){
                        
                        id = data[0];
                        len = stoul(data[1]);
                        type = data[2];
                        rpkm = data[3];
                        scalingFactor = data[4];
                        for(auto it=data.cbegin()+5; it!=data.cend(); it++) rtstops.push_back( stod(*it) );
                        //rtstops.assign(data.cbegin()+5, data.cend());

                        StringArray oldData;
                        split(lineOld, '\t', oldData);
                        for(auto it=oldData.cbegin()+5; it!=oldData.cend(); it++) baseDensities.push_back( stod(*it) );
                        //baseDensities.assign(oldData.cbegin()+5, oldData.cend());
                    }
                }
            }else{
                    id = data[0];
                    len = stoul(data[1]);
                    rpkm = data[2];
                    for(auto it=data.cbegin()+3; it!=data.cend(); it++) rtstops.push_back( stod(*it) );
                    //baseDensities.assign(data.cbegin()+4, data.cend());

                    StringArray oldData;
                    split(lineOld, '\t', oldData);
                    for(auto it=oldData.cbegin()+3; it!=oldData.cend(); it++) baseDensities.push_back( stod(*it) );
                    //rtstops.assign(oldData.cbegin()+4, oldData.cend());
            }

            if( trans_len.find(id) == trans_len.end() )
            //if(true)
            {
                trans_len[id] = len;
                trans_rpkm[id] = rpkm;
                if(format == "normalized")
                    trans_scalingFactor[id] = scalingFactor;
                trans_rtstop[id] = rtstops;
                trans_bd[id] = baseDensities;
            }else{
                if(len != trans_len[id])
                    throw Unexpected_Error("FATAL Error: "+id+" different transcript length "+to_string(len)+"/"+to_string(trans_len[id]));

                trans_rpkm[id] += ","+rpkm;
                if(format == "normalized")
                    trans_scalingFactor[id] += ","+scalingFactor;

                DoubleArray &ref_trans_rtstop = trans_rtstop[id];
                DoubleArray &ref_trans_bd = trans_bd[id];
                for(int idxPos=0; idxPos<=len; idxPos++)
                {
                    ref_trans_rtstop[idxPos] += rtstops[idxPos];
                    ref_trans_bd[idxPos] += baseDensities[idxPos];
                }
            }
        }else{
            idOld = data[0];
            lineOld = line;
        }
    }
    IN.close();
}

void outputSignals(const string &rt_file, const MapStringuLONG &trans_len, const MapStringString &trans_rpkm, 
                    const MapStringString &trans_scalingFactor, const MapStringT<DoubleArray> &trans_rtstop, 
                    const MapStringT<DoubleArray> &trans_bd, const string &format)
{
    cerr << "output signal to " << rt_file << "\n\t" << currentDateTime() << endl;
    ofstream OUT(rt_file, ofstream::out);
    if(not OUT)
    {
        cerr << "FATAL Error: cannot write to " << rt_file << "..." << endl;
        exit(-1);
    }
    //OUT << std::fixed << std::setprecision(3);
    string line; uLONG lineCount(0);

    cerr << "Will write " << trans_len.size() << " transcripts" << endl;
    for(auto item=trans_len.cbegin(); item!=trans_len.cend(); item++)
    {
        const string &id = item->first;
        
        OUT << id << "\t" << trans_len.at(id) << "\t" << trans_rpkm.at(id);
        if(format == "normalized") OUT << "\t" << trans_scalingFactor.at(id);
        for(const double &bd: trans_bd.at(id))
            OUT << "\t" << bd;
        OUT << "\n";

        OUT << id << "\t" << trans_len.at(id) << "\t" << trans_rpkm.at(id);
        if(format == "normalized") OUT << "\t" << trans_scalingFactor.at(id);
        for(const double &rt: trans_rtstop.at(id))
            OUT << "\t" << rt;
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

    MapStringuLONG trans_len;
    MapStringString trans_rpkm;
    MapStringString trans_scalingFactor;
    MapStringT<DoubleArray> trans_rtstop;
    MapStringT<DoubleArray> trans_bd;

    cout << "start to read..." << endl;
    for(auto file_name: param.input_rt_list)
        combineSignals(file_name, trans_len, trans_rpkm, trans_scalingFactor, trans_rtstop, trans_bd, param.format);

    cout << "start to write..." << endl;
    outputSignals(param.output_rt, trans_len, trans_rpkm, trans_scalingFactor, trans_rtstop, trans_bd, param.format);

    return 0;
}




