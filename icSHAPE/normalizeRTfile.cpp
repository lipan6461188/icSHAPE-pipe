
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

Color::Modifier RED(Color::FG_RED);
Color::Modifier DEF(Color::FG_DEFAULT);
Color::Modifier YELLOW(Color::FG_YELLOW);

void print_usage()
{
    char buff[8000];
    const char *help_info = 
            "## --------------------------------------\n"
            "Calculate average signals\n\n"

            "Command:\n"
            "%s -i input_baseDensity_file -o output_normalized_baseDensity_file \n"
            " # what it is:\n"
            " -i     input_baseDensity_file\n"
            " -o     output_normalized_baseDensity_file\n\n"

            "# more options:\n"
            " -d     head to skip (default: 32)\n"
            " -l     tail to skip (default: 32)\n"
            " -f     scaling form (the benchmarked value will be scaled to this) (default: 100)\n\n"

            " -m     normalize method (default: mean:vigintile2)\n\n"
            
            " -v     verbose mode \n"
            " -r     raw mode (100%% identical with perl icSHAPE pipeline, because I correct a bug from it)\n\n";

    sprintf(buff, help_info, "normalizeRTfile");
    cout << buff << endl;
}

enum NormSampleMethod { SMART, UPPER, QUARTILE, DECILE, VIGINTILE };
enum NormCalcMethiod { MEDIAN, MEAN, PEAK };

struct Param
{
    string input_file;
    string output_file;

    uINT head_skip = 32;
    uINT tail_skip = 32;
    
    double scalling_form = 100.0;

    string norm_method = "mean:vigintile2";

    NormSampleMethod norm_sample_method = DECILE;
    NormCalcMethiod norm_calc_method = MEDIAN;
    uINT norm_sample_factor = 2;
    bool log_op = false;

    bool verbose = false;
    bool raw_mode = false;

    operator bool() const { 
        if( input_file.empty() or output_file.empty() ) return false;
        return true;
    }

    string show_norm_method()
    {
        string norm_method;

        if(norm_calc_method==MEDIAN)
            norm_method += "median";
        else if(norm_calc_method==MEAN)
            norm_method += "mean";
        else
            norm_method += "peak";

        norm_method += ":";

        if(norm_sample_method==SMART)
            norm_method += "smart";
        else if(norm_sample_method==UPPER)
            norm_method += "upper";
        else if(norm_sample_method==QUARTILE)
            norm_method += "quartile";
        else if(norm_sample_method==DECILE)
            norm_method += "decile";
        else if(norm_sample_method==VIGINTILE)
            norm_method += "vigintile";
        
        norm_method += to_string(norm_sample_factor);

        if(log_op)
        {
            norm_method += ":log";
        }

        return norm_method;
    }
};


void parse_norm_method(Param &param)
{
    StringArray data;
    split(param.norm_method, ':', data);
    for(const string &item: data)
    {
        if(item == "smart")
            param.norm_sample_method = SMART;
        else if(item == "upper")
            param.norm_sample_method = UPPER;
        else if(item.substr(0,8) == "quartile")
        {
            param.norm_sample_method = QUARTILE;
            param.norm_sample_factor = stoul(item.substr(8));
        }else if(item.substr(0,6) == "decile")
        {
            param.norm_sample_method = DECILE;
            param.norm_sample_factor = stoul(item.substr(6));
        }
        else if(item.substr(0,9) == "vigintile")
        {
            param.norm_sample_method = VIGINTILE;
            param.norm_sample_factor = stoul(item.substr(9));
        }
        else if(item == "median")
            param.norm_calc_method = MEDIAN;
        else if(item == "mean")
            param.norm_calc_method = MEAN;
        else if(item == "peak")
            param.norm_calc_method = PEAK;
        else if(item == "log")
            param.log_op = true;
        else{
            cerr << RED << "Error: invalid -m parameter" << DEF << endl;
            exit(-1);
        }
    }
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
            }else if(not strcmp(argv[i]+1, "d"))
            {
                has_next(argc, i);
                param.head_skip = stoul(argv[i+1]);
                i++;
            }else if(not strcmp(argv[i]+1, "l"))
            {
                has_next(argc, i);
                param.tail_skip = stoul(argv[i+1]);
                i++;
            }else if(not strcmp(argv[i]+1, "f"))
            {
                has_next(argc, i);
                param.scalling_form = stod(argv[i+1]);
                i++;
            }else if(not strcmp(argv[i]+1, "m"))
            {
                has_next(argc, i);
                param.norm_method = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "v"))
            {
                //has_next(argc, i);
                param.verbose = true;
                //i++;
            }else if(not strcmp(argv[i]+1, "r"))
            {
                //has_next(argc, i);
                param.raw_mode = true;
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

    parse_norm_method(param);
    return param;
}

double calcScalingFactor(const DoubleArray &data_array, const uINT &start, 
                        const uINT &end, const Param &param)
{
    DoubleArray sorted_data(data_array.cbegin()+start, data_array.cbegin()+end+1);
    sort(sorted_data.rbegin(), sorted_data.rend());

    //cout << start << "\t" << end << endl;
    //for(double v: sorted_data)
    //    cout << " " << v;
    //cout << " " << sorted_data.size() << "\n";

    uINT len = sorted_data.size();
    DoubleArray rangeOfSelection;
    if(param.norm_sample_method == SMART)
    {
        for(int idx=0; idx<len; idx++)
        {
            if(sorted_data[idx] <= 0) break;
            rangeOfSelection.push_back(sorted_data[idx]);
        }
    }else{
        uINT selectStart = 0;
        uINT selectEnd = len - 1;

        double f = param.norm_sample_factor;

        if(param.norm_sample_method == UPPER)
            selectEnd = len / 2 - 1;
        else if(param.norm_sample_method == QUARTILE)
        {
            selectStart = (f-1)*len/4;
            selectEnd = int(f*len/4) - 1;
        }else if(param.norm_sample_method == DECILE)
        {
            selectStart = (f-1)*len / 10;
            selectEnd = int(f*len/10) - 1;
        }else if(param.norm_sample_method == VIGINTILE)
        {
            selectStart = (f-1)*len / 20;
            selectEnd = int(f*len/20) - 1;
        }else{
            cerr << RED << "FATAL Error: invalid norm_method" << DEF << endl;
            exit(-1);
        }

        for(auto it=sorted_data.cbegin()+selectStart; it!=sorted_data.cbegin()+selectEnd+1; it++)
            rangeOfSelection.push_back(*it);
    }

    len = rangeOfSelection.size();

    double scalling_factor = 1.0;
    if(param.norm_calc_method == MEDIAN)
    {
        double median = 1;
        if ( len % 2 == 0 ) {  median = (rangeOfSelection[len/2-1] + rangeOfSelection[len/2]) /2;  }
        else {  median = rangeOfSelection[(len-1)/2];  }
        scalling_factor = median;

    }else if(param.norm_calc_method == MEAN)
    {
        double Sum;
        for(const double &v: rangeOfSelection) Sum += v;
        scalling_factor = Sum/rangeOfSelection.size();

    }else if(param.norm_calc_method == PEAK)
    {
        scalling_factor = rangeOfSelection[0];
    }

    return scalling_factor;
}
/*
void read_full_file(const Param &param, vector<StringArray> &contents)
{
    contents.clear();

    const string input_file(param.input_file);

    cerr << "Read " << input_file << " now " << currentDateTime() << endl;

    ifstream IN(input_file, ifstream::in);

    if(not IN)
    {
        cerr << RED << "FATAL Error: cannot open " << input_file << DEF << endl;
        exit(-1);
    }

    string line;
    while(getline(IN, line))
    {
        if(line[0] == '#') continue;
        StringArray data;
        trim(line);
        split(line, '\t', data);
        contents.push_back(data);
    }
}

void process_full_file(const Param &param, StringMatrix &contents)
{
    cerr << "Normalize base density from file $baseDensityFile...\n\t" << currentDateTime() << endl;

    const string output_file(param.output_file);
    const string norm_method(param.norm_method);
    uINT headToSkip(param.head_skip);
    uINT tailToSkip(param.tail_skip);
    const double scalling_form(param.scalling_form);

    headToSkip++;
    ofstream OUT(output_file, ofstream::out);

        if(not OUT)
    {
        cerr << RED << "FATAL Error: cannot write to " << output_file << DEF << endl;
        exit(-1);
    }
    OUT << "#transcript\tlength\ttype\tbase_frequency, start from position 1.\n";

    string transcript;
    uLONG len = 0;
    string rpkm;
    
    uLONG lineCount = 0;
    for(uLONG i=0; i<contents.size(); i++)
    {
        ++lineCount;
        if(lineCount % 1000 == 0) 
            cerr << "  line: " << lineCount << endl;

        DoubleArray baseDensity, rtstop;
        double scalling_factor = 1.0;
        //trim(line);

        const StringArray &data = contents[i];
        //StringArray data;
        //split(line, '\t', data);
        transcript = data[0];
        len = stoul(data[1]);
        rpkm = data[2];
        for(auto it=data.cbegin()+3; it!=data.cend(); it++) baseDensity.push_back( stod(*it) );

        uINT trimed_last = len - tailToSkip;
        while(trimed_last < headToSkip+40)
        {
            headToSkip /= 2;
            tailToSkip /= 2;
            trimed_last = len - tailToSkip;
            cerr << "Warning! Transcript $transcript too short. update headToSkip: " << headToSkip << "; tailToSkip: " << tailToSkip << endl;
            if(headToSkip==0 and tailToSkip==0) break;
        }

        scalling_factor = calcScalingFactor(baseDensity, headToSkip, trimed_last, param);

        if(scalling_factor > 1)
        {
            scalling_factor = scalling_factor / scalling_form;
            OUT << transcript << "\t" << len << "\tbaseDensity\t" << rpkm << "\t" << scalling_factor;
            if(norm_method == "log")
            {
                for(int i=1; i<=len; i++)
                {
                    char buff[400];
                    sprintf(buff, "%.3f", log2(baseDensity[i]/scalling_factor+1));
                    OUT << "\t" << buff;
                }
            }else{
                for(int i=1; i<=len; i++)
                {
                    char buff[400];
                    sprintf(buff, "%.3f", baseDensity[i]/scalling_factor);
                    OUT << "\t" << buff;
                }
            }
            OUT << "\n";
        }

        i++;
        const StringArray &rt_data = contents[i];

        //getline(IN, line);
        //trim(line);
        //split(line, '\t', data);
        for(auto it=rt_data.cbegin()+3; it!=rt_data.cend(); it++) rtstop.push_back( stod(*it) );
        scalling_factor = calcScalingFactor(rtstop, headToSkip, trimed_last, param);

        if(scalling_factor > 1)
        {
            scalling_factor = scalling_factor / scalling_form;

            OUT << transcript << "\t" << len << "\tRTstop\t" << rpkm << "\t" << scalling_factor;
            if(norm_method == "log")
            {
                for(int i=1; i<=len; i++)
                {
                    char buff[400];
                    sprintf(buff, "%.3f", log2(rtstop[i]/scalling_factor+1));
                    OUT << "\t" << buff;
                }
            }else{
                for(int i=1; i<=len; i++)
                {
                    char buff[400];
                    sprintf(buff, "%.3f", rtstop[i]/scalling_factor);
                    OUT << "\t" << buff;
                }
            }
            OUT << "\n";
        }
    }

    //IN.close();
    OUT.close();

}
*/

void normalizeBaseDensity(const Param &param)
{
    cerr << "Normalize base density from file $baseDensityFile...\n\t" << currentDateTime() << endl;

    const string input_file(param.input_file);
    const string output_file(param.output_file);
    const bool log_op(param.log_op);
    uINT headToSkip(param.head_skip+1);
    uINT tailToSkip(param.tail_skip);
    const double scalling_form(param.scalling_form);

    ifstream IN(input_file, ifstream::in);
    ofstream OUT(output_file, ofstream::out);

    if(not IN)
    {
        cerr << RED << "FATAL Error: cannot open " << input_file << DEF << endl;
        exit(-1);
    }
    if(not OUT)
    {
        cerr << RED << "FATAL Error: cannot write to " << output_file << DEF << endl;
        exit(-1);
    }
    OUT << "#transcript\tlength\ttype\tbase_frequency, start from position 1.\n";

    string transcript;
    uLONG len = 0;
    string rpkm;
    
    string line;
    uLONG lineCount = 0;
    while(getline(IN, line))
    {
        if(line[0] == '#') continue;

        if(not param.raw_mode)
        {
            headToSkip = param.head_skip+1;
            tailToSkip = param.tail_skip;
        }

        ++lineCount;
        if(lineCount % 1000 == 0) 
            cerr << "  line: " << lineCount << endl;

        DoubleArray baseDensity, rtstop;
        double scalling_factor = 1.0;
        trim(line);

        StringArray data;
        split(line, '\t', data);
        transcript = data[0];
        len = stoul(data[1]);
        rpkm = data[2];
        for(auto it=data.cbegin()+3; it!=data.cend(); it++) baseDensity.push_back( stod(*it) );

        uINT trimed_last = len - tailToSkip;
        while(trimed_last < headToSkip+40)
        {
            headToSkip /= 2;
            tailToSkip /= 2;
            trimed_last = len - tailToSkip;
            if(not param.verbose)
                cerr << "Warning! Transcript $transcript too short. update headToSkip: " << headToSkip << "; tailToSkip: " << tailToSkip << endl;
            if(headToSkip==0 and tailToSkip==0) break;
        }

        scalling_factor = calcScalingFactor(baseDensity, headToSkip, trimed_last, param);

        if(scalling_factor > 1)
        {
            scalling_factor = scalling_factor / scalling_form;
            OUT << transcript << "\t" << len << "\tbaseDensity\t" << rpkm << "\t" << scalling_factor;
            
            auto ss = OUT.precision();
            OUT << std::fixed << std::setprecision(3);
            if(log_op)
            {
                for(int i=1; i<=len; i++)
                {
                    char buff[400];
                    sprintf(buff, "%.3f", log2(baseDensity[i]/scalling_factor+1));
                    OUT << "\t" << buff;
                }
            }else{
                for(int i=1; i<=len; i++)
                {
                    //char buff[400];
                    //cout << baseDensity[i] << "/" << scalling_factor << endl;
                    //sprintf(buff, "%.3f", baseDensity[i]/scalling_factor);
                    OUT << "\t" << baseDensity[i]/scalling_factor;
                }
                
            }
            OUT.precision(ss);
            OUT << "\n";
        }

        getline(IN, line);
        trim(line);
        split(line, '\t', data);
        for(auto it=data.cbegin()+3; it!=data.cend(); it++) rtstop.push_back( stod(*it) );
        scalling_factor = calcScalingFactor(rtstop, headToSkip, trimed_last, param);

        if(scalling_factor > 1)
        {
            scalling_factor = scalling_factor / scalling_form;

            OUT << transcript << "\t" << len << "\tRTstop\t" << rpkm << "\t" << scalling_factor;
            if(log_op)
            {
                for(int i=1; i<=len; i++)
                {
                    char buff[400];
                    sprintf(buff, "%.3f", log2(rtstop[i]/scalling_factor+1));
                    OUT << "\t" << buff;
                }
            }else{
                for(int i=1; i<=len; i++)
                {
                    char buff[400];
                    sprintf(buff, "%.3f", rtstop[i]/scalling_factor);
                    OUT << "\t" << buff;
                }
            }
            OUT << "\n";
        }else{
            if(param.verbose)
                cerr << "Filter RTstop of " << transcript << " for small scalling_factor=" << scalling_factor << endl;
        }
    }

    IN.close();
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

    cout << "normalize method: " << param.show_norm_method() << endl;
    normalizeBaseDensity(param);

    //StringMatrix rt_data;
    //read_full_file(param, rt_data);
    //process_full_file(param, rt_data);

    return 0;
}


