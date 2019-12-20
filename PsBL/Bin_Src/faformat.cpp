/*

    A Programe to format fasta file(sort, format output, sub-sample)

*/


//#include "paris_plot.h"
#include "param.h"
#include "fasta.h"
#include "version.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <functional>
#include <cstring>
#include <numeric>
#include <chrono>       // std::chrono::system_clock
#include <random>
#include <iomanip>
#include <regex>
#include <list>

using namespace std;
using namespace pan;

#define CALL_FAFORMAT_VERSION "1.000"
#define VERSION_DATE "2017-11-25"
#define COMPILE_DATE __DATE__

using std::accumulate;

Color::Modifier RED(Color::FG_RED);
Color::Modifier DEF(Color::FG_DEFAULT);

void print_usage()
{
    char buff[2000];
    const char *help_info = 
            "faformat - format fasta file(sort, format output, sub-sample)\n"
            "=============================================================\n"
            "\e[1mUSAGE:\e[0m\n"
            "\tfaformat -in input_fasta -out output(stdout) [ -sort len|chr_id -reverse -num 60 -all -sto -append\n"
            "\t         -remove chr_id_1,chr_id_2... -sample number -fetch chr_id_1,chr_id_2... ]\n"
            "\e[1mHELP:\e[0m\n"
            "\t[order]\n"
            "\t-sort: len|chr_id, output fasta sort by sequence length or chr_id(default: no sort)\n"
            "\t-reverse: sort all sequence in reverse order(default: no)\n\n"

            "\t[output format]\n"
            "\t-num: base number of each line(default: 60)\n"
            "\t-all: output all base in a line(default: no)\n"
            "\t-append: append the output to the file(default: no)\n"
            "\t-sto: output sto format(default: no)\n\n"

            "\t[sub-sample]\n"
            "\t-remove: remove some chr from input file(default: no remove)\n"
            "\t-sample: sample some chr from input file(default: no sample)\n"
            "\t-fetch: get some chr from input file(default: no fetch)\n\n"
            
            "\t[sub-sample pattern]\n"
            "\t-rp_chrid: <regex> remove some chr from input file whose chr_id match the regex \n"
            "\t-rp_anno: <regex> get some chr from input file whose chr annotation match the regex\n"
            "\t-fp_chrid: <regex> remove some chr from input file whose chr_id match the regex\n"
            "\t-fp_anno: <regex> get some chr from input file whose chr annotation match the regex\n\n"

            "\e[1mHELP:\e[0m\n"
            "\e[1mVERSION:\e[0m\n\t%s\n"
            "\e[1mLIB VERSION:\e[0m\n\t%s\n"
            "\e[1mVERSION DATE:\e[0m\n\t%s\n"
            "\e[1mCOMPILE DATE:\e[0m\n\t%s\n"
            "\e[1mAUTHOR:\e[0m\n\t%s\n";
    sprintf(buff, help_info, CALL_FAFORMAT_VERSION, VERSION, VERSION_DATE, COMPILE_DATE, "Li Pan");
    cout << buff << endl;
}

struct Param
{
    enum SORT_METHOD{ NO_SORT, LENGTH, CHR_ID };

    string input_fasta;
    string output_fasta;
    
    SORT_METHOD sort_method = NO_SORT;
    bool reverse = false;

    uLONG number_each_line = 60;
    bool write_sto = false;

    bool append = false;

    // sample
    StringArray remove_list;
    StringArray fetch_list;
    uINT sample_num = 0;

    // sub-sample pattern
    string rev_chrid_pattern;
    string rev_anno_pattern;
    string fet_chrid_pattern;
    string fet_anno_pattern;

    string param_string;

    operator bool()
    { 
        if(input_fasta.empty())
            return false;

        if( (not remove_list.empty()) + (not fetch_list.empty()) + (sample_num != 0) >= 2 )
        {
            //cerr << RED << "FATAL Error: -remove/-sample/-fetch are exclusive " << endl; 
            return false;
        }

        return true;
    }
};

void has_next(int argc, int current)
{
    if(current + 1 >= argc)
    {
        cerr << RED << "FATAL ERROR: Parameter Error" << DEF << DEF << endl;
        print_usage();
        exit(-1);
    }
}

Param read_param(int argc, char *argv[])
{
    Param param;

    bool find_remove(false), find_sample(false), find_fetch(false), find_num(false), find_all(false), find_sto(false);
    bool match_chrid(false), match_anno(false), remove_chrid(false), remove_anno(false);

    for(size_t i=1; i<argc; i++)
    {
        if( argv[i][0] == '-' )
        {
            if(not strcmp(argv[i]+1, "in"))
            {
                has_next(argc, i);
                param.input_fasta = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "out"))
            {
                has_next(argc, i);
                param.output_fasta = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "sort"))
            {
                has_next(argc, i);
                if( not strcmp(argv[i+1], "len") )
                    param.sort_method = Param::LENGTH;
                else if( not strcmp(argv[i+1], "chr_id") )
                    param.sort_method = Param::CHR_ID;
                else{
                    cerr << RED << "FATAL ERROR: unknown sort type: " << argv[i+1] << DEF << endl;
                    print_usage();
                    exit(-1);
                }
                i++;
            }else if(not strcmp(argv[i]+1, "reverse"))
            {
                //has_next(argc, i);
                param.reverse = true;
                //i++;
                //find_num = true;
            }else if(not strcmp(argv[i]+1, "num"))
            {
                has_next(argc, i);
                param.number_each_line = stoul(argv[i+1]);
                i++;
                find_num = true;
            }else if(not strcmp(argv[i]+1, "all"))
            {
                //has_next(argc, i);
                param.number_each_line = -1UL;
                //i++;
                find_all = true;
            }else if(not strcmp(argv[i]+1, "sto"))
            {
                //has_next(argc, i);
                param.write_sto = true;
                //i++;
                find_sto = true;
            }else if(not strcmp(argv[i]+1, "append"))
            {
                param.append = true;
            }else if(not strcmp(argv[i]+1, "remove"))
            {
                has_next(argc, i);
                const string remove_list = argv[i+1];
                param.remove_list = split(remove_list, ',');
                i++;
                find_remove = true;
            }else if(not strcmp(argv[i]+1, "fetch"))
            {
                has_next(argc, i);
                const string fetch_list = argv[i+1];
                param.fetch_list = split(fetch_list, ',');
                i++;
                find_fetch = true;
            }else if(not strcmp(argv[i]+1, "sample"))
            {
                has_next(argc, i);
                //const string remove_list = argv[i+1];
                param.sample_num = stoul(argv[i+1]);
                i++;
                find_sample = true;
            }else if(not strcmp(argv[i]+1, "rp_chrid"))
            {
                has_next(argc, i);
                param.rev_chrid_pattern = argv[i+1];
                i++;
                match_chrid = true;
            }else if(not strcmp(argv[i]+1, "rp_anno"))
            {
                has_next(argc, i);
                param.rev_anno_pattern = argv[i+1];
                i++;
                match_anno = true;
            }else if(not strcmp(argv[i]+1, "fp_chrid"))
            {
                has_next(argc, i);
                param.fet_chrid_pattern = argv[i+1];
                i++;
                remove_chrid = true;
            }else if(not strcmp(argv[i]+1, "fp_anno"))
            {
                has_next(argc, i);
                param.fet_anno_pattern = argv[i+1];
                i++;
                remove_anno = true;
            }else if(not strcmp(argv[i]+1, "h"))
            {
                print_usage();
                exit(0);
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

    if( find_num + find_all + find_sto >= 2 )
    {
        cerr << RED << "FATAL Error: -num/-all/-sto are exclusive " << endl; 
        exit(-1);
    }

    if( find_remove + find_sample + find_fetch + match_chrid + match_anno + remove_chrid + remove_anno >= 2 )
    {
        cerr << RED << "FATAL Error: -remove/-sample/-fetch/-rp_chrid/-rp_anno/-fp_chrid/-fp_anno are exclusive " << endl; 
        exit(-1);
    }

    param.param_string = param_string(argc, argv);
    return param;
}

void faformat(const Param &param)
{
    Fasta fasta(param.input_fasta);
    StringArray chr_ids = fasta.get_chr_ids();
    list<string> chr_ids_list(chr_ids.cbegin(), chr_ids.cend());

    //ofstream OUT(param.output_fasta, param.append ? ofstream::app : ofstream::out);
    ofstream OUT;

    if(!param.output_fasta.empty())
        OUT.open(param.output_fasta, param.append ? ofstream::app : ofstream::out);
    else
        OUT.basic_ios<char>::rdbuf(std::cout.rdbuf());

    if(not OUT)
    {
        cerr << "FATAL Error: " << param.output_fasta << " is unwritable" << endl;
        exit(-1);
    }

    // -remove
    if(not param.remove_list.empty())
        for(auto iter=param.remove_list.cbegin(); iter!=param.remove_list.cend(); iter++)
        {
            auto pos = find(chr_ids_list.cbegin(), chr_ids_list.cend(), *iter);
            chr_ids_list.erase( pos );
        }

    // -fetch
    if(not param.fetch_list.empty())
    {
        chr_ids_list.clear();
        for(auto iter=param.fetch_list.cbegin(); iter!=param.fetch_list.cend(); iter++)
            if( fasta.has_chr(*iter) )
                chr_ids_list.push_back(*iter);
    }

    // -sample
    if(param.sample_num != 0)
    {
        if(param.sample_num >= chr_ids_list.size())
        {
            cerr << "FATAL Error: -sample larger than fasta sequence number" << endl;
            exit(-1);
        }
        StringArray old_trans_list(chr_ids);
        chr_ids.clear();

        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::shuffle(old_trans_list.begin(), old_trans_list.end(), std::default_random_engine(seed));
        chr_ids_list.assign(old_trans_list.cbegin(), old_trans_list.cbegin()+param.sample_num);
    }

    //if (std::regex_search ("KU501215.1", std::regex("^KU") ))
    //    std::cout << "string literal matched\n";

    // -rp_chrid
    if(not param.rev_chrid_pattern.empty())
    {
        //cout << "rev_chrid_pattern: " << param.rev_chrid_pattern << endl;
        regex reg(param.rev_chrid_pattern);
        for(auto iter=chr_ids_list.begin(); iter!=chr_ids_list.end(); )
        {
            //cout << "regex_match(" << *iter << ", reg)" << endl; 
            if( regex_search(*iter, reg) )
            {
                //cout << "rem: " << *iter << endl;
                iter = chr_ids_list.erase( iter );
            }else{
                iter++;
            }
        }
    }

    // -rp_anno
    if(not param.rev_anno_pattern.empty())
    {
        //cout << "rev_anno_pattern" << endl;
        regex reg(param.rev_anno_pattern);
        for(auto iter=chr_ids_list.begin(); iter!=chr_ids_list.end(); )
        {
            if( regex_search(fasta.get_chr_anno(*iter), reg) )
            {
                iter = chr_ids_list.erase( iter );
            }else{
                iter++;
            }
        }
    }

    // -fp_chrid
    if(not param.fet_chrid_pattern.empty())
    {
        //cout << "fet_chrid_pattern" << endl;
        regex reg(param.fet_chrid_pattern);
        for(auto iter=chr_ids_list.begin(); iter!=chr_ids_list.end(); )
        {
            if( not regex_search(*iter, reg) )
            {
                iter = chr_ids_list.erase( iter );
            }else{
                iter++;
            }
        }
    }

    // -fp_anno
    if(not param.fet_anno_pattern.empty())
    {
        //cout << "fet_anno_pattern" << endl;
        regex reg(param.fet_anno_pattern);
        for(auto iter=chr_ids_list.begin(); iter!=chr_ids_list.end(); )
        {
            if( not regex_search(fasta.get_chr_anno(*iter), reg) )
            {
                iter = chr_ids_list.erase( iter );
            }else{
                iter++;
            }
        }
    }

    //cout << chr_ids << endl;
    chr_ids.assign(chr_ids_list.cbegin(), chr_ids_list.cend());

    // output
    if(param.sort_method == Param::CHR_ID)
    {
        if(not param.reverse)
            sort(chr_ids.begin(), chr_ids.end());
        else
            sort(chr_ids.rbegin(), chr_ids.rend());
    }else if(param.sort_method == Param::LENGTH)
    {
        vector<pair<string, uLONG>> chr_lens;
        for(const string &chr_id: chr_ids)
            chr_lens.push_back( make_pair( chr_id, fasta.get_chr_len(chr_id) ) );

        if(not param.reverse)
            sort(chr_lens.begin(), chr_lens.end(), [](const pair<string, uLONG> &p1, const pair<string, uLONG> &p2){ return p1.second < p2.second; });
        else
            sort(chr_lens.begin(), chr_lens.end(), [](const pair<string, uLONG> &p1, const pair<string, uLONG> &p2){ return p1.second > p2.second; });
        
        chr_ids.clear();
        for_each(chr_lens.cbegin(), chr_lens.cend(), [&chr_ids](const pair<string, uLONG> &p){ chr_ids.push_back(p.first); });
    }else{
        //cerr << "Unknown sort method..." << endl;
        // no sort
    }

    uLONG longest_name(0);
    if(param.write_sto)
        for(auto iter=chr_ids.cbegin(); iter!=chr_ids.cend(); iter++)
            longest_name = max(longest_name, static_cast<uLONG>(iter->size()));

    if(param.write_sto)
    {
        OUT << "# STOCKHOLM 1.0\n\n";
        bool has_out(false);
        for(auto iter=chr_ids.cbegin(); iter!=chr_ids.cend(); iter++)
        {
            const string &anno = fasta.get_chr_anno(*iter);
            if(not anno.empty())
            {
                has_out = true;
                OUT << "#=GS " << setw(longest_name+5) << std::left << *iter << setw(0) << "DE " << anno << "\n";
            }
        }
        if(has_out)
            OUT << "\n";
    }

    for(auto iter=chr_ids.cbegin(); iter!=chr_ids.cend(); iter++)
    {
        if(param.write_sto)
        {
            OUT << setw(longest_name+5) << std::left << *iter << setw(0) << fasta.get_chr_seq(*iter) << "\n";
        }else{
            OUT << ">" << *iter;
            const string &anno = fasta.get_chr_anno(*iter);
            if(not anno.empty())
                OUT << "\t" << anno;
            OUT << "\n" << flat_seq( fasta.get_chr_seq(*iter), param.number_each_line );
        }
    }

    if(param.write_sto)
        OUT << "\n//\n";
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

    faformat(param);

    return 0;
}


















