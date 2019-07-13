/*

    A Programe to cluster PARIS dg

*/

#include "paris.h"
#include "param.h"
#include "fold.h"
#include "fasta.h"
#include "version.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <sstream>

using namespace std;
using namespace pan;

#define SAM2DG_VERSION "1.000"
#define DATE __DATE__

Color::Modifier RED(Color::FG_RED);
Color::Modifier DEF(Color::FG_DEFAULT);
Color::Modifier YELLOW(Color::FG_YELLOW);

void print_usage()
{
    char buff[2000];
    const char *help_info = 
            "dg_cluster - cluster PARIS dg\n"
            "=========================\n"
            "\e[1mUSAGE:\e[0m\n"
            "\tsam2dg -in input_dg -out output_dg [-out_tab file_name -min_overlap 5 -multiDG no   \n"
            "\t                    -max_gap 10 -max_total 30 -check_reads no                       \n"
            "\t                    -tag_sam input_sam_1,input_sam_2...,output_sam                  \n"
            "\t                    -min_support 2 -uniqDG no -fasta input_fasta ]                  \n"
            "\e[1mHELP:\e[0m\n"
            "\t-input_dh: input unclustered duplex group file from sam2dg\n"
            "\t-output_dg: output clustered duplex group\n"
            "\t-fasta: genome fasta file\n\n"

            "\t[Step 1 -- Culster duplex groups]\n"
            "\t-min_overlap: minimum overlap between 2 reads to cluster (default: 5)\n"
            "\t-multiDG: yes/no. allow a read cluster to multiple duplex group, which\n"
            "\t          will make tagged read in multiple duplex group ambigous (default: no)\n"
            "\t-uniqDG: yes/no. remove duplicate duplex group (defult: no)\n\n"

            "\t[Step 2 -- Collapse duplex groups]\n"
            "\t-max_gap: max gap between 2 clustered duplex group (default: 10)\n"
            "\t-max_total: max length of left/right arm of duplex group (default: 30)\n"
            "\t-check_reads: yes/no. if check reads when collappsed 2 group (default: no)\n\n"

            "\t[Step 3 -- tag sam reads]\n"
            "\t-tag_sam: input,output sam of sam2dg, multiple input sam files is allowed \n\n"
            
            "\t[Step 4 -- output]\n"
            "\t-min_support: mininum support reads for each duplex group (default: 2)\n\n"

            "\e[1mVERSION:\e[0m\n\t%s\n"
            "\e[1mLIB VERSION:\e[0m\n\t%s\n"
            "\e[1mDATE:\e[0m\n\t%s\n"
            "\e[1mAUTHOR:\e[0m\n\t%s\n";

    //ostringstream warning;
    //warning << YELLOW << WARNING << DEF;

    sprintf(buff, help_info, SAM2DG_VERSION, VERSION, DATE, "Li Pan");
    cout << buff << endl;
}

struct Param
{
    enum SORT_TYPE{ no_sort, balance, left_priority };

    string input_dg;
    string genome_fasta;
    string output_dg;
    string output_tab;

    uLONG min_overlap = 5;
    bool multiDG = false;
    bool uniqDG = false;

    uLONG max_gap = 10;
    uLONG max_total = 30;
    bool check_reads = false;

    StringArray input_tag_sam;
    string output_tag_sam;

    uLONG min_support = 2;

    operator bool(){ return input_dg.empty() or output_dg.empty() ? false : true; }
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
            if(not strcmp(argv[i]+1, "min_overlap"))
            {
                has_next(argc, i);
                param.min_overlap = stoul(string(argv[i+1]));
                i++;
            }else if(not strcmp(argv[i]+1, "multiDG"))
            {
                has_next(argc, i);

                if(not strcmp(argv[i+1], "yes"))
                    param.multiDG = true;
                else if(not strcmp(argv[i+1], "no"))
                    param.multiDG = false;
                else{
                    cerr << RED << "FATAL ERROR: unknown -multiDG option: " << argv[i+1] << DEF << endl;
                    exit(-1);   
                }
                i++;
            }else if(not strcmp(argv[i]+1, "uniqDG"))
            {
                has_next(argc, i);

                if(not strcmp(argv[i+1], "yes"))
                    param.uniqDG = true;
                else if(not strcmp(argv[i+1], "no"))
                    param.uniqDG = false;
                else{
                    cerr << RED << "FATAL ERROR: unknown -uniqDG option: " << argv[i+1] << DEF << endl;
                    exit(-1);   
                }
                i++;
            }else if(not strcmp(argv[i]+1, "in"))
            {
                has_next(argc, i);
                param.input_dg = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "fasta"))
            {
                has_next(argc, i);
                param.genome_fasta = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "out"))
            {
                has_next(argc, i);
                param.output_dg = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "out_tab"))
            {
                has_next(argc, i);
                param.output_tab = argv[i+1];
                i++;
            }else if(not strcmp(argv[i]+1, "h"))
            {
                print_usage();
                exit(-1);
            }else if(not strcmp(argv[i]+1, "check_reads"))
            {
                has_next(argc, i);

                if(not strcmp(argv[i+1], "yes"))
                    param.check_reads = true;
                else if(not strcmp(argv[i+1], "no"))
                    param.check_reads = false;
                else{
                    cerr << RED << "FATAL ERROR: unknown -check_reads option: " << argv[i+1] << DEF << endl;
                    exit(-1);   
                }
                i++;
            }else if(not strcmp(argv[i]+1, "max_gap"))
            {
                has_next(argc, i);
                param.max_gap = stoul(string(argv[i+1]));
                i++;
            }else if(not strcmp(argv[i]+1, "max_total"))
            {
                has_next(argc, i);
                param.max_total = stoul(string(argv[i+1]));
                i++;
            }else if(not strcmp(argv[i]+1, "tag_sam"))
            {
                has_next(argc, i);
                auto items = split( argv[i+1] , ',');
                if(items.size() < 2)
                {
                    cerr << "FATAL ERROR: bad -tag_sam option: " << argv[i+1] << DEF << endl;
                    exit(-1);
                }
                param.output_tag_sam = items.back();
                items.pop_back();
                param.input_tag_sam = items;
                i++;
            }else if(not strcmp(argv[i]+1, "min_support"))
            {
                has_next(argc, i);
                param.min_support = stoul(string(argv[i+1]));
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

struct Duplex_Group;

struct Ext_Duplex_Hang: public Duplex_Hang
{
    vector<Duplex_Group *> dgs;
    vector<uLONG> dg_ids;

    void finalize();
};

struct Duplex_Group
{
    uLONG start_1;
    uLONG end_1;
    uLONG start_2;
    uLONG end_2;
    vector< sp<Ext_Duplex_Hang> > reads;

    uLONG left_cov = 0;
    uLONG right_cov = 0;
    uLONG dg_id = -1UL;
    double score = 0.0;

    uLONG support() const { return reads.size(); }

    const string &chr_id_1()const { return reads.front()->chr_id_1; };
    const string &chr_id_2()const { return reads.front()->chr_id_2; };
    char strand_1()const { return reads.front()->strand_1; };
    char strand_2()const { return reads.front()->strand_2; };

    Duplex_Group( sp<Ext_Duplex_Hang> sp_dh );
    void add_duplex_hang( sp<Ext_Duplex_Hang> dh );
    long check_overlap(const Ext_Duplex_Hang &query_dh) const;
    long check_overlap(const Duplex_Group &query_dg, const uINT max_gap, const uINT max_total, const bool check_reads) const;
    void merge_duplex_group(const Duplex_Group &dg);

    Duplex_Group(const Duplex_Group& dg) = delete;
    Duplex_Group &operator=(Duplex_Group &) = delete;

    ~Duplex_Group()
    {
        for_each(reads.cbegin(), reads.cend(), [&](sp<Ext_Duplex_Hang> p_dh) -> void
        {
            p_dh->dgs.erase( find(p_dh->dgs.cbegin(), p_dh->dgs.cend(), this) );
        });
    }

    friend bool operator<(const Duplex_Group& dg_1, const Duplex_Group& dg_2);
};

ostream &operator<<(ostream &OUT, const Duplex_Group& dp)
{
    OUT << dp.chr_id_1() << "\t" << dp.strand_1() << "\t" << dp.start_1 << "\t" << dp.end_1 << "\t" 
        << dp.chr_id_2() << "\t" << dp.strand_2() << "\t" << dp.start_2 << "\t" << dp.end_2 << "\t"
        << "support:" << dp.support() << "\n"; 
    return OUT;
}

void Ext_Duplex_Hang::finalize()
{
    dg_ids.clear();
    for(Duplex_Group * p_dg: dgs)
    {
        if(p_dg->dg_id == -1UL)
        {
            cerr << "Warning: A undefined duplex group\n";
            continue;
        }
        dg_ids.push_back(p_dg->dg_id);
    }
}

inline Duplex_Group::Duplex_Group( sp<Ext_Duplex_Hang> sp_dh ): 
            start_1(sp_dh->start_1), end_1(sp_dh->end_1), 
            start_2(sp_dh->start_2), end_2(sp_dh->end_2)
            {  reads.push_back(sp_dh); sp_dh->dgs.push_back(this); }

void Duplex_Group::add_duplex_hang( sp<Ext_Duplex_Hang> sp_dh )
{
    //support++;
    start_1 = max(start_1, sp_dh->start_1);
    end_1 = min(end_1, sp_dh->end_1);
    start_2 = max(start_2, sp_dh->start_2);
    end_2 = min(end_2, sp_dh->end_2);
    reads.push_back(sp_dh);
    sp_dh->dgs.push_back(this);
}

/* return 0, -1, overlapped bases */
long Duplex_Group::check_overlap(const Ext_Duplex_Hang &query_dh) const
{
    long overlap = 0;
    if(end_1 < query_dh.start_1 or chr_id_1() != query_dh.chr_id_1
        or chr_id_2() != query_dh.chr_id_2 or strand_1() != query_dh.strand_1
        or strand_2() != query_dh.strand_2 )
        return -1;
    
    if(start_2 <= query_dh.end_2 and query_dh.start_2 <= end_2)
    {
        uINT left_overlap = min(end_1, query_dh.end_1) - query_dh.start_1 + 1;
        uINT right_overlap = min(end_2, query_dh.end_2) - max(start_2, query_dh.start_2) + 1;
        overlap = min(left_overlap, right_overlap);
    }

    return overlap;
}

bool quick_read_overlap(const Ext_Duplex_Hang &dh_1, const Ext_Duplex_Hang &dh_2)
{
    if(dh_1.start_1 < dh_2.end_1 and dh_2.start_1 < dh_1.end_1 and dh_1.start_2 < dh_2.end_2 and dh_2.start_2 < dh_2.end_2)
        return true;
    else
        return false;
}

long Duplex_Group::check_overlap(const Duplex_Group &query_dg, 
                            const uINT max_gap, 
                            const uINT max_total, 
                            const bool check_reads) const
{

    //cout << max_gap << "\t" << max_total << endl;
    //clog << "Show Self: \n\t" << *this;

    long overlap = 0;
    if( end_1 + max_gap < query_dg.start_1 or 
        chr_id_2() != query_dg.chr_id_2() or strand_2() != query_dg.strand_2() or
        chr_id_1() != query_dg.chr_id_1() or strand_1() != query_dg.strand_1() )
    {
        //clog << end_1 + max_gap << "<" << query_dg.start_1 << endl;
        return -1;
    }

    uINT gap_2(0);
    if(query_dg.start_2 > end_2)
        gap_2 = query_dg.start_2 - end_2 - 1;
    else if(query_dg.end_2 < start_2)
        gap_2 = start_2 - query_dg.end_2 - 1;

    uINT total_1 = max(end_1, query_dg.end_1) - start_1 + 1;
    uINT total_2 = max(end_2, query_dg.end_2) - min(start_2, query_dg.start_2) + 1;

    if(total_1 <= max_total and total_2 <= max_total and gap_2 <= max_gap)
        overlap = 1;
    
    if(overlap == 1 and check_reads)
    {
        //const uINT min_read_overlap(10);
        uINT overlap_count = 0;
        for( const sp<Ext_Duplex_Hang> dh_1: reads)
        {
            for(const sp<Ext_Duplex_Hang> dh_2: query_dg.reads)
            {
                if( quick_read_overlap(*dh_1, *dh_2) )
                {
                    if(++overlap_count >= 2)
                        return 1;
                }
            }
        }
        return 0;
    }

    return overlap;
}


void Duplex_Group::merge_duplex_group(const Duplex_Group &dg)
{
    for(sp<Ext_Duplex_Hang> sp_dh: dg.reads)
        //if(find(reads.cbegin(), reads.cend(), sp_dh) == reads.cend())
    {
        sp_dh->dgs.push_back(this);
        reads.push_back(sp_dh);
    }

    start_1 = min(start_1, dg.start_1);
    end_1 = max(end_1, dg.end_1);
    start_2 = min(start_2, dg.end_2);
    end_2 = max(end_2, dg.end_2);

    sort(reads.begin(), reads.end());
    auto back = unique(reads.begin(), reads.end());
    reads.erase(back, reads.end());
}

bool operator<(const Duplex_Group& dg_1, const Duplex_Group& dg_2)
{
    // chr 1
    if(dg_1.chr_id_1() < dg_2.chr_id_1())
        return true;
    else if(dg_1.chr_id_1() > dg_2.chr_id_1())
        return false;
    else{
        // strand 1
        if(dg_1.strand_1() < dg_2.strand_1())
            return true;
        else if(dg_1.strand_1() > dg_2.strand_1())
            return false;
        else{
            // chr 2
            if(dg_1.chr_id_2() < dg_2.chr_id_2())
                return true;
            else if(dg_1.chr_id_2() > dg_2.chr_id_2())
                return false;
            else{
                // strand 2
                if(dg_1.strand_2() < dg_2.strand_2())
                    return true;
                else if(dg_1.strand_2() > dg_2.strand_2())
                    return false;
                else{
                    // start 1
                    if(dg_1.start_1 < dg_2.start_1)
                        return true;
                    else if(dg_1.start_1 > dg_2.start_1)
                        return false;
                    else{
                        // end 1
                        if(dg_1.end_1 < dg_2.end_1)
                            return true;
                        else if(dg_1.end_1 > dg_2.end_1)
                            return false;
                        else
                        {
                            // start 2
                            if(dg_1.start_2 < dg_2.start_2)
                                return true;
                            else if(dg_1.start_2 > dg_2.start_2)
                                return false;
                            else
                            {
                                // end 2
                                if(dg_1.end_2 < dg_2.end_2)
                                    return true;
                                else if(dg_1.end_2 > dg_2.end_2)
                                    return false;
                                else{
                                    return false;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


class Genome_Covarege
{
public:
    Genome_Covarege( const vector< sp<Ext_Duplex_Hang> > &dh_list);
    uLONG max_cov(const string &chr_id, uLONG start, uLONG end) const;

private:
    //void build_index(const uINT interval);

    //using size_type=vector<Region>::size_type;

    MapStringT< vector<uLONG> > chr_regions;
    //MapStringT< vector<size_type> > start_index;
    //MapStringT< vector<size_type> > end_index;

    //uINT interval;
};

Genome_Covarege::Genome_Covarege(const vector< sp<Ext_Duplex_Hang> > &dh_list)
{
    // find the max length of each chr
    struct INIT_ULONG { uLONG Value = 0; };

    MapStringT< INIT_ULONG > max_chr_length;
    for(const sp<Ext_Duplex_Hang> sp_dh: dh_list)
    {
        max_chr_length[sp_dh->chr_id_1].Value = max( max_chr_length[sp_dh->chr_id_1].Value, sp_dh->end_1 );
        max_chr_length[sp_dh->chr_id_2].Value = max( max_chr_length[sp_dh->chr_id_2].Value, sp_dh->end_2 );
    }

    for(const auto &chr_length: max_chr_length)
        chr_regions[chr_length.first].assign(max_chr_length[chr_length.first].Value, 0);

    for(auto iter=dh_list.cbegin(); iter!=dh_list.cend(); iter++)
    {
        for(uLONG index=(*iter)->start_1; index<(*iter)->end_1; index++)
            ++chr_regions[(*iter)->chr_id_1][index-1];
        for(uLONG index=(*iter)->start_2; index<(*iter)->end_2; index++)
            ++chr_regions[(*iter)->chr_id_2][index-1];
    }

    //build_index(interval);
}

uLONG Genome_Covarege::max_cov(const string &chr_id, uLONG start, uLONG end) const
{
    //size_type s_idx = start / interval;
    //size_type e_idx = end / interval;



    try{
        //const vector<Region> &regions = chr_regions.at(chr_id);
        //vector<uLONG> coverage(end-start+1, 0);

       // cout << start << "-" << end << "\t||\t" << start_index.at(chr_id)[s_idx] << "-" << end_index.at(chr_id)[e_idx] << endl;

        /*
        for(auto iter=regions.cbegin()+start_index.at(chr_id)[s_idx]; iter!=regions.cbegin()+end_index.at(chr_id)[e_idx]; iter++)
        {
            if( iter->first <= end and start <= iter->second )
            {
                uLONG lower_border = max(start, iter->first);
                uLONG upper_border = min(end, iter->second);
                for(uLONG index=lower_border; index!=upper_border; index++)
                    ++coverage[ index-start ];
            }
        }
        */

        if(start >= chr_regions.at(chr_id).size())
            return 0;
        end = min(end, chr_regions.at(chr_id).size());

        uLONG cur_max_cov = 0;
        for_each(chr_regions.at(chr_id).cbegin()+start, chr_regions.at(chr_id).cbegin()+end, [&](uLONG value){ if(cur_max_cov < value){ cur_max_cov = value; } });
        return cur_max_cov;

    }catch(out_of_range)
    {
        return 0;
    }
}

/*
void Genome_Covarege::build_index(const uINT interval)
{
    this->interval = interval;
    start_index.clear();
    end_index.clear();

    for(auto iter=chr_regions.cbegin(); iter!=chr_regions.cend(); iter++)
    {
        const string &chr_name = iter->first;
        //uLONG last_start_pos = interval;
        //uLONG last_end_pos = interval;

        start_index[chr_name].push_back(0);

        const vector<Region> &ref_region = iter->second;
        
        for(size_t idx=0; idx<ref_region.size(); idx++)
            while(ref_region[idx].first > (start_index[chr_name].size())*interval and idx<ref_region.size())
                start_index[chr_name].push_back(idx);

        uLONG largest_end = 0;
        for(size_t idx=0; idx<ref_region.size(); idx++)
        {
            largest_end = max(largest_end, ref_region[idx].second);
            while(largest_end > (end_index[chr_name].size()+1)*interval and idx<ref_region.size())
                end_index[chr_name].push_back(idx);
        }

        end_index[chr_name].push_back(ref_region.size());
    }

    //cout << start_index << endl;
    //cout << end_index << endl;

}
*/


bool read_dh(istream &IN, Ext_Duplex_Hang &dh)
{
    string sep;
    IN >> dh.read_id >> dh.chr_id_1 >> dh.strand_1 >> dh.flag_1 >> dh.cigar_1 >> dh.start_1 >> dh.end_1;
    IN >> sep;
    IN >>  dh.chr_id_2 >> dh.strand_2 >> dh.flag_2 >> dh.cigar_2 >> dh.start_2 >> dh.end_2;
    return static_cast<bool>(IN);
}

inline bool duplicate_read_hang(const Ext_Duplex_Hang &dh_1, const Ext_Duplex_Hang &dh_2)
{
    return (dh_1.end_1 != dh_2.end_1 or dh_1.end_2 != dh_2.end_2 or dh_1.start_1 != dh_2.start_1 or dh_1.start_2 != dh_2.start_2 or
        dh_1.cigar_1 != dh_2.cigar_1 or dh_1.cigar_2 != dh_2.cigar_2) ? false : true;
}

void intersect_paris_reads(const string &input_dg, 
                            vector< sp<Duplex_Group> > &dg_list,
                            vector< sp<Ext_Duplex_Hang> > &read_list,
                            const uINT min_overlap,
                            const bool multimapDG,
                            const bool uniqDG)
{
    using size_type = vector<Duplex_Group *>::size_type;

    ifstream IN(input_dg, ifstream::in);
    if(not IN)
    {
        cerr << "FATAL Error: " << input_dg << " is unreadable" << endl;
        exit(-1);
    }

    dg_list.clear();
    size_type firstPossible = 0;

    //vector< sp<Ext_Duplex_Hang> > read_list;

    sp<Ext_Duplex_Hang> sp_dh = make_shared<Ext_Duplex_Hang>();// new Duplex_Hang;
    read_dh(IN, *sp_dh);
    sp<Duplex_Group> sp_dg = make_shared<Duplex_Group>(sp_dh);// new Duplex_Group(sp_dh);
    dg_list.push_back(sp_dg);
    read_list.push_back(sp_dh);
    //sp<Ext_Duplex_Hang> last_dg = 

    sp_dh = make_shared<Ext_Duplex_Hang>();
    while(read_dh(IN, *sp_dh))
    {
        if( uniqDG and duplicate_read_hang(*sp_dh, *read_list.back()) )
            continue;

        read_list.push_back(sp_dh);

        bool lastDGoverlapped(false), nonOverlapped(true);
        for(size_type idx=firstPossible; idx<dg_list.size(); idx++)
        {
            long overlap = dg_list[idx]->check_overlap(*sp_dh);
            if( overlap >= min_overlap )
            {
                nonOverlapped = false;
                lastDGoverlapped = true;
                dg_list[idx]->add_duplex_hang(sp_dh);
                if(not multimapDG)
                    break;
                //else
                //    dh = new Duplex_Hang(*dh);
            }else if(overlap == -1)
            {
                if(not lastDGoverlapped)
                    firstPossible = idx + 1;
            }else{
                lastDGoverlapped = true;
            }
        }
        //if(not nonOverlapped or multimapDG)
        //    delete dh;

        if(nonOverlapped)
        {
            //Duplex_Group* dg = new Duplex_Group(dh);
            dg_list.push_back( make_shared<Duplex_Group>(sp_dh) );
        }
        sp_dh = make_shared<Ext_Duplex_Hang>();
        //dh = new Duplex_Hang;
    }
    //delete dh;
    IN.close();

    //clog << "Sort..." << endl;
    sort(dg_list.begin(), dg_list.end(), [](sp<Duplex_Group> sp_dg_1, sp<Duplex_Group> sp_dg_2){ return *sp_dg_1 < *sp_dg_2; });
/*
    ofstream OUT("/tmp/tmp_dg.txt", ios::out);
    for_each(dg_list.cbegin(), dg_list.cend(), [&](const sp<Duplex_Group> &sp_dg){ OUT << *sp_dg; });
   // OUT << dg_list << endl;
    OUT.close();
*/
    //return dg_list;
}


void collapse_DG(const vector< sp<Duplex_Group> > &dg_list, 
                vector< sp<Duplex_Group> > &dg_array,
                const uINT max_gap, 
                const uINT max_total,
                const bool check_reads)
{
    using size_type = vector< sp<Duplex_Group> >::size_type;

    uINT merged_dg_count = 0;
    size_type firstPossible = 0;
    dg_array.clear();

    dg_array.push_back(dg_list[0]);
    double point = 0.00;
    for(size_type idx=1; idx<dg_list.size(); idx++)
    {
        if(1.0*idx/dg_list.size() > point)
        {
            clog << "\tProcess " << point*100 << "% " << idx << "\t" << "firstPossible: " << firstPossible << endl;
            point += 0.05;
        }

        //if(point > 0.699)
        //    cout << "\t\t" << idx << "\t" << *dg_list[idx] << "\n";

        bool lastDGoverlapped(false);
        long overlap;

        //uLONG total = 
        for(size_type idy=firstPossible; idy<dg_array.size(); idy++)
        {
           // if(point == 0.7 and (idy-firstPossible) % () )

            overlap = dg_array[idy]->check_overlap(*dg_list[idx], max_gap, max_total, check_reads);

            if(overlap == -1)
            {
                if(not lastDGoverlapped)
                    firstPossible = idy + 1;
                lastDGoverlapped = false;
                //dg_array.push_back(dg_list[idx]);
                //break;
            }else if(overlap > 0)
            {
                lastDGoverlapped = true;
                dg_array[idy]->merge_duplex_group(*dg_list[idx]);
                ++merged_dg_count;
                break;
            }
            else{
                lastDGoverlapped = true;
            }
        }
        if(overlap == 0 or overlap == -1)
            dg_array.push_back(dg_list[idx]);
    }

    for(uINT idx=0; idx<dg_array.size(); idx++)
    {
        dg_array[idx]->dg_id = idx;
    }

   // sort(dg_list.begin(), dg_list.end(), [](sp<Duplex_Group> sp_dg_1, sp<Duplex_Group> sp_dg_2){ return *sp_dg_1 < *sp_dg_2; });
    clog << "\tmerged_dg_count: " << merged_dg_count << endl;
}

/*  */
void finalize_reads(vector< sp<Duplex_Group> > &dg_array)
{
    for(auto iter=dg_array.begin(); iter!=dg_array.end(); iter++)
        for(auto read_iter=(*iter)->reads.begin(); read_iter!=(*iter)->reads.end(); read_iter++)
            (*read_iter)->finalize();
}

void print_read(const vector< sp<Duplex_Group> > &dg_array, const string &out_file_name)
{
    ofstream OUT(out_file_name, ios::out);
    if(not OUT)
    {
        cerr << "FATAL ERROR: " << out_file_name << " is unreadable" << endl;
    }
    for(auto iter=dg_array.cbegin(); iter!=dg_array.cend(); iter++)
    {
        for(auto read_iter=(*iter)->reads.cbegin(); read_iter!=(*iter)->reads.cend(); read_iter++)
        {
            if((*read_iter)->chr_id_1 != (*read_iter)->chr_id_2 or (*read_iter)->strand_1 != (*read_iter)->strand_2)
            {
                OUT << (*read_iter)->chr_id_1 << "\t" << (*read_iter)->start_1 << "\t" << (*read_iter)->end_1 
                    << "\t" << (*read_iter)->read_id << "\t1\t" << (*read_iter)->strand_1 << "\n";
                OUT << (*read_iter)->chr_id_2 << "\t" << (*read_iter)->start_2 << "\t" << (*read_iter)->end_2
                    << "\t" << (*read_iter)->read_id << "\t2\t" << (*read_iter)->strand_2 << "\n";
            }else{
                OUT << (*read_iter)->chr_id_1 << "\t" << (*read_iter)->start_1 << "\t" << (*read_iter)->end_2
                    << "\t" << (*read_iter)->read_id << "\t1\t" << (*read_iter)->strand_1 << "\n";
            }
        }
    }
    OUT.close();
}

void build_reads_dg_map(const vector< sp<Duplex_Group> > &dg_array,
                        MapStringT<uLONG> &reads_dg_map,
                        uLONG min_support)
{
    reads_dg_map.clear();
    for(auto iter=dg_array.cbegin(); iter!=dg_array.cend(); iter++)
        if((*iter)->support() >= min_support)
            for(auto read_iter=(*iter)->reads.cbegin(); read_iter!=(*iter)->reads.cend(); read_iter++)
                reads_dg_map[ (*read_iter)->read_id ] = (*read_iter)->dg_ids.front();
}

void tag_sam_dg(const StringArray &input_sam_files,
                const string &output_sam_file,
                const MapStringT<uLONG> &reads_dg_map)
{
    ifstream IN;
    ofstream OUT(output_sam_file, ofstream::out);
    if(not OUT)
    {
        cerr << "FATAL Error: " << output_sam_file << " is unwritable" << endl;
        exit(-1);
    }

    bool success;
    vector<Sam_Record> read_records;
    Sam_Head sam_head;
    
    for(StringArray::size_type idx=0; idx<input_sam_files.size(); idx++)
    {
        IN.open(input_sam_files[idx]);
        if(not IN)
        {
            cerr << "FATAL Error: " << input_sam_files[idx] << " is unreadable" << endl;
            exit(-1);
        }

        if(read_sam_head(IN, sam_head) and idx == 0)
            write_sam_head(OUT, sam_head);

        while(not IN.eof())
        {
            success = read_a_read_record(IN, read_records);
            
            if(not success)
                break;
            try{
                const uLONG dg_id = reads_dg_map.at(read_records.front().read_id);
                for_each(read_records.begin(), read_records.end(), [&](Sam_Record &record){ record.attributes.push_back("DG:i:"+to_string(dg_id)); });
                OUT << read_records;
            }catch(out_of_range)
            {
                continue;
            }
        }
        IN.close();
    }
    OUT.close();
}

void calc_DG_Score(vector< sp<Duplex_Group> > &dg_array, const Genome_Covarege &coverage)
{
    for(sp<Duplex_Group> sp_dg: dg_array)
    {
        //clog << sp_dg->chr_id_1() << "\t" << sp_dg->start_1 << "\t" << sp_dg->end_1 << endl;
        sp_dg->left_cov = coverage.max_cov(sp_dg->chr_id_1(), sp_dg->start_1-1, sp_dg->end_1);
        //clog << sp_dg->chr_id_1() << "\t" << sp_dg->start_1 << "\t" << sp_dg->end_1 << endl;
        sp_dg->right_cov = coverage.max_cov(sp_dg->chr_id_2(), sp_dg->start_2-1, sp_dg->end_2);
        sp_dg->score = 1.0 * sp_dg->support() / 2 * ( 1.0/(sp_dg->left_cov+1) + 1.0/(sp_dg->right_cov+1) );
    }
}

void write_dg( const string &out_file_name,
        const vector< sp<Duplex_Group> > &dg_array,
        const uLONG min_support,
        const string &genome_fasta="")
{
    char group_head[4000];
    char read_line[4000];
    char structure_line[4000];

    ofstream OUT(out_file_name, ofstream::out);

    if(not genome_fasta.empty())
    {
        Fasta genome(genome_fasta);
        
        double counter[] = { 0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95 };
        const vector< sp<Duplex_Group> >::size_type total_number = dg_array.size();
        vector< sp<Duplex_Group> >::size_type cur_count = 0;
        size_t counter_idx = 0;

        for(auto iter=dg_array.cbegin(); iter!=dg_array.cend(); iter++)
        {
            cur_count++;
            if( 1.0 * cur_count / total_number > counter[counter_idx] )
            {
                clog << "\tProcess " << counter[counter_idx] * 100 << "% " << cur_count << endl;
                counter_idx++;
            }

            if((*iter)->support() >= min_support)
            {
                try
                {
                    uLONG len_1 = (*iter)->end_1 - (*iter)->start_1 + 1;
                    string seq_1 = genome.get_chr_subbseq((*iter)->chr_id_1(), (*iter)->start_1-1, len_1, (*iter)->strand_1()=='+' ? POSITIVE : NEGATIVE);
                    uLONG len_2 = (*iter)->end_2 - (*iter)->start_2 + 1;
                    string seq_2 = genome.get_chr_subbseq((*iter)->chr_id_2(), (*iter)->start_2-1, len_2, (*iter)->strand_2()=='+' ? POSITIVE : NEGATIVE);

                    string whole_seq = seq_1 + "III" + seq_2;
                    string whole_structure = fold_two_seq(seq_1, seq_2, false)->at(1);

                    sprintf(group_head, "Group %lu == position %s(%c):%lu-%lu|%s(%c):%lu-%lu, support %lu, left %lu, right %lu, score %.3f.",
                        (*iter)->dg_id, (*iter)->chr_id_1().c_str(), (*iter)->strand_1(), (*iter)->start_1, (*iter)->end_1, 
                                     (*iter)->chr_id_2().c_str(), (*iter)->strand_2(), (*iter)->start_2, (*iter)->end_2, 
                                     (*iter)->support(), (*iter)->left_cov, (*iter)->right_cov, (*iter)->score);
                    sprintf(structure_line, "\n---%s\n---%s\n", whole_seq.c_str(), whole_structure.c_str());
                    OUT << group_head << structure_line;

                    for(auto read_iter=(*iter)->reads.cbegin(); read_iter!=(*iter)->reads.cend(); read_iter++)
                    {
                        sprintf(read_line, "\t%s\t%s|%c:%lu-%lu<=>%s|%c:%lu-%lu", (*read_iter)->read_id.c_str(), 
                            (*read_iter)->chr_id_1.c_str(), (*read_iter)->strand_1, (*read_iter)->start_1, (*read_iter)->end_1,
                            (*read_iter)->chr_id_2.c_str(), (*read_iter)->strand_2, (*read_iter)->start_2, (*read_iter)->end_2 );
                        OUT << read_line << "\n";
                    }
                }catch(out_of_range){
                    cerr << "Warning: " << (*iter)->chr_id_1() << " or " << (*iter)->chr_id_2() << " not in genome file" << endl;
                }
            }
        }
    }else{
        for(auto iter=dg_array.cbegin(); iter!=dg_array.cend(); iter++)
        {
            if((*iter)->support() >= min_support)
            {
                sprintf(group_head, "Group %lu == position %s(%c):%lu-%lu|%s(%c):%lu-%lu, support %lu, left %lu, right %lu, score %.3f.",
                    (*iter)->dg_id, (*iter)->chr_id_1().c_str(), (*iter)->strand_1(), (*iter)->start_1, (*iter)->end_1, 
                                 (*iter)->chr_id_2().c_str(), (*iter)->strand_2(), (*iter)->start_2, (*iter)->end_2, 
                                 (*iter)->support(), (*iter)->left_cov, (*iter)->right_cov, (*iter)->score);
                OUT << group_head << "\n---\n";

                for(auto read_iter=(*iter)->reads.cbegin(); read_iter!=(*iter)->reads.cend(); read_iter++)
                {
                    sprintf(read_line, "\t%s\t%s|%c:%lu-%lu<=>%s|%c:%lu-%lu", (*read_iter)->read_id.c_str(), 
                        (*read_iter)->chr_id_1.c_str(), (*read_iter)->strand_1, (*read_iter)->start_1, (*read_iter)->end_1,
                        (*read_iter)->chr_id_2.c_str(), (*read_iter)->strand_2, (*read_iter)->start_2, (*read_iter)->end_2 );
                    OUT << read_line << "\n";
                }
            }
        }
        
    }
    OUT.close();
}

void write_tab( const string &out_file_name, 
                const vector< sp<Duplex_Group> > &dg_array,
                const uLONG min_support)
{
    ofstream OUT(out_file_name, ofstream::out);

    OUT << "#" << "\t" << "dg_id\t" << "chr_id_1\t" << "strand_1\t" << "start_1\t" << "end_1\t" << "chr_id_2\t" << "strand_2\t" << "start_2\t" << "end_2\t"
        << "support\t" << "left_cov\t" << "right_cov\t" << "score\n";

    for(auto iter=dg_array.cbegin(); iter!=dg_array.cend(); iter++)
    {
        if((*iter)->support() < min_support)
            continue;

        OUT << ">" << "\t" << (*iter)->dg_id << "\t" << (*iter)->chr_id_1() << "\t" << (*iter)->strand_1() << "\t" << (*iter)->start_1 << "\t" << (*iter)->end_1
                                             << "\t" << (*iter)->chr_id_2() << "\t" << (*iter)->strand_2() << "\t" << (*iter)->start_2 << "\t" << (*iter)->end_2 
                                             << "\t" << (*iter)->support() << "\t" << (*iter)->left_cov << "\t" << (*iter)->right_cov << "\t" << (*iter)->score << "\n";

        for(auto read_iter=(*iter)->reads.cbegin(); read_iter!=(*iter)->reads.cend(); read_iter++)
            OUT << (*read_iter)->read_id << "\t" << (*read_iter)->chr_id_1 << "\t" << (*read_iter)->strand_1 << "\t" << (*read_iter)->start_1 << "\t" << (*read_iter)->end_1
                                         << "\t" << (*read_iter)->chr_id_2 << "\t" << (*read_iter)->strand_2 << "\t" << (*read_iter)->start_2 << "\t" << (*read_iter)->end_2 << "\n";
    }

    OUT.close();
}

void check_dg(const vector< sp<Duplex_Group> > &dgs)
{
    for(const sp<Duplex_Group> sp_dg: dgs)
    {
        if(sp_dg->start_1 > sp_dg->end_1 or sp_dg->start_2 > sp_dg->end_2)
        {
            cerr << *sp_dg << endl;
            for(auto iter=sp_dg->reads.cbegin(); iter!=sp_dg->reads.cend(); iter++)
                cerr << **iter;
            exit(-1);
        }
    }
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

    vector< sp<Duplex_Group> > dg_list, dg_array;
    vector< sp<Ext_Duplex_Hang> > read_list;
    MapStringT<uLONG> reads_dg_map;

    clog << "start to intersect_paris_reads..." << endl;
    intersect_paris_reads(param.input_dg, dg_list, read_list, param.min_overlap, param.multiDG, param.uniqDG);

    //check_dg(dg_list);

    clog << "start to build genome coverage index..." << endl;
    Genome_Covarege coverage( read_list );

    clog << "start to collapse_DG..." << endl;
    collapse_DG(dg_list, dg_array, param.max_gap, param.max_total, param.check_reads);
    
    //check_dg(dg_array);
    while(dg_list != dg_array)
    {
        //sort(dg_array.begin(), dg_array.end(), [](sp<Duplex_Group> sp_dg_1, sp<Duplex_Group> sp_dg_2){ return *sp_dg_1 < *sp_dg_2; });
        dg_list.clear();
        dg_list = dg_array;
        collapse_DG(dg_list, dg_array, param.max_gap, param.max_total, param.check_reads);
        //check_dg(dg_array);
    }
    dg_list.clear();
    //sort(dg_array.begin(), dg_array.end(), [](sp<Duplex_Group> sp_dg_1, sp<Duplex_Group> sp_dg_2){ return *sp_dg_1 < *sp_dg_2; });

    clog << "start to finalize_reads..." << endl;
    finalize_reads(dg_array);

    /*
    clog << "start to print_read..." << endl;
    print_read(dg_array, "/tmp/59.reads");
    */
    
    clog << "start to calc_DG_Score..." << endl;
    calc_DG_Score(dg_array, coverage);

    clog << "start to write_dg..." << endl;
    write_dg(param.output_dg, dg_array, param.min_support, param.genome_fasta);

    if(not param.output_tab.empty())
    {
        clog << "start to write_tab..." << endl;
        write_tab( param.output_tab, dg_array, param.min_support);
    }
    
    if(not param.input_tag_sam.empty())
    {
        clog << "start to build_reads_dg_map..." << endl;
        build_reads_dg_map(dg_array, reads_dg_map, param.min_support);
        
        clog << "start to tag_sam_dg..." << endl;
        tag_sam_dg(param.input_tag_sam, param.output_tag_sam, reads_dg_map);
    }

    return 0;
}












