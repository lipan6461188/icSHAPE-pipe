
#include "paris.h"
#include "exceptions.h"

using namespace std;

namespace pan{

/*

    string read_id;

    string chr_id_1;
    char strand_1;
    uINT flag_1;
    string cigar_1;
    uLONG start_1;
    uLONG end_1;
    
    string chr_id_2;
    char strand_2;
    uINT flag_2;
    string cigar_2;
    uLONG start_2;
    uLONG end_2;

*/


ostream &operator<<(ostream &OUT, const Duplex_Hang &dh)
{
    //OUT << dh.read_id << "\t";

    OUT << dh.read_id << "\t"
        << dh.chr_id_1 << "\t"
        << dh.strand_1 << "\t" 
        << dh.flag_1 << "\t" 
        << dh.cigar_1 << "\t" 
        << dh.start_1 << "\t" 
        << dh.end_1 << "\t"
        << "||" << "\t"
        << dh.chr_id_2 << "\t"
        << dh.strand_2 << "\t" 
        << dh.flag_2 << "\t" 
        << dh.cigar_2 << "\t" 
        << dh.start_2 << "\t" 
        << dh.end_2 << "\n";

    return OUT;
}

ostream &operator<<(ostream &OUT, const vector<Duplex_Hang> &dh_array)
{
    for(const auto &dh: dh_array)
        OUT << dh;

    return OUT;
}

bool operator<(const Duplex_Hang &dh_1, const Duplex_Hang &dh_2)
{
    // chr 1
    if(dh_1.chr_id_1 < dh_2.chr_id_1)
        return true;
    else if(dh_1.chr_id_1 > dh_2.chr_id_1)
        return false;
    else{
        // chr 2
        if(dh_1.chr_id_2 < dh_2.chr_id_2)
            return true;
        else if(dh_1.chr_id_2 > dh_2.chr_id_2)
            return false;
        else{
            // strand 1
            if(dh_1.strand_1 < dh_2.strand_1)
                return true;
            else if(dh_1.strand_1 > dh_2.strand_1)
                return false;
            else{
                // strand 2
                if(dh_1.strand_2 < dh_2.strand_2)
                    return true;
                else if(dh_1.strand_2 > dh_2.strand_2)
                    return false;
                else{
                    // start 1
                    if(dh_1.start_1 < dh_2.start_1)
                        return true;
                    else if(dh_1.start_1 > dh_2.start_1)
                        return false;
                    else{
                        // start 2
                        if(dh_1.start_2 < dh_2.start_2)
                            return true;
                        else if(dh_1.start_2 > dh_2.start_2)
                            return false;
                        else
                            return false;
                    }
                }
            }
        }
    }
}


/* return chromosome length */
uLONG get_chromosome_hang(const string &sam_file_name, 
                            vector<Duplex_Hang> &dh_array, 
                            uINT min_gap,
                            uINT min_hang,
                            const string &chr_id,
                            const char strand)
{
    ifstream IN(sam_file_name, ifstream::in);
    if(not IN)
    {
        throw runtime_error( "Bad_Input_File: "+sam_file_name );
    }

    dh_array.clear();
    
    Sam_Head sam_head;
    read_sam_head(IN, sam_head);
    if(sam_head.trans_len.find(chr_id) == sam_head.trans_len.end())
    {
        throw runtime_error(chr_id+" not find in "+sam_file_name);
    }

    uLONG valid_count = 0;

    Sam_Record read_record;
    while(not IN.eof())
    {
        bool success = read_a_sam_record(IN, read_record);
        if(not success)
            break;
        if(read_record.chr_id != chr_id)
            continue;
        if(read_record.strand() != strand)
            continue;
        RegionArray matchRegion;
        get_global_match_region(read_record.cigar, read_record.pos, matchRegion);
        if( matchRegion.size() == 2 )
        {
            auto gap = matchRegion.at(1).first - matchRegion.at(0).second - 1;
            auto left_hang = matchRegion.at(0).second - matchRegion.at(0).first + 1;
            auto right_hang = matchRegion.at(0).second - matchRegion.at(0).first + 1;
            if( gap >= min_gap and left_hang >= min_hang and right_hang >= min_hang )
            {
                Duplex_Hang dh( read_record.read_id, 

                                read_record.chr_id, 
                                read_is_reverse(read_record)?'-':'+', 
                                read_record.flag,
                                read_record.cigar,
                                matchRegion.at(0).first,
                                matchRegion.at(0).second,

                                read_record.chr_id, 
                                read_is_reverse(read_record)?'-':'+', 
                                read_record.flag,
                                read_record.cigar,
                                matchRegion.at(1).first,
                                matchRegion.at(1).second);
                dh_array.push_back(dh);
                valid_count++;
            }
        }
    }

/*

 Duplex_Hang(const string &read_id, 
        const string &chr_id_1, const char &strand_1, const uINT &flag_1, const string &cigar_1, const uLONG &start_1, const uLONG &end_1,
        const string &chr_id_2, const char &strand_2, const uINT &flag_2, const string &cigar_2, const uLONG &start_2, const uLONG &end_2)
*/


/*

    string read_id;

    string chr_id_1;
    char strand_1;
    uINT flag_1;
    string cigar_1;
    uLONG start_1;
    uLONG end_1;
    
    string chr_id_2;
    char strand_2;
    uINT flag_2;
    string cigar_2;
    uLONG start_2;
    uLONG end_2;

*/


    IN.close();

    sort(dh_array.begin(), dh_array.end());
/*
    sort(dh_array.begin(), dh_array.end(), [](const Duplex_Hang&d_1, const Duplex_Hang&d_2)->bool{
            if(d_1.chr_id < d_2.chr_id){
                return true;
            }
            else if(d_1.chr_id == d_2.chr_id)
            { 
                return d_1.start_1<d_2.start_1 ? true : false;
            }
            else{ 
                return false; 
            }
        });
*/

    //clog << "Useful Reads Numbers: " << valid_count << endl;
    return sam_head.trans_len.at(chr_id);
}

void read_dh_from_sam(const string &sam_file_name, vector<Duplex_Hang> &dh_array)
{
    ifstream IN(sam_file_name, ifstream::in);
    if(not IN)
    {
        throw Bad_IO( "Bad_Input_File: "+sam_file_name, true);
    }

    dh_array.clear();
    
    Sam_Head sam_head;
    read_sam_head(IN, sam_head);
    
    /*
    if(sam_head.trans_len.find(chr_id) == sam_head.trans_len.end())
    {
        throw Unexpected_Error(chr_id+" not find in "+sam_file_name, true);
    }
    */

    uLONG valid_count = 0;

    vector<Sam_Record> read_records;
    while(not IN.eof())
    {
        bool success = read_a_read_record(IN, read_records);
        if(not success)
            break;
        if(read_records.size() == 1)
        {
            auto &read_record = read_records[0];
            RegionArray matchRegion;
            get_global_match_region(read_record.cigar, read_record.pos, matchRegion);
            if( matchRegion.size() == 2 )
            {
                Duplex_Hang dh( read_record.read_id, 

                                read_record.chr_id, 
                                read_is_reverse(read_record)?'-':'+', 
                                read_record.flag,
                                read_record.cigar,
                                matchRegion.at(0).first,
                                matchRegion.at(0).second,

                                read_record.chr_id, 
                                read_is_reverse(read_record)?'-':'+', 
                                read_record.flag,
                                read_record.cigar,
                                matchRegion.at(1).first,
                                matchRegion.at(1).second);
                dh_array.push_back(dh);
                valid_count++;
            }
        }else if(read_records.size() == 2)
        {
            RegionArray matchRegion_1, matchRegion_2;
            get_global_match_region(read_records[0].cigar, read_records[0].pos, matchRegion_1);
            get_global_match_region(read_records[1].cigar, read_records[1].pos, matchRegion_2);
            
            if( matchRegion_1.size() == 1 and matchRegion_2.size() == 1 )
            {
                Duplex_Hang dh( read_records[0].read_id, 

                                read_records[0].chr_id, 
                                read_is_reverse(read_records[0])?'-':'+', 
                                read_records[0].flag,
                                read_records[0].cigar,
                                matchRegion_1.at(0).first,
                                matchRegion_1.at(0).second,

                                read_records[1].chr_id,
                                read_is_reverse(read_records[1])?'-':'+', 
                                read_records[1].flag,
                                read_records[1].cigar,
                                matchRegion_2.at(1).first,
                                matchRegion_2.at(1).second);
                dh_array.push_back(dh);
                valid_count++;
            }
        }else{
            // bad input
        }
    }

    IN.close();
    sort(dh_array.begin(), dh_array.end());
    //return sam_head.trans_len.at(chr_id);
}




ostream& operator<<(ostream& OUT, const InterRegion& iter_reg)
{
    OUT << iter_reg.region.first.first << "-" << iter_reg.region.first.second << "\t" 
        << iter_reg.region.second.first << "-" << iter_reg.region.second.second << "\t" 
        << iter_reg.max_cov << "\t" << iter_reg.max_point << "\n";
    return OUT;
}

ostream& operator<<(ostream& OUT, const vector<InterRegion>& iter_reg)
{
    for(const InterRegion& reg: iter_reg)
        OUT << reg;
    return OUT;
}

bool operator<(const InterRegion &inter_1, const InterRegion &inter_2)
{
    if(inter_1.region.first < inter_2.region.first)
        return true;
    else if(inter_1.region.first == inter_2.region.first and inter_1.region.second < inter_2.region.second)
        return true;
    return false;
}

pair<uLONG, uLONG> overlap(const InterRegion &iter_1, const InterRegion &iter_2)
{
    uLONG over_left = overlap(iter_1.region.first, iter_2.region.first);
    uLONG over_right = overlap(iter_1.region.second, iter_2.region.second);
    return make_pair(over_left, over_right);
}
/*

struct InterRegion {
    InterRegion(Region region, double max_cov): region(region),max_cov(max_cov){}

    Region region;
    uLONG max_cov = 0;
};

*/



vector<Point> search_TT_cross_linking(const SStructure &structure)
{
    vector<Point> cross_link_points;

    const string &sequence = structure.getSequence();
    const BpArray &bps = structure.get_base_pairs();

    for(auto iter=bps.begin(); iter!=bps.cend(); iter++)
    {
        if(iter->left < iter->right)
        {
            char top_left_base('.'), top_center_base('.'), top_right_base('.');
            char bot_left_base('.'), bot_center_base('.'), bot_right_base('.');

            if(iter->right != sequence.size()-1)
                top_left_base = sequence[iter->right];
            top_center_base = sequence[iter->right-1];
            top_right_base = sequence[iter->right-2];

            if(iter->left != 0)
                bot_left_base = sequence[iter->left-2];
            bot_center_base = sequence[iter->left-1];
            bot_right_base = sequence[iter->left];

            bool top_U = (top_center_base == 'T' or top_center_base == 'U');
            bool bot_U = (bot_center_base == 'T' or bot_center_base == 'U');

            bool top_left_flanking_U = (top_left_base == 'T' or top_left_base == 'U');
            bool top_right_flanking_U = (top_right_base == 'T' or top_right_base == 'U');
            
            bool bot_left_flanking_U = (bot_left_base == 'T' or bot_left_base == 'U');
            bool bot_right_flanking_U = (bot_right_base == 'T' or bot_right_base == 'U');

            if( top_U and bot_left_flanking_U )
                cross_link_points.push_back(Point(iter->left-1, iter->right));
            
            if( top_U and bot_right_flanking_U )
                cross_link_points.push_back(Point(iter->left+1, iter->right));
            
            if( bot_U and top_left_flanking_U )
                cross_link_points.push_back(Point(iter->left, iter->right+1));
            
            if( bot_U and top_right_flanking_U )
                cross_link_points.push_back(Point(iter->left, iter->right-1));   
        }
    }

    sort(cross_link_points.begin(), cross_link_points.end(), [](const Point &p_1, const Point &p_2){ return p_1.first < p_2.first; });
    return cross_link_points;
}


void compress_regions(const RegionArray &raw_regions,
                    RegionArray &target_regions,
                    const uLONG raw_size,
                    const uLONG target_size)
{
    target_regions.clear();

    auto round = [](double raw)->uLONG{ return (raw-floor(raw) >= 0.5) ? ceil(raw) : floor(raw); };

    if(sorted_no_overlap(raw_regions))
    {
        double step = 1.0 * raw_size / target_size;
        //cout << "raw: " << raw_regions << endl;
        for_each(raw_regions.cbegin(), raw_regions.cend(), [&](const Region &region){ target_regions.push_back( Region( round(region.first/step), round(region.second/step) ) ); });
        //cout << "tranfor" << target_regions << endl;
        target_regions.front().first = max(target_regions.front().first, 1UL);
        target_regions.back().second = min(target_size, target_regions.back().second);
        //cout << "tranfor" << target_regions << endl;

        for(auto iter=target_regions.begin()+1; iter!=target_regions.end(); )
        {
            if( (iter-1)->second == iter->first )
            {
                if( iter->first == iter->second )
                {
                    iter = target_regions.erase(iter);
                    continue;
                }
                else{
                    iter->first++;
                }
            }
            iter++;
        }
    }else{
        ostringstream OUT;
        OUT << raw_regions << "\n";
        OUT << "Error: unsorted or overlapped regions in compress_regions\n";
        throw Unexpected_Error(OUT.str());
    }
}


void read_domain_file(const string &domain_file_name, RegionArray &regions)
{
    regions.clear();
    ifstream IN(domain_file_name, ifstream::in);
    
    string cur_line;
    while(getline(IN, cur_line))
    {
        StringArray items;
        split(cur_line, items);
        if(items.size() == 1)
        {
            continue;
        }else if(items.size() == 2)
        {
            regions.push_back( Region(stoul(items[0]), stoul(items[1])) );
        }else if(items[2] != "skip")
        {
            regions.push_back( Region(stoul(items[0]), stoul(items[1])) );
        }
    }
    IN.close();
    sort(regions.begin(), regions.end());
}




}