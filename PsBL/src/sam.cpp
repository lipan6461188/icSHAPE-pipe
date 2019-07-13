
#include "sam.h"
#include "fasta.h"

using namespace std;
//using namespace pan;

namespace pan{

char buffer_1[500], buffer_2[500];

bool read_is_mapped(const Sam_Record &read_record)
{
    if( read_record.flag & 4 )
    {
        //printRead(read);
        return false;
    }
    else
        return true;
}

bool read_is_gapped(const Sam_Record &read_record)
{
    if( read_record.cigar.find('N') != string::npos )
    {
        //printRead(read);
        return true;
    }
    else
        return false;
}

bool read_is_reverse(const Sam_Record &read_record)
{
    if ( read_record.flag & 16 )
    {
        //printRead(read);
        return true;
    }
    else
        return false;
}

bool read_is_primary(const Sam_Record &read_record)
{
    if ( read_record.flag & 256 or read_record.flag & 2048 )
    {
        //printRead(read);
        return false;
    }
    else
        return true;
}

bool read_is_multimap(const vector<Sam_Record> &read_records)
{
    if(read_records.size() == 0)
        return false;
    else if( not read_is_paired(read_records[0]) )
    {
        return read_records.size() > 1;
    }else{
        //if(read_records.size() % 2 == 1) //奇数
        //    throw Unexpected_Error("paired-end reads shouldn't be odd number");
        if(read_records.size() > 2)
            return true;
        else if( is_paired_reads(read_records[0], read_records[1]) )
            return false;
        else
            return true;
    }
}


string Sam_Record::get_attr(const string &attr_name)
{
    uINT attr_name_len = attr_name.size();
    for(const string &attr_item: attributes)
    {
        if( attr_item.substr(0, attr_name_len) == attr_name )
        {
            return split(attr_item, ':').at(2);
        }
    }
    return "";
}

bool read_a_sam_record(istream &IN, Sam_Record &read_record)
{
    // clear read_record
    read_record.attributes.clear();

    string cur_line;
    //cout << "read_a_sam_record..." << endl;
    while(IN and cur_line.empty())
        getline(IN, cur_line);
    
    if(not IN or IN.eof())
    {
        return false;
    }
    //cout << cur_line << endl;
    if(cur_line.at(0) == '@')
    {
        return false;
    }
    istringstream string_in(cur_line);
    string_in >> read_record.read_id >> read_record.flag >> read_record.chr_id >> read_record.pos >> 
        read_record.map_quanlity >> read_record.cigar >> read_record.read_id_next >> read_record.pos_next >> 
        read_record.temp_len >> read_record.read_seq >> read_record.read_quality;

    string cur_attributes;
    while( string_in.good() )
    {
    	//cout << cur_attributes << endl;;
        string_in >> cur_attributes;
        read_record.attributes.push_back(cur_attributes);
    }

    return true;
}



void write_a_sam_record(ostream &OUT, const Sam_Record &read_record)
{
    OUT << read_record.read_id << "\t" << read_record.flag << "\t" << read_record.chr_id << "\t" << 
        read_record.pos << "\t" << read_record.map_quanlity << "\t" << read_record.cigar << "\t" << 
        read_record.read_id_next << "\t" << read_record.pos_next << "\t" << read_record.temp_len << "\t" << 
        read_record.read_seq << "\t" << read_record.read_quality;

    for(auto iter=read_record.attributes.cbegin(); iter!=read_record.attributes.cend(); iter++)
        OUT << "\t" << *iter;
    OUT << "\n";
}

bool read_sam_head(istream &IN, Sam_Head &sam_head)
{
    //clear sam_head
    sam_head.trans_len.clear();
    sam_head.head_list.clear();

    if(not IN)
    {
        return false;
    }
    auto last_pos = IN.tellg();
    string cur_line;
    while(getline(IN, cur_line))
    {
        if(cur_line.at(0) != '@')
        {
        	//cout << cur_line << endl;
            IN.seekg(last_pos);
            break;
        }
        last_pos = IN.tellg();
        if( cur_line.substr(0, 3) == "@SQ" )
        {
            string cur_item;
            string trans_id;
            uLONG trans_len = 0;
            istringstream string_in(cur_line);
            while( string_in.good() )
            {
                string_in >> cur_item;
                if(cur_item.substr(0, 3) == "SN:")
                    trans_id = cur_item.substr(3);
                else if(cur_item.substr(0, 3) == "LN:")
                    trans_len = stoul(cur_item.substr(3));
            }
            sam_head.trans_len[ trans_id ] = trans_len;
        }else{
            sam_head.head_list.push_back( cur_line );
        }
    }
    return true;
}


bool read_sam_head(bam_hdr_t *hdr, Sam_Head &sam_head)
{
    getBamHead(hdr, sam_head.trans_len);
    for(auto it=sam_head.trans_len.cbegin(); it!=sam_head.trans_len.cend(); it++)
        sam_head.head_list.push_back(it->first);

    return true;
}

bool read_a_sam_record(BGZF* fn_hd, bam_hdr_t *hdr, Sam_Record &read_record)
{
    // clear read_record
    read_record.attributes.clear();

    bam1_t *record = bam_init1();
    int ret = bam_read1(fn_hd, record);
    if(ret < 0)
        return false;

    read_record.read_id = getBamQName(record);
    read_record.flag = getBamFlag(record);
    if(read_record.flag & 4)
    {
       // unmap
        read_record.chr_id = "*";
        read_record.pos = 0;
        read_record.map_quanlity = 0;
    }else{
        read_record.chr_id = getBamRef(record, hdr);
        read_record.pos = getBamRefPos(record);
        read_record.map_quanlity = getBamMapQuanlity(record);
        read_record.cigar = getBamCigar(record);
        read_record.read_id_next = getBamMateRef(record, hdr);
        read_record.pos_next = getBamMateRefPos(record);
        read_record.temp_len = 0;
        read_record.read_seq = getBamSeq(record);
        read_record.read_quality = getBamQuanlity(record);

        string cur_attributes;
        stringstream string_in(getBamTag(record));
        while( string_in.good() )
        {
            string_in >> cur_attributes;
            read_record.attributes.push_back(cur_attributes);
        }
    }

    bam_destroy1(record);
    return true;
}

bool read_a_read_record(BGZF* fn_hd, bam_hdr_t *hdr, vector<Sam_Record> &read_records, Sam_Record* &p_cache)
{
    // clear read_records
    read_records.clear();

    if(p_cache)
    {
        read_records.push_back( *p_cache );
        delete p_cache;
        p_cache = nullptr;
    }

    Sam_Record read_record;
    while( read_a_sam_record(fn_hd, hdr, read_record) )
    {
        if(read_records.size() == 0)
        {
            read_records.push_back(read_record);
        }else{
            if( read_records.front().read_id != read_record.read_id )
            {
                p_cache = new Sam_Record;
                *p_cache = read_record;
                break;
            }
            else{
                read_records.push_back(read_record);
            }
        }
    }

    if(read_records.size() == 0)
        return false;

    return true;
}


void write_sam_head(ostream &OUT, const Sam_Head &sam_head)
{
	for(auto iter=sam_head.head_list.cbegin(); iter!=sam_head.head_list.cend(); iter++)
	{
		OUT << *iter << "\n";
	}
	for(auto iter=sam_head.trans_len.cbegin(); iter!=sam_head.trans_len.cend(); iter++)
	{
		OUT << "@SQ\tSN:" << iter->first << "\tLN:" << iter->second << "\n";
	}
}

ostream& operator<<(ostream &OUT, const Sam_Head &read_records)
{
    write_sam_head(OUT, read_records);
    return OUT;
}

bool operator<(const Sam_Record& record_1, const Sam_Record& record_2)
{
    if(record_1.chr_id < record_2.chr_id)
        return true;
    else if(record_1.chr_id > record_2.chr_id)
        return false;
    else{
        if(record_1.pos < record_2.pos)
            return true;
        else if(record_1.pos > record_2.pos)
            return false;
        else
            return false;
    }
}

bool read_a_read_record(istream &IN, vector<Sam_Record> &read_records)
{
    // clear read_records
    read_records.clear();

    if(not IN)
    {
        return false;
    }
    auto last_pos = IN.tellg();

    Sam_Record read_record;
    //bool success = true;

    bool success = true;
    while( not IN.eof() )
    {

        success = read_a_sam_record(IN, read_record);
        /*
        try{
            success = read_a_sam_record(IN, read_record);
        }catch(std::out_of_range){
            cerr << "read_a_sam_record Error in read_a_read_record..." << endl;
            exit(0);
        }*/
    
        if(not success)
            break;

        if(read_records.size() == 0)
        {
            read_records.push_back(read_record);
        }else{
            if( read_records.front().read_id != read_record.read_id )
            {
                IN.seekg(last_pos);
                break;
            }
            else{
                read_records.push_back(read_record);
            }
        }
        last_pos = IN.tellg();
    }
    if(read_records.size() == 0)
        return false;
    return true;
}

void write_a_read_record(ostream &OUT, const vector<Sam_Record> &read_records)
{
    for(const Sam_Record &read_record: read_records)
    {
        write_a_sam_record(OUT, read_record);
    }
}

void remove_read_mutimap(Sam_Record &read_record)
{
    read_record.flag &= ~256UL;
}

void add_read_mutimap(Sam_Record &read_record)
{
    read_record.flag |= 256UL;
}

void remove_read_reversed(Sam_Record &read_record)
{
    read_record.flag &= ~16UL;
}

void add_read_reversed(Sam_Record &read_record)
{
    read_record.flag |= 16UL;
}

ostream& operator<<(ostream &OUT, const Sam_Record &read_record)
{
    write_a_sam_record(OUT, read_record);
    return OUT;
}

ostream& operator<<(ostream &OUT, const vector<Sam_Record> &read_records)
{
    write_a_read_record(OUT, read_records);
    return OUT;
}

ostream& operator<<(ostream &OUT, const ReadPair &read_pair)
{
    OUT << read_pair.first << read_pair.second;
    return OUT;
}

void split_cigar(const string &cigar, vector<int> &cigarLen, vector<char> &cigarAlpha)
{
    string lastNum;
    cigarLen.clear();
    cigarAlpha.clear();
    for(size_t i=0; i<cigar.size(); i++)
    {
        char character = cigar.at(i);
        if(character >= '0' and character <= '9')
        {
            lastNum.push_back(character);
        }else{
            if(not lastNum.empty())
            {
                cigarLen.push_back( stoi(lastNum) );
                lastNum.clear();
            }
            cigarAlpha.push_back(character);
        }
    }
}

/*
 string cigar;
 while(cin >> cigar)
 {
 vector<pair<uLONG, uLONG>> matchRegion;
 getMatchRegion(cigar, 1, matchRegion);
 for(auto region: matchRegion)
 cout << region.first << "-" << region.second << endl;
 cout << "\n";
 }
 */

/*
 *  M 1 match/mismatch
 *  D 1 deletion(very short skip)
 *  I 0 insertion
 *  N 0 1 skipped(genome region is included)
 *  S 0 soft-clip
 *  H 0 hard-clip
 *  P 0 padding
 *  = 1 match
 *  X 1 mismatch
 *
 *  We can test these cigars:
 *  1M2N2M 1M2N2N2M 1M2N2N2M2M 1M2N2N2M2M10S 10S1M2N2N2M2M10S 10S1M2N1X1=2N2M2M10S
 */

void get_global_match_region(const string &cigar, 
                            uLONG startPos, 
                            RegionArray &matchRegion)
{
    vector<int> cigarLen;
    vector<char> cigarAlpha;
    matchRegion.clear();
    split_cigar(cigar, cigarLen, cigarAlpha);

    uLONG lastStartPos = startPos;
    uLONG curGenomePos = startPos;
    for(size_t i=0; i<cigarLen.size(); ++i)
    {
        switch(cigarAlpha[i])
        {
            case 'M': case 'D': case '=': case 'X':
                curGenomePos += cigarLen[i];
                break;
            case 'I': case 'S': case 'H':
                if(lastStartPos != curGenomePos)
                {
                    matchRegion.push_back(make_pair(lastStartPos, curGenomePos-1));
                    lastStartPos = curGenomePos;
                }
                break;
            case 'P':
                break;
            case 'N':
                if(lastStartPos != curGenomePos)
                    matchRegion.push_back(make_pair(lastStartPos, curGenomePos-1));
                lastStartPos = curGenomePos = curGenomePos + cigarLen[i];
                break;
            default:
                cerr << "Unrecognized Cigar Alpha: " << cigar << endl;
                matchRegion.clear();
                return;
        }
    }

    if(lastStartPos != curGenomePos)
        matchRegion.push_back(make_pair(lastStartPos, curGenomePos-1));
}

string::size_type remove_left_D(string &raw_cigar)
{
    vector<int> cigarLen;
    vector<char> cigarAlpha;
    split_cigar(raw_cigar, cigarLen, cigarAlpha);

    if(cigarAlpha.at(0) == 'D')
    {
        raw_cigar = "";
        for(string::size_type idx=1; idx<cigarLen.size(); idx++)
            raw_cigar += to_string(cigarLen.at(idx)) + cigarAlpha.at(idx);
        return cigarLen.at(0);
    }else{
        return 0;
    }
}

string::size_type remove_right_D(string &raw_cigar)
{
    vector<int> cigarLen;
    vector<char> cigarAlpha;
    split_cigar(raw_cigar, cigarLen, cigarAlpha);

    if(cigarAlpha.back() == 'D')
    {
        raw_cigar = "";
        for(string::size_type idx=0; idx<cigarLen.size()-1; idx++)
            raw_cigar += to_string(cigarLen.at(idx)) + cigarAlpha.at(idx);
        return cigarLen.back();
    }else{
        return 0;
    }
}


/*
    Cigar:

ref     √   x   √   x
read    √   √   x   x
Cigar   M   I   D   P
        X   S   N
            H

*/

void get_local_match_region(const string &cigar, 
                            RegionArray &matchRegion)
{
    vector<int> cigarLen;
    vector<char> cigarAlpha;
    matchRegion.clear();
    split_cigar(cigar, cigarLen, cigarAlpha);

   // cout << cigarLen;
    //cout << cigarAlpha;

    uLONG lastStartPos = 1;
    uLONG curGenomePos = 1;
    for(size_t i=0; i<cigarLen.size(); ++i)
    {
        switch(cigarAlpha[i])
        {
            case 'M': case '=': case 'X':
                curGenomePos += cigarLen[i];
                break;
            case 'I': case 'S': case 'H':
                if(lastStartPos != curGenomePos)
                    matchRegion.push_back(make_pair(lastStartPos, curGenomePos-1));
                lastStartPos = curGenomePos = curGenomePos + cigarLen[i];
                break;
            case 'D': case 'N':
                break;
            case 'P':
                break;
            default:
                throw Unexpected_Error("Unrecognized Cigar Alpha: "+cigar);
                //cerr << "Unrecognized Cigar Alpha: " << cigar << endl;
                //matchRegion.clear();
                //return;
        }
    }

    if(lastStartPos != curGenomePos)
        matchRegion.push_back(make_pair(lastStartPos, curGenomePos-1));
}

// reverse cigar code
string reverse_cigar(const string &cigar)
{
    vector<int> cigarLen;
    vector<char> cigarAlpha;
    string reversed_cigar;

    split_cigar(cigar, cigarLen, cigarAlpha);

    for(int idx=cigarLen.size()-1; idx>=0; idx--)
    {
        reversed_cigar += to_string(cigarLen[idx]) + cigarAlpha[idx];
    }

    return reversed_cigar;
}



// get the attributes from a sam record
uINT read_max_softclip(const string &cigar)
{
    vector<int> cigarLen;
    vector<char> cigarAlpha;
    split_cigar(cigar, cigarLen, cigarAlpha);

    uINT max_soft(0);
    
    if(cigarAlpha.front() == 'S')
        max_soft = cigarLen.front();

    if(cigarAlpha.back() == 'S')
        max_soft = std::max((uINT)cigarLen.back(), max_soft);

    return max_soft;
}



MATCH_PATTERN cigar_match_pattern(const string &raw_cigar)
{
    vector<int> cigarLen;
    vector<char> cigarAlpha;
    split_cigar(raw_cigar, cigarLen, cigarAlpha);

    //MATCH_PATTERN match_pattern;
    if(cigarAlpha.front() == 'S' and cigarAlpha.back() != 'S')
        return MATCH_RIGHT;
    else if(cigarAlpha.front() != 'S' and cigarAlpha.back() == 'S')
        return MATCH_LEFT;
    else if(cigarAlpha.front() != 'S' and cigarAlpha.back() != 'S')
        return MATCH_TOTAL;
    else{
        if(cigarLen.front() == cigarLen.back()) 
            return MATCH_CENTER;
        else
            return (cigarLen.front() > cigarLen.back()) ? MATCH_RIGHT : MATCH_LEFT;
    }
}

/** 
 * \param: raw_cigar raw cigar code
 * \return: cigar code of matched part
 */
string get_match_cigar(const string &raw_cigar)
{
    vector<int> cigarLen;
    vector<char> cigarAlpha;

    split_cigar(raw_cigar, cigarLen, cigarAlpha);

    MATCH_PATTERN match_pattern = cigar_match_pattern(raw_cigar);
    if(match_pattern == MATCH_TOTAL)
        return raw_cigar;

    string match_cigar;
    if(match_pattern == MATCH_LEFT)
    {
        for(size_t idx=0; idx<cigarLen.size()-1; idx++)
            match_cigar += to_string(cigarLen[idx]) + cigarAlpha[idx];
        return match_cigar;
    }else{
        for(size_t idx=1; idx<cigarLen.size(); idx++)
            match_cigar += to_string(cigarLen[idx]) + cigarAlpha[idx];
        return match_cigar;
    }

    return match_cigar;
}


/** 
 * \param: raw_cigar raw cigar code
 * \return: cigar code of given region
 */


string get_region_cigar(const string &raw_cigar, const Region &region)
{
    string match_cigar;

    vector<int> cigarLen;
    vector<char> cigarAlpha;
    split_cigar(raw_cigar, cigarLen, cigarAlpha);

    uLONG last_pos = 1;
    for(size_t i=0; i<cigarLen.size(); ++i)
    {
        bool in_region = (last_pos >= region.first and last_pos <= region.second);
        if(not in_region and (cigarAlpha[i] == 'D' or cigarAlpha[i] == 'N' or cigarAlpha[i] == 'P')) {  }
        else if(in_region and cigarAlpha[i] == 'D')
            match_cigar += to_string(cigarLen[i]) + cigarAlpha[i];
        else if(in_region and cigarAlpha[i] == 'N'){  }
        else if(in_region and cigarAlpha[i] == 'P')
            match_cigar += to_string(cigarLen[i]) + cigarAlpha[i];
        else{
            uLONG cur_pos = last_pos + cigarLen[i];
            //cerr << last_pos << "-" << cur_pos << endl;
            if(cur_pos <= region.first)
            {}
            else if(last_pos > region.second)
            {}
            else{
                if(region.first <= last_pos and cur_pos-1 <= region.second)
                {
                    match_cigar += to_string(cigarLen[i]) + cigarAlpha[i];
                }else{
                    uLONG min_end = min(region.second, cur_pos-1);
                    uLONG max_start = max(region.first, last_pos);
                    if(min_end < max_start)
                    {
                        //cerr << "A Accodent Error" << endl;
                        //exit(-1);
                        throw runtime_error("A Accodent Error");
                    }
                    uLONG overlap = min_end - max_start + 1;
                    match_cigar += to_string(overlap) + cigarAlpha[i];
                }
            }
            last_pos = cur_pos;
        }
    }
    return match_cigar;
}


/*
enum GROUP_STATUS{
    GROUP_SUCCESS=0,    // success
    GROUP_DIFF_ID,      // map different read id
    GROUP_UNMMAPED,     // unmmaped read
    GROUP_DIFF_STRAND,  // map to different strand
    GROUP_DIFF_CHR,     // map to different chromosome
    GROUP_MULTI_MATCH,  // read map to too many regions
    GROUP_OVERLAP_GLOBAL, // two map to overlaped global region
    GROUP_ABNORMAL_MATCH_PATTERN // two reads are not MATCH_LEFT/MATCH_RIGHT pair
};
*/
GROUP_STATUS group_read_pair(const ReadPair &raw_read_pair, Sam_Record &grouped_sam_record)
{

    if(raw_read_pair.first.read_id != raw_read_pair.second.read_id)
        return GROUP_DIFF_ID;
    if(not read_is_mapped(raw_read_pair.first) or not read_is_mapped(raw_read_pair.second))
        return GROUP_UNMMAPED;
    // all map to sense or anti-sense
    if( read_is_reverse(raw_read_pair.first) ^ read_is_reverse(raw_read_pair.second) )
        return GROUP_DIFF_STRAND;

    if(raw_read_pair.first.chr_id != raw_read_pair.second.chr_id)
        return GROUP_DIFF_CHR;

    ReadPair read_pair = raw_read_pair;

    RegionArray global_matchRegion_1, global_matchRegion_2;
    get_global_match_region(read_pair.first.cigar, read_pair.first.pos, global_matchRegion_1);
    get_global_match_region(read_pair.second.cigar, read_pair.second.pos, global_matchRegion_2);
    if(global_matchRegion_1.size() != 1 or global_matchRegion_2.size() != 1)
        return GROUP_MULTI_MATCH;

    uLONG next_start = max(global_matchRegion_1.at(0).first, global_matchRegion_2.at(0).first);
    uLONG first_end = min(global_matchRegion_1.at(0).second, global_matchRegion_2.at(0).second);
    // overlap
    if( first_end >= next_start )
        return GROUP_OVERLAP_GLOBAL;
    
    uLONG gap = next_start - first_end - 1;
    if(global_matchRegion_1.at(0).first >= global_matchRegion_2.at(0).first )
        std::swap(read_pair.first, read_pair.second);

    MATCH_PATTERN first_match_pattern = cigar_match_pattern(read_pair.first.cigar);
    MATCH_PATTERN second_match_pattern = cigar_match_pattern(read_pair.second.cigar);
    if(first_match_pattern == MATCH_TOTAL or second_match_pattern == MATCH_TOTAL)
        return GROUP_ABNORMAL_MATCH_PATTERN;
    else if(first_match_pattern == MATCH_RIGHT and second_match_pattern == MATCH_LEFT)
        { /* normal */ }
    else if(first_match_pattern == MATCH_LEFT and second_match_pattern == MATCH_RIGHT)
        { /* normal */ }
    else{
        //char buffer[500];
#ifdef DEBUG_GROUP
        sprintf(buffer_1, "Unexpected Cigar: %s\t%s in %s, Skip it...\n", read_pair.first.cigar.c_str(), read_pair.second.cigar.c_str(), read_pair.first.read_id.c_str());
        cerr << buffer_1;
#endif
        return GROUP_ABNORMAL_MATCH_PATTERN;
    }

    // local-mix
    RegionArray local_matchRegion_1, local_matchRegion_2;
    get_local_match_region(read_pair.first.cigar, local_matchRegion_1);
    get_local_match_region(read_pair.second.cigar, local_matchRegion_2);
    local_matchRegion_1.at(0).second = local_matchRegion_1.back().second;
    local_matchRegion_2.at(0).second = local_matchRegion_2.back().second;

#ifdef DEBUG_GROUP
        sprintf(buffer_1, "left_match_region: %lu-%lu", local_matchRegion_1.front().first, local_matchRegion_1.front().second);
        sprintf(buffer_2, "right_match_region: %lu-%lu", local_matchRegion_2.front().first, local_matchRegion_2.front().second);
        cout << buffer_1 << endl;
        cout << buffer_2 << endl;
#endif

    uLONG local_next_start = max(local_matchRegion_1.at(0).first, local_matchRegion_2.at(0).first);
    uLONG local_first_end = min(local_matchRegion_1.at(0).second, local_matchRegion_2.at(0).second);

    string new_cigar, psuedo_cigar, new_seq, new_qaul;
    uLONG pos; uLONG gap_ambious = 0; uLONG ambiguous = 0;
    if(local_next_start > local_first_end)
    {
        uLONG left_start = local_matchRegion_1.at(0).first;
        uLONG left_end = local_matchRegion_1.at(0).second;
        uLONG right_start = local_matchRegion_2.at(0).first;
        uLONG right_end = local_matchRegion_2.at(0).second;

        // no overlap
        new_cigar += get_region_cigar(read_pair.first.cigar, Region(left_start, left_end));
        new_cigar += to_string(gap) + 'N';
        new_cigar += get_region_cigar(read_pair.second.cigar, Region(right_start, right_end));

        new_seq += read_pair.first.read_seq.substr(left_start-1, left_end-left_start+1);
        new_seq += read_pair.first.read_seq.substr(right_start-1, right_end-right_start+1);

        new_qaul += read_pair.first.read_quality.substr(left_start-1, left_end-left_start+1);
        new_qaul += read_pair.first.read_quality.substr(right_start-1, right_end-right_start+1);

        //psuedo_cigar = new_cigar;
        pos = read_pair.first.pos;
    }else{
        ambiguous = local_first_end - local_next_start + 1;
        uLONG left_start = local_matchRegion_1.at(0).first;
        uLONG left_end = local_matchRegion_1.at(0).second;
        uLONG right_start = local_matchRegion_2.at(0).first;
        uLONG right_end = local_matchRegion_2.at(0).second;

        uLONG ambiguous_left_start, ambiguous_left_end, ambiguous_right_start, ambiguous_right_end;

        if(first_match_pattern == MATCH_LEFT)
        {
            left_end -= ambiguous;
            right_start += ambiguous;
            pos = read_pair.first.pos;
            gap_ambious = ambiguous;

            ambiguous_left_start = left_end+1;
            ambiguous_left_end = left_end+ambiguous;
            ambiguous_right_start = right_start-ambiguous;
            ambiguous_right_end = right_start-1;
        }else{
            left_start += ambiguous;
            right_end -= ambiguous;
            pos = read_pair.first.pos + ambiguous;
            gap_ambious = 0;

            ambiguous_left_start = left_start-ambiguous;
            ambiguous_left_end = left_start-1;
            ambiguous_right_start = right_end+1;
            ambiguous_right_end = right_end+ambiguous;
        }

       // cout << ambiguous << "\t" << left_start << "-" << left_end << "\t" << right_start << "-" << right_end << "\t";

        if(left_end < left_start or right_end < right_start)
        {
#ifdef DEBUG_GROUP
            cerr << read_pair.first.read_id << endl;
            cerr << "\tleft_end < left_start " << left_end << " < " << left_start << endl;
            cerr << "\tright_end < right_start " << right_end << " < " << right_start << endl;
#endif
            return GROUP_INVALID_BORDER;
        }
        //assert(left_end >= left_start);
        //assert(right_end >= right_start);

        //cout << left_start << "\t" <<  left_end << endl;
        new_cigar += get_region_cigar(read_pair.first.cigar, Region(left_start, left_end));
        pos += remove_left_D(new_cigar);
        //psuedo_cigar += new_cigar;

        // comput gap length


        string left_ambiguous_cigar = get_region_cigar(read_pair.first.cigar, Region(ambiguous_left_start, ambiguous_left_end));
        string right_ambiguous_cigar = get_region_cigar(read_pair.second.cigar, Region(ambiguous_right_start, ambiguous_right_end));

#ifdef DEBUG_GROUP
        sprintf(buffer_1, "left_ambiguous_cigar(%lu-%lu): %s", ambiguous_left_start, ambiguous_left_end, left_ambiguous_cigar.c_str());
        sprintf(buffer_2, "right_ambiguous_cigar(%lu-%lu): %s", ambiguous_right_start, ambiguous_right_end, right_ambiguous_cigar.c_str());
        cout << buffer_1 << endl;
        cout << buffer_2 << endl;
#endif

        vector<int> cigarLen_1, cigarLen_2; vector<char> cigarAlpha_1, cigarAlpha_2;
        split_cigar(left_ambiguous_cigar, cigarLen_1, cigarAlpha_1);
        split_cigar(right_ambiguous_cigar, cigarLen_2, cigarAlpha_2);
        uLONG in_genome_not_in_reads = 0;
        for(size_t i=0; i<cigarAlpha_1.size(); i++)
            if(cigarAlpha_1.at(i) == 'D')
            {
                if(first_match_pattern != MATCH_LEFT)
                    pos += cigarLen_1.at(i);
                else
                    in_genome_not_in_reads += cigarLen_1.at(i);
            }
        for(size_t i=0; i<cigarAlpha_2.size(); i++)
            if(cigarAlpha_2.at(i) == 'D')
            {
                if(first_match_pattern == MATCH_LEFT)
                    in_genome_not_in_reads += cigarLen_2.at(i);
            }

        //psuedo_cigar += to_string(gap+in_genome_not_in_reads) + 'N' + to_string(ambiguous) + 'B';
        new_cigar += to_string(gap+2*gap_ambious+in_genome_not_in_reads) + 'N';
        //cout << gap << "\t" << gap_ambious << "\t" << in_genome_not_in_reads << endl;

        new_cigar += get_region_cigar(read_pair.second.cigar, Region(right_start, right_end));
        //psuedo_cigar += get_region_cigar(read_pair.second.cigar, Region(right_start, right_end));

        new_seq += read_pair.first.read_seq.substr(left_start-1, left_end-left_start+1);
        new_seq += read_pair.second.read_seq.substr(right_start-1, right_end-right_start+1);

        new_qaul += read_pair.first.read_quality.substr(left_start-1, left_end-left_start+1);
        new_qaul += read_pair.second.read_quality.substr(right_start-1, right_end-right_start+1);
    }

    grouped_sam_record = read_pair.first;
    grouped_sam_record.pos = pos; //min(read_pair.first.pos, read_pair.second.pos);
    grouped_sam_record.cigar = new_cigar;
    grouped_sam_record.flag = read_is_reverse(raw_read_pair.first) ? 16 : 0;
    grouped_sam_record.map_quanlity = 1;
    grouped_sam_record.read_seq = new_seq;
    grouped_sam_record.read_quality = new_qaul;
    grouped_sam_record.attributes.clear();
    grouped_sam_record.attributes.push_back("AB:i:"+to_string(ambiguous)); //ambigous base number
    grouped_sam_record.attributes.push_back( string("AT:Z:")+ ((first_match_pattern == MATCH_LEFT) ? "CENTER" : "FLANKING") ); //ambigous base type

    return GROUP_SUCCESS;
}

/*
enum TRIM_STATUS{
    TRIM_SUCCESS=0,    // success
    TRIM_DIFF_ID,      // map different read id
    TRIM_UNMMAPED,     // unmmaped read
    TRIM_MULTI_MATCH,  // read map to too many regions
    TRIM_OVERLAP_GLOBAL, // two map to overlaped global region
    TRIM_ABNORMAL_MATCH_PATTERN, // two reads are not MATCH_LEFT/MATCH_RIGHT pair
    TRIM_INVALID_BORDER,         // unexpected start-end border
    TRIM_OVERLAP_INCLUDE         // too much overlap
};
*/

TRIM_STATUS trim_read_pair(const ReadPair &raw_read_pair, 
                            ReadPair &trimmed_read_pair, 
                            pair<uLONG, uLONG> &bias)
{

    if(raw_read_pair.first.read_id != raw_read_pair.second.read_id)
        return TRIM_DIFF_ID;
    if(not read_is_mapped(raw_read_pair.first) or not read_is_mapped(raw_read_pair.second))
        return TRIM_UNMMAPED;

    uLONG left_should_bias = 0;
    uLONG right_should_bias = 0;

    bool read_1_reverse = read_is_reverse(raw_read_pair.first);
    bool read_2_reverse = read_is_reverse(raw_read_pair.second);
    bool same_strand = not ( read_1_reverse ^ read_2_reverse );

    bool same_chr = ( raw_read_pair.first.chr_id == raw_read_pair.second.chr_id );

    ReadPair read_pair = raw_read_pair;
    
    RegionArray global_matchRegion_1, global_matchRegion_2;
    get_global_match_region(read_pair.first.cigar, read_pair.first.pos, global_matchRegion_1);
    get_global_match_region(read_pair.second.cigar, read_pair.second.pos, global_matchRegion_2);
    if(global_matchRegion_1.size() != 1 or global_matchRegion_2.size() != 1)
        return TRIM_MULTI_MATCH;

    uLONG next_start = max(global_matchRegion_1.at(0).first, global_matchRegion_2.at(0).first);
    uLONG first_end = min(global_matchRegion_1.at(0).second, global_matchRegion_2.at(0).second);

    if( same_chr )
    {
        if( same_strand )
        {
            // overlap
            if( first_end >= next_start )
                return TRIM_OVERLAP_GLOBAL;
        }else{
            // preserve
        }
    }

    //clog << "Error 1" << endl;
    MATCH_PATTERN first_match_pattern = cigar_match_pattern(read_pair.first.cigar);
    MATCH_PATTERN second_match_pattern = cigar_match_pattern(read_pair.second.cigar);
    
    if(first_match_pattern == MATCH_TOTAL or second_match_pattern == MATCH_TOTAL)
        return TRIM_ABNORMAL_MATCH_PATTERN;

    bool has_swap = false;

    if( same_strand )
    {
        if(first_match_pattern == MATCH_RIGHT and second_match_pattern == MATCH_LEFT)
            { /* normal */ std::swap(read_pair.first, read_pair.second); has_swap = true; }
        else if(first_match_pattern == MATCH_LEFT and second_match_pattern == MATCH_RIGHT)
            { /* normal */ }
        else{
            //char buffer[500];
            return TRIM_ABNORMAL_MATCH_PATTERN;
        }
    }else{
        if(first_match_pattern == MATCH_RIGHT and second_match_pattern == MATCH_RIGHT)
            { /* normal */
                if( read_2_reverse )
                {    std::swap(read_pair.first, read_pair.second); has_swap = true; }
             }
        else if(first_match_pattern == MATCH_LEFT and second_match_pattern == MATCH_LEFT)
            { /* normal */ 
                if( read_1_reverse )
                {    std::swap(read_pair.first, read_pair.second); has_swap = true; }
            }
        else{
            return TRIM_ABNORMAL_MATCH_PATTERN;
        }
    }

    //clog << "Error 2" << endl;
    trimmed_read_pair = read_pair;
    trimmed_read_pair.first.attributes.clear();
    trimmed_read_pair.second.attributes.clear();

    // local-mix
    RegionArray local_matchRegion_1, local_matchRegion_2;
    get_local_match_region(read_pair.first.cigar, local_matchRegion_1);
    get_local_match_region(read_pair.second.cigar, local_matchRegion_2);
    local_matchRegion_1.at(0).second = local_matchRegion_1.back().second;
    local_matchRegion_2.at(0).second = local_matchRegion_2.back().second;

    uLONG ambiguous = 0;

    if( same_strand )
    {
        if( local_matchRegion_1.at(0).first <= local_matchRegion_2.at(0).first and local_matchRegion_2.at(0).second <= local_matchRegion_1.at(0).second )
            return TRIM_OVERLAP_INCLUDE;
        if( local_matchRegion_2.at(0).first <= local_matchRegion_1.at(0).first and local_matchRegion_1.at(0).second <= local_matchRegion_2.at(0).second )
            return TRIM_OVERLAP_INCLUDE;

        uLONG local_next_start = max(local_matchRegion_1.at(0).first, local_matchRegion_2.at(0).first);
        uLONG local_first_end = min(local_matchRegion_1.at(0).second, local_matchRegion_2.at(0).second);
       
        if(local_next_start > local_first_end)
        {

            uLONG start_1 = local_matchRegion_1.at(0).first;
            uLONG end_1 = local_matchRegion_1.at(0).second;
            uLONG start_2 = local_matchRegion_2.at(0).first;
            uLONG end_2 = local_matchRegion_2.at(0).second;

            trimmed_read_pair.first.cigar = get_region_cigar(read_pair.first.cigar, Region(start_1, end_1));
            trimmed_read_pair.first.cigar += to_string(end_2 - start_2 + 1) + "S";

            trimmed_read_pair.second.cigar = to_string(end_1 - start_1 + 1) + "S";
            trimmed_read_pair.second.cigar += get_region_cigar(read_pair.second.cigar, Region(start_2, end_2));

            trimmed_read_pair.first.read_seq = read_pair.first.read_seq.substr(start_1-1, end_1-start_1+1) + read_pair.first.read_seq.substr(start_2-1, end_2-start_2+1);
            trimmed_read_pair.second.read_seq = read_pair.second.read_seq.substr(start_1-1, end_1-start_1+1) + read_pair.second.read_seq.substr(start_2-1, end_2-start_2+1);

            trimmed_read_pair.first.read_quality = read_pair.first.read_quality.substr(start_1-1, end_1-start_1+1) + read_pair.first.read_quality.substr(start_2-1, end_2-start_2+1);
            trimmed_read_pair.second.read_quality = read_pair.second.read_quality.substr(start_1-1, end_1-start_1+1) + read_pair.second.read_quality.substr(start_2-1, end_2-start_2+1);

        }else{
            ambiguous = local_first_end - local_next_start + 1;
            uLONG left_start = local_matchRegion_1.at(0).first;
            uLONG left_end = local_matchRegion_1.at(0).second;
            uLONG right_start = local_matchRegion_2.at(0).first;
            uLONG right_end = local_matchRegion_2.at(0).second;

            uLONG ambiguous_left_start, ambiguous_left_end, ambiguous_right_start, ambiguous_right_end;

            left_end -= ambiguous;
            right_start += ambiguous;

            ambiguous_left_start = left_end+1;
            ambiguous_left_end = left_end+ambiguous;
            ambiguous_right_start = right_start-ambiguous;
            ambiguous_right_end = right_start-1;


            string tmp_left_cigar = get_region_cigar(read_pair.first.cigar, Region(left_start, left_end));
            auto pos_D_left_bias = remove_right_D(tmp_left_cigar);
            trimmed_read_pair.first.cigar = tmp_left_cigar;
            //trimmed_read_pair.first.cigar = get_region_cigar(read_pair.first.cigar, Region(left_start, left_end));
            trimmed_read_pair.first.cigar += to_string(right_end - right_start + 1) + "S";

            trimmed_read_pair.second.cigar = to_string(left_end - left_start + 1) + "S";
            string tmp_right_cigar = get_region_cigar(read_pair.second.cigar, Region(right_start, right_end));
            auto pos_D_right_bias = remove_left_D(tmp_right_cigar);
            trimmed_read_pair.second.cigar += tmp_right_cigar;

            trimmed_read_pair.first.read_seq = read_pair.first.read_seq.substr(left_start-1, left_end-left_start+1) + read_pair.first.read_seq.substr(right_start-1, right_end-right_start+1);
            trimmed_read_pair.second.read_seq = read_pair.second.read_seq.substr(left_start-1, left_end-left_start+1) + read_pair.second.read_seq.substr(right_start-1, right_end-right_start+1);

            trimmed_read_pair.first.read_quality = read_pair.first.read_quality.substr(left_start-1, left_end-left_start+1) + read_pair.first.read_quality.substr(right_start-1, right_end-right_start+1);
            trimmed_read_pair.second.read_quality = read_pair.second.read_quality.substr(left_start-1, left_end-left_start+1) + read_pair.second.read_quality.substr(right_start-1, right_end-right_start+1);

            right_should_bias = ambiguous + pos_D_right_bias;
            left_should_bias = ambiguous + pos_D_left_bias;

            string left_ambiguous_cigar = get_region_cigar(read_pair.first.cigar, Region(ambiguous_left_start, ambiguous_left_end));
            string right_ambiguous_cigar = get_region_cigar(read_pair.second.cigar, Region(ambiguous_right_start, ambiguous_right_end));
            vector<int> cigarLen_left, cigarLen_right; vector<char> cigarAlpha_left, cigarAlpha_right;
            split_cigar(right_ambiguous_cigar, cigarLen_left, cigarAlpha_left);
            split_cigar(right_ambiguous_cigar, cigarLen_right, cigarAlpha_right);
            
            //uLONG in_genome_not_in_reads_left = 0;

            for(size_t i=0; i<cigarAlpha_right.size(); i++)
                if(cigarAlpha_right.at(i) == 'D')
                {
                    right_should_bias += cigarLen_right.at(i);
                }

            for(size_t i=0; i<cigarAlpha_left.size(); i++)
                if(cigarAlpha_left.at(i) == 'D')
                {
                    left_should_bias += cigarLen_left.at(i);
                }

            trimmed_read_pair.second.pos = read_pair.second.pos + right_should_bias;

            bias.first = left_should_bias;
            bias.second = right_should_bias;


            if(left_end < left_start or right_end < right_start)
            {
                return TRIM_INVALID_BORDER;
            }
        }

        trimmed_read_pair.first.attributes.push_back("AB:i:"+to_string(ambiguous)); 
        trimmed_read_pair.second.attributes.push_back("AB:i:"+to_string(ambiguous)); 
        trimmed_read_pair.first.map_quanlity = 1;
        trimmed_read_pair.second.map_quanlity = 1;

        if(has_swap == true)
            std::swap( trimmed_read_pair.first, trimmed_read_pair.second );

    }else{

        //clog << "Error 3" << endl;
        if( not read_is_reverse(read_pair.first) )
        {
            swap( read_pair.first, read_pair.second );
        }

        ReadPair in_read_pair = read_pair;
        ReadPair out_read_pair;

        in_read_pair.first.read_seq = reverse_comp( in_read_pair.first.read_seq );
        std::reverse( in_read_pair.first.read_quality.begin(), in_read_pair.first.read_quality.end() );
        in_read_pair.first.cigar = reverse_cigar(in_read_pair.first.cigar);
        in_read_pair.first.flag &= ~16UL;
        in_read_pair.first.pos = in_read_pair.second.pos + 1000; // a random number

        pair<uLONG, uLONG> bias(0, 0);

        TRIM_STATUS trim_result = trim_read_pair(in_read_pair, out_read_pair, bias);
        if(trim_result != TRIM_SUCCESS)
            return trim_result;

        out_read_pair.first.read_seq = reverse_comp( out_read_pair.first.read_seq );
        std::reverse( out_read_pair.first.read_quality.begin(), out_read_pair.first.read_quality.end() );
        out_read_pair.first.cigar = reverse_cigar(out_read_pair.first.cigar);
        out_read_pair.first.flag |= 16UL;

        if(first_match_pattern == MATCH_LEFT)
        {
            out_read_pair.first.pos = read_pair.first.pos;
            out_read_pair.second.pos = read_pair.second.pos;
        }else{
            out_read_pair.first.pos = read_pair.first.pos + bias.first;
            out_read_pair.second.pos = read_pair.second.pos + bias.second;
        }

        trimmed_read_pair = out_read_pair;
    }

    return TRIM_SUCCESS;
}






/*
void print_match_region(const RegionArray &matchRegion)
{
    for(auto iter=matchRegion.cbegin(); iter!=matchRegion.cend(); iter++)
        cout << iter->first << "-" << iter->second << "\t";
}
*/

void filter_unmapped_record(vector<Sam_Record> &read_records)
{
    auto end = std::remove_if(read_records.begin(), 
        read_records.end(), 
        [](const Sam_Record &record){ return (not read_is_mapped(record)); });
    read_records.erase(end, read_records.end());
}

void filter_ungapped_record(vector<Sam_Record> &read_records)
{
    auto end = std::remove_if(read_records.begin(), 
        read_records.end(), 
        [](const Sam_Record &record){ return (not read_is_gapped(record)); });
    read_records.erase(end, read_records.end());
}

void filter_gapped_record(vector<Sam_Record> &read_records)
{
    auto end = std::remove_if(read_records.begin(), 
        read_records.end(), 
        [](const Sam_Record &record){ return read_is_gapped(record); });
    read_records.erase(end, read_records.end());
}

void filter_reversed_record(vector<Sam_Record> &read_records)
{
    auto end = std::remove_if(read_records.begin(), 
        read_records.end(), 
        [](const Sam_Record &record){ return read_is_reverse(record); });
    read_records.erase(end, read_records.end());
}

void filter_secondary_record(vector<Sam_Record> &read_records)
{
    auto end = std::remove_if(read_records.begin(), 
        read_records.end(), 
        [](const Sam_Record &record){ return (not read_is_primary(record)); });
    read_records.erase(end, read_records.end());
}



/*
void filter_second_in_pair(vector<Sam_Record> &read_records)
{
    auto end = std::remove_if(read_records.begin(), 
        read_records.end(), 
        [](const Sam_Record &record){ return (not read_is_mate_reversed(record) and read_is_first_in_pair(record)); });
    read_records.erase(end, read_records.end());
}
*/

void statistic_sam_abundance(const string &sam_file, 
        MapStringT<uLONG> &uniq_mapped_reads, 
        MapStringT<double> &multi_mapped_reads,
        MapStringT<uLONG> &ref_len,
        const RPKM_PARAM &param)
{
    uniq_mapped_reads.clear(); multi_mapped_reads.clear();

    ifstream IN(sam_file, ifstream::in);
    if(not IN)
    {
        throw runtime_error( "Bad_Input_File: "+sam_file );
    }
    
    Sam_Head sam_head;
    read_sam_head(IN, sam_head);

    ref_len = sam_head.trans_len;
    for_each(ref_len.cbegin(), ref_len.cend(), [&](const pair<string, uLONG> &cur_len){ uniq_mapped_reads[cur_len.first]=0; multi_mapped_reads[cur_len.first]=0; });

    uLONG cur_read_nums = 0;
    //Sam_Record read_record;
    vector<Sam_Record> read_records;
    while(not IN.eof())
    {
        bool success = read_a_read_record(IN, read_records);
        //bool success = read_a_sam_record(IN, read_record);
        if(not success)
            break;

        cur_read_nums++;
        if(param.p_verbose_out and cur_read_nums % 100000 == 0)
            *param.p_verbose_out << "\tLog -- Currently Read " << cur_read_nums << " Reads..." << endl;

        filter_unmapped_record(read_records);

        // if not preserve reversed read then remove it
        if(not param.preserve_reverse_map)
            filter_reversed_record(read_records);

        if(read_records.size() == 0)
            continue;

        if(read_is_multimap(read_records) and not param.preserve_multimap)
            continue;

        if(read_records.size() == 1)
            // unique map: single-end
        {
            ++uniq_mapped_reads.at(read_records[0].chr_id);
        }else
        {
            if(read_is_paired(read_records[0]))
            {
                if(read_is_multimap(read_records))
                {
                    // map map: pair-end
                    double ave_rpkm = 0.5/read_records.size();
                    for_each(read_records.cbegin(), read_records.cend(), [&](const Sam_Record &record){ multi_mapped_reads.at(record.chr_id) += ave_rpkm; });
                }else{
                    // unique map: pair-end
                    if(read_records[0].chr_id == read_records[1].chr_id)
                        ++uniq_mapped_reads.at(read_records[0].chr_id);
                    else
                    {
                        ++uniq_mapped_reads.at(read_records[0].chr_id);
                        ++uniq_mapped_reads.at(read_records[1].chr_id);
                    }
                }
            }else{
                // multi-map pair-end
                double ave_rpkm = 1.0/read_records.size();
                for_each(read_records.cbegin(), read_records.cend(), [&](const Sam_Record &record){ multi_mapped_reads.at(record.chr_id) += ave_rpkm; });
            }
        }
    }
    IN.close();
}

void calc_rpkm(const string &sam_file, MapStringT<double> &rpkm, const RPKM_PARAM &param)
{
    rpkm.clear();

    MapStringT<uLONG> uniq_mapped_reads; MapStringT<double> multi_mapped_reads; MapStringT<uLONG> ref_len;
    statistic_sam_abundance(sam_file, uniq_mapped_reads, multi_mapped_reads, ref_len, param);

    double total_mapped_reads = 0;
    for_each(ref_len.cbegin(), ref_len.cend(),[&](const pair<string, uLONG> &cur_len){ total_mapped_reads += uniq_mapped_reads[cur_len.first] + multi_mapped_reads[cur_len.first]; });
    cout << "total_mapped_reads: " << total_mapped_reads << endl;

    for(auto iter=ref_len.cbegin(); iter!=ref_len.cend(); iter++)
    {
        rpkm[iter->first] = rpkm_func( total_mapped_reads, multi_mapped_reads[iter->first] + uniq_mapped_reads[iter->first], iter->second );
    }
}




};

