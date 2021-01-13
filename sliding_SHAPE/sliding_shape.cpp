
#include "sliding_shape.h"
#include <math.h>

Color::Modifier RED(Color::FG_RED);
Color::Modifier DEF(Color::FG_DEFAULT);
Color::Modifier YELLOW(Color::FG_YELLOW);


Map_Record::Map_Record(const StringArray &line_data)
{
    if(line_data[1] == "+")
        strand = POSITIVE;
    else
        strand = NEGATIVE;

    for(auto it=line_data.cbegin()+2; it!=line_data.cend(); it+=2)
        regions.emplace_back( stoul(*it), stoul(*(it+1)) );

    if(regions.size() == 0)
    {
        string str;
        for(string v: line_data)
            str += v;
        throw Unexpected_Error("FATAL Error: invliad line: "+str);
    }
}

/**** Junctions ****/

// read sjdbList.fromGTF.out.tab file (in STAR index directory)
void load_junctions(const string &file_name, MapStringT<JunctionArray> &junctions, const MapStringuLONG &chr_size)
{
    junctions.clear();
    ifstream IN(file_name, ifstream::in);
    if(not IN)
    {
        cerr << RED << "FATAL Error: cannot read " << file_name << DEF << endl;
        exit(-1);
    }

    string line;
    while(getline(IN, line))
    {
        StringArray data;
        split(line, data);

        const string chr_id = data[0]+data[3];
        uLONG s = stoul(data[1]);
        uLONG e = stoul(data[2]);

        if( chr_size.find(chr_id) != chr_size.end() )
        {
            if( chr_size.at(chr_id)>e )
                junctions[ chr_id ].emplace_back(s, e);
            else
                cerr << RED << "Warning: " << line << " in juntion file exceed the length of chromosome " << chr_size.at(chr_id) << " , skip it" << DEF << endl;
        }else{
            cerr << RED << "Warning: " << chr_id << " in juntion file not in size file, skip it" << DEF << endl;
        }
    }
    IN.close();

    for(auto it=junctions.begin(); it!=junctions.end(); it++)
        sort(it->second.begin(), it->second.end());
}



// build junction index, left means build with junction start; right means build with junction end
void buildLeftJunctionMap(JunctionArray &junctions, map<uLONG, vector<Junction*>> &junc_map, const uLONG &binsize)
{
    junc_map.clear();
    for(auto it=junctions.begin(); it!=junctions.end(); it++)
    {
        uLONG index = it->first / binsize;
        junc_map[index].push_back( &(*it) );
    }
}

void buildRightJunctionMap(JunctionArray &junctions, map<uLONG, vector<Junction*>> &junc_map, const uLONG &binsize)
{
    junc_map.clear();
    for(auto it=junctions.begin(); it!=junctions.end(); it++)
    {
        uLONG index = it->second / binsize;
        junc_map[index].push_back( &(*it) );
    }
}

// count junction reads
void build_junction_support(const vector<Map_Record> &record_array, JunctionArray &junctions)
{
    map<uLONG, vector<Junction*>> junc_map;
    const uLONG binsize = 100000;
    buildLeftJunctionMap(junctions, junc_map, binsize);

    for(auto it=record_array.cbegin(); it!=record_array.cend(); it++)
    {
        for(auto it2=it->regions.cbegin(); it2!=it->regions.cend()-1; it2++)
        {
            Region cur_junc(it2->second+1, (it2+1)->first-1);
            uLONG index = cur_junc.first / binsize;
            if(junc_map.find(index) == junc_map.end())
                continue;
            for(Junction* const p: junc_map.at(index))
            {
                if(p->first == cur_junc.first and p->second == cur_junc.second)
                {
                    p->support++;
                    break;
                }
                // TODO: if it is a unannotated junction...
            }
        }
    }
}


// preserve junctions with max supported reads when junctions overlap
void combine_junction(JunctionArray &junctions)
{
    JunctionArray tmp_junc;

    tmp_junc.push_back(junctions[0]);
    for(uLONG i=1; i<junctions.size(); i++)
    {
        uLONG s = tmp_junc.back().first;
        uLONG e = tmp_junc.back().second;
        if(s <= junctions[i].second and junctions[i].first <= e)
        {
            if(junctions[i].support > tmp_junc.back().support)
                tmp_junc.back() = junctions[i];
        }else{
            tmp_junc.push_back(junctions[i]);
        }
    }
    junctions = tmp_junc;

    sort(junctions.begin(), junctions.end());
}


// check junction overlap and combine junction when overklap
void check_overlap(JunctionArray &junctions, const uLONG &size)
{
    for(uLONG i=0; i<junctions.size()-1; )
    {
        uLONG cs = junctions[i].first;
        uLONG ce = junctions[i].second;
        uLONG ns = junctions[i+1].first;
        uLONG ne = junctions[i+1].second;

        if(ne > size+1)
        {
            cerr << RED << "FATAL Error: junction exceed chromsome length(" << size << "): " << ns << "-" << ne << "; remove it" << DEF << endl;
            junctions.erase( junctions.begin()+i+1 );
            continue;
        }

        if(cs<=ne+1 and ns<=ce+1)
        {
            junctions[i].first = cs>ns ? ns : cs;
            junctions[i].second = ce>ne ? ne : ce;
            if(junctions[i+1].support > junctions[i].support)
                junctions[i].support = junctions[i+1].support;
            
            junctions.erase(junctions.begin()+i+1);
            cerr << RED << "Warning: combine overlap: " << cs << "-" << ce << ";" << ns << "-" << ne << DEF << endl;
        }

        ++i;
    }
}

// write junctions into a file
void write_junctions(const MapStringT<JunctionArray> &junctions, const string &outFile)
{
    ofstream OUT(outFile, ofstream::out);
    if(not OUT)
    {
        cerr << RED << "FATAL Error: cannot write " << outFile << DEF << endl;
        exit(-1);
    }

    for(auto it=junctions.cbegin(); it!=junctions.cend(); it++)
    {
        const string &chr_id = it->first;
        string chr = chr_id.substr(0, chr_id.size()-1);
        char strand = chr_id[chr_id.size()-1];
        for( const Junction &ja: it->second)
            OUT << chr << "\t" << strand << "\t" << ja.first << "\t" << ja.second << "\t" << ja.support << "\n";
    }

    OUT.close();
}


/**** Read Tab file (from sam file) ****/

// read a single chromosome data, input file must be sorted!!
bool read_chr(vector<Map_Record> &record_array, ifstream *hander, string &chr_id)
{
    if(hander->eof())
        return false;

    record_array.clear();

    string line; StringArray data;
    getline(*hander, line);
    split(line, '\t', data);
    chr_id = data[0]+data[1];
    record_array.emplace_back(data);

    auto last_pos = hander->tellg();
    while(getline(*hander, line))
    {
        split(line, '\t', data);
        if(chr_id != data[0]+data[1])
        {
            hander->seekg(last_pos);
            break;
        }
        record_array.emplace_back(data);
        last_pos = hander->tellg();
    }

    return true;
}


/**** Read chromosome size file ****/

// load chromosome size file, it can be chrNameLength.txt from STAR index directory
void load_chr_size(const string &file_name, MapStringuLONG &chr_size)
{
    ifstream IN(file_name, ifstream::in);
    if(not IN)
    {
        cerr << RED << "FATAL Error: cannot read " << file_name << DEF << endl;
        exit(-1);
    }

    chr_size.clear();

    string line;
    while(getline(IN, line))
    {
        StringArray data;
        split(line, data);
        chr_size[ data[0] ] = stoul(data[1]);
        chr_size[ data[0]+"+" ] = stoul(data[1]);
        chr_size[ data[0]+"-" ] = stoul(data[1]);
    }
    IN.close();
}



/**** calculate RT and BD ****/

// calculate chromsome RT and BD in positive strand (input record_array must be in positive strand)
void calc_chr_BDRT_Pos(uIntArray &BD, uIntArray &RT, 
        const vector<Map_Record> &record_array, 
        JunctionArray &junctions,
        const uLONG &size,
        const uLONG &BD_ext,
        const uLONG &binsize)
{
    // build junction map
    map<uLONG, vector<Junction* >> junc_map;
    buildRightJunctionMap(junctions, junc_map, binsize);

    // scan
    for(auto it=record_array.crbegin(); it!=record_array.crend(); it++)
    {
        // add BD
        for(uLONG r_i=0; r_i<it->regions.size(); r_i++)
            for(uLONG i=it->regions[r_i].first; i<=it->regions[r_i].second; i++)
                if(i<=size)
                    ++BD[i];

        uLONG s = it->regions[0].first;
        uLONG bin_id = s/binsize;

        // add RT
        if(junc_map.find(bin_id) == junc_map.end())
        {
            if(s-1 <= size)
                ++RT[s-1];
        }else{
            bool find = false;
            for(const Junction* const jp: junc_map[bin_id])
            {
                if(jp->second+1 == s)
                {
                    find = true;
                    if(jp->first-1 <= size)
                        ++RT[jp->first-1];
                    break;
                }
            }
            if(not find)
                if(s-1 <= size)
                    ++RT[s-1];
        }

        // ext BD
        if(BD_ext > 0)
        {
            uLONG ext = BD_ext;
            if(BD_ext > s)
                ext = s;
            if(junc_map.find(bin_id) == junc_map.end())
            {
                for(uLONG i=1; i<=ext; i++)
                    if(s-i <= size)
                        ++BD[s-i];
            }else{
                const vector<Junction *> &cur_junc_map = junc_map[bin_id];
                
                bool find = false;
                for(uINT i=0; i<cur_junc_map.size(); i++)
                {
                    const Junction* const cur_jp = cur_junc_map[i];
                    if(cur_jp->first <= s and s <= cur_jp->second)
                    {
                        // in junction
                        uLONG ext_start;
                        if(s-cur_jp->first > ext)
                            ext_start = s - ext;
                        else
                            ext_start = cur_jp->first;
                        
                        for(uLONG i=ext_start; i<s; i++)
                            if(i <= size)
                                ++BD[i];

                        find = true;
                        break;
                    }else
                    {
                        Junction* nex_jp = nullptr;
                        if(i != cur_junc_map.size()-1)
                            nex_jp = cur_junc_map[i+1];
                        else
                            nex_jp = new Junction( (bin_id+1)*binsize, (bin_id+1)*binsize );

                        if(cur_jp->second < s and s < nex_jp->first)
                        {
                            // between junction
                            uLONG ext_left = ext;
                            uLONG start = cur_jp->second;
                            uLONG cur_pos = s - 1;
                            long cur_jun_index = i;
                            while(ext_left > 0)
                            {
                                if(start==cur_pos)
                                {
                                    if(cur_jun_index >= 1)
                                    {
                                        --cur_jun_index;
                                        start = cur_junc_map[cur_jun_index]->second;
                                        cur_pos = cur_junc_map[cur_jun_index+1]->first - 1;
                                    }
                                    else{
                                        cur_pos = cur_junc_map[0]->first - 1;
                                        start = 0;
                                    }
                                }
                                if(cur_pos <= size)
                                    ++BD[cur_pos];
                                --ext_left;
                                --cur_pos;
                            }
                            find = true;
                           // break;
                        }
                        if(i == cur_junc_map.size()-1)
                            delete nex_jp;
                        if(find)
                            break;
                    }
                }
                if(not find)
                    // before the first junction
                    for(uLONG i=1; i<=ext; i++)
                        if(s-i <= size)
                            ++BD[s-i];
            }
        }
    }
}

// calculate chromsome RT and BD in negative strand (input record_array must be in negative strand)
void calc_chr_BDRT_Neg(uIntArray &BD, uIntArray &RT, 
        const vector<Map_Record> &record_array, 
        JunctionArray &junctions,
        const uLONG &size,
        const uLONG &BD_ext,
        const uLONG &binsize)
{
    // build junction map
    map<uLONG, vector<Junction* >> junc_map;
    buildLeftJunctionMap(junctions, junc_map, binsize);

    //uLONG index = 0;

    // scan
    for(auto it=record_array.cbegin(); it!=record_array.cend(); it++)
    {
        //clog << ++index << it->regions.back().second << endl;

        // add BD
        for(uLONG r_i=0; r_i<it->regions.size(); r_i++)
            for(uLONG i=it->regions[r_i].first; i<=it->regions[r_i].second; i++)
                if(i <= size)
                    ++BD[i];

        uLONG s = it->regions.back().second;
        uLONG bin_id = s/binsize;

        // add RT
        if(junc_map.find(bin_id) == junc_map.end())
        {
            if(s+1 <= size)
                ++RT[s+1];
        }else{
            bool find = false;
            for(const Junction* const jp: junc_map[bin_id])
            {
                if(jp->first-1 == s)
                {
                    find = true;
                    if(jp->second+1 <= size)
                        ++RT[jp->second+1];
                    break;
                }
            }
            if(not find)
                if(s+1 <= size)
                    ++RT[s+1];
        }

        // ext BD
        if(BD_ext > 0)
        {
            uLONG ext = BD_ext;
            if(BD_ext+s > size)
                ext = size-s;
            if(junc_map.find(bin_id) == junc_map.end())
            {
                for(uLONG i=1; i<=ext; i++)
                    if(s+i<=size)
                        ++BD[s+i];
            }else{
                const vector<Junction *> &cur_junc_map = junc_map[bin_id];
                
                bool find = false;
                for(uLONG i=0; i<cur_junc_map.size(); i++)
                {
                    const Junction* const cur_jp = cur_junc_map[i];
                    if(cur_jp->first <= s and s <= cur_jp->second)
                    {
                        // in junction
                        uLONG ext_end;
                        if(cur_jp->second-s > ext)
                            ext_end = s+ext;
                        else
                            ext_end = cur_jp->second;
                        
                        for(uLONG i=s+1; i<=ext_end; i++)
                            if(i<=size)
                                ++BD[i];

                        find = true;
                        break;
                    }else
                    {
                        Junction* last_jp = nullptr;
                        if(i != 0)
                            last_jp = cur_junc_map[i-1];
                        else
                            last_jp = new Junction( bin_id*binsize, bin_id*binsize );

                        if(last_jp->second < s and s < cur_jp->first)
                        {
                            // between junction
                            uLONG ext_left = ext;
                            uLONG end = cur_jp->first;
                            uLONG cur_pos = s + 1;
                            uLONG cur_jun_index = i;
                            while(ext_left > 0)
                            {
                                if(end==cur_pos)
                                {
                                    if(cur_jun_index < cur_junc_map.size()-1)
                                    {
                                        ++cur_jun_index;
                                        end = cur_junc_map[cur_jun_index]->first;
                                        cur_pos = cur_junc_map[cur_jun_index-1]->second + 1; // cur_jp->first-1;
                                    }
                                    else{
                                        cur_pos = cur_junc_map.back()->second + 1;
                                        end = size; // max value
                                    }
                                }
                                if(cur_pos<=size)
                                    ++BD[cur_pos];
                                --ext_left;
                                ++cur_pos;
                            }
                            find = true;
                        }
                        if(i == 0)
                            delete last_jp;
                        if(find)
                            break;
                    }
                }
                if(not find)
                    // before the first junction
                    for(uLONG i=1; i<=ext; i++)
                        if(s+i<=size)
                            ++BD[s+i];
            }
        }
    }

    //clog << "Finish...FIsh" << endl;
}

/**** Cutoff functions ****/

// calculate averaged score from a float array, if the valid score number is less than half of number, it will return null
float mean_score(const FloatArray &data_array, unsigned short &valid_count)
{
    float sum = 0;
    valid_count = 0;

    if(data_array.size() < 2)
        return null;

    for(const float &v: data_array)
        if(v != null)
        {
            sum += v;
            ++valid_count;
        }
    if(valid_count > data_array.size()/2)
        return sum/valid_count;
    else
        return null;
}

// discriminate if a region is low expressed, is the number of base exceed minBD exceed bd.size()/5, then it is not a low exp region
bool low_exp_region(const deque<float> &bd, const uLONG &minBD)
{
    uLONG nums = 0;
    for(const float& d: bd)
        if(d >= minBD)
            ++nums;

    if(nums >= bd.size()/5)
        return false;
    else
        return true;
}


// if bd/rt is below a cutoff then skip a calculate
bool valid_cov(const deque<float> &bd, const deque<float> &rt, const float &min_bd, const float &min_rt)
{
    uLONG big_bd_num = count_if(bd.cbegin(), bd.cend(), [&min_bd](const float &v){ return v>=min_bd; });
    float ave_rt = std::accumulate(rt.cbegin(), rt.cend(), 0.0) / rt.size();

    if(big_bd_num > bd.size()/2 and ave_rt >= min_rt)
        return true;
    else
        return false;
}


//###############################
//###############################
//    Insert 1 /**** icSHAPE fansion (with NAI+DMSO) ****/
//###############################
//###############################

/**** icSHAPE fansion (with NAI+DMSO) ****/

// calculate SHAPE score with NAI and DMSO RT/BD in non-junction regions
void sliding_non_junction(  const uIntArray &NAI_RT, const uIntArray &DMSO_RT,
                            const uIntArray &NAI_BD, const uIntArray &DMSO_BD,
                            const JunctionArray &junctions, FloatArray &score, 
                            const string &chr_id, const STRAND &strand, 
                            ofstream &OUT, const icSHAPE_Param &param,
                            const vector<bool> &chr_mask, const string &chr_seq)
{
    bool use_mask = chr_mask.empty() ? false : true;
    if(use_mask and NAI_RT.size() != chr_mask.size())
        throw runtime_error("NAI_RT.size() != chr_mask.size() -- "+to_string(NAI_RT.size())+" != "+to_string(chr_mask.size()));

    const char strand_char = (strand == POSITIVE) ? '+' : '-';
    const uINT out_min_cov = param.out_min_cov;
    const uINT min_cov = param.min_cov;
    const uINT wsize = param.wsize;
    const uINT org_wstep = param.wstep;
    uINT wstep = org_wstep;

    deque<bool> mask;
    deque<float> nai_rt;
    deque<float> nai_bd;
    deque<float> dmso_rt;
    deque<float> dmso_bd;

    const uLONG chr_size = NAI_RT.size() - 1;
    if(DMSO_RT.size() != chr_size+1 or NAI_BD.size() != chr_size+1 or DMSO_BD.size() != chr_size+1)
    {
        cerr << RED << "FATAL Error: RT,BD different length" << DEF << endl;
        exit(-1);
    }

    //////////////////////////////////////////Only useful when non-junction///////////////////////////////////////////////////
    if(param.no_sliding)
    {

        if(chr_size < 50)
        {
            cerr << RED << "Warning: Chromosome " << chr_id << "(" << strand_char << ")" << " is too short. Skip it" << chr_size;
            return;
        }

        FloatArray scores;
        for(uLONG i=0; i<chr_size-30; i++)
        {
            if(use_mask)
                mask.push_back(chr_mask[i]);
            nai_rt.push_back(NAI_RT[i]);
            nai_bd.push_back(NAI_BD[i]);
            dmso_rt.push_back(DMSO_RT[i]);
            dmso_bd.push_back(DMSO_BD[i]);
        }

        if(use_mask)
            calculate_mask_score(nai_rt, nai_bd, dmso_rt, dmso_bd, scores, mask, param);
        else
            calculate_score(nai_rt, nai_bd, dmso_rt, dmso_bd, scores, param);

        //for(uLONG i=6; i<chr_size-35; i++)
        //    score[i] = scores[i];

        for(uLONG i=6; i<chr_size-30; i++)
            if(DMSO_BD[i] < min_cov)
                score[i] = null;
            else
                score[i] = scores[i];

        for(uLONG index=1; index<chr_size; index++)
        {

            string base;
            bool masked = true;
            if(use_mask)
            {
                base = base+chr_seq[index-1]+"\t";
                masked = chr_mask[index];
            }

            if(masked and (NAI_RT[index] >= 1 or DMSO_RT[index] >= 1 or DMSO_BD[index] >= out_min_cov))
            {
                OUT << chr_id << '\t' << strand_char << '\t' << index << '\t' << base
                    << NAI_RT[index] << '\t' << NAI_BD[index] << '\t'
                    << DMSO_RT[index] << '\t' << DMSO_BD[index] << '\t'
                    << score[index] << '\t' << (score[index]==null ? 0 : 1) << '\t' 
                    << score[index] << '\n';
            }
        }
        return;
    }

    if(chr_size < wsize*2)
    {
        cerr << RED << "Warning: Chromosome " << chr_id << "(" << strand_char << ")" << " is too short. Skip it" << chr_size;
        cerr << "Not use sliding window (-non-sliding) or use a shorter window size (-wsize) is recommend for this one" << DEF << endl;
        return;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    uINT ji = 0;

    deque<uLONG> index_array;

    bool has_junction = (junctions.size() >= 1);
    uLONG juncNum = junctions.size();

    uLONG cur_i = 1;
    // move to start
    while(NAI_RT[cur_i] == 0 and cur_i<=chr_size)
        ++cur_i;
    if(has_junction and ji<juncNum)
        while(junctions[ji].first < cur_i)
            ++ji;

    uLONG de_len = 0;
    while(de_len < wsize and cur_i<=chr_size)
    {
        if(has_junction and ji<juncNum)
            if(junctions[ji].first == cur_i)
                cur_i = junctions[ji++].second+1;

        nai_rt.push_back(NAI_RT[cur_i]);
        nai_bd.push_back(NAI_BD[cur_i]);
        dmso_rt.push_back(DMSO_RT[cur_i]);
        dmso_bd.push_back(DMSO_BD[cur_i]);
        if(use_mask)
            mask.push_back(chr_mask[cur_i]);

        index_array.push_back(cur_i);

        ++cur_i;
        ++de_len;
    }

    //clog << "\tFirst calculate" << endl;
    deque<FloatArray> precalculated;
    FloatArray scores;

    if(use_mask)
    {
        calculate_mask_score(nai_rt, nai_bd, dmso_rt, dmso_bd, scores, mask, param);
    }else{
        calculate_score(nai_rt, nai_bd, dmso_rt, dmso_bd, scores, param);
    }

    //calculate_score(nai_rt, nai_bd, dmso_rt, dmso_bd, scores, param);
    for(const float &score: scores)
    {
        FloatArray fa;
        fa.push_back(score);
        precalculated.push_back(fa);
    }
    //clog << "\tprecalculated size:" <<  << endl;

    //clog << "\twindow step" << endl;
    while(cur_i<=chr_size)
    {
        for(long i=0; i<wstep and cur_i<=chr_size; i++)
        {
            if(cur_i % 1000000 == 0)
                clog << "\t lines: " << cur_i << "(" << 100.0*cur_i/chr_size << "%) wStep: " << wstep << endl;

            uLONG index = index_array.front();
            unsigned short valid_count;
            double shape_score = mean_score(precalculated.front(), valid_count);
            //c_count[ index ] = valid_count;
            if(DMSO_BD[index] < min_cov)
                shape_score = null;
            score[ index ] = shape_score;

            string base;
            bool masked = true;
            if(use_mask)
            {
                base = base+chr_seq[index-1]+"\t";
                masked = chr_mask[index];
            }

            if(masked and (NAI_RT[index] >= 1 or DMSO_RT[index] >= 1 or DMSO_BD[index] >= out_min_cov))
            {
                OUT << chr_id << '\t' << strand_char << '\t' << index << '\t' << base
                    << NAI_RT[index] << '\t' << NAI_BD[index] << '\t'
                    << DMSO_RT[index] << '\t' << DMSO_BD[index] << '\t'
                    << shape_score << '\t' << valid_count << '\t';

                for(float v: precalculated.front())
                    OUT << v << ',';

                OUT << '\n';
            }

            index_array.pop_front();
            precalculated.pop_front();

            nai_rt.pop_front();
            nai_bd.pop_front();
            dmso_rt.pop_front();
            dmso_bd.pop_front();
            if(use_mask)
                mask.pop_front();

            if(has_junction and ji<juncNum)
                if(junctions[ji].first == cur_i)
                    cur_i = junctions[ji++].second+1;

            nai_rt.push_back(NAI_RT[cur_i]);
            nai_bd.push_back(NAI_BD[cur_i]);
            dmso_rt.push_back(DMSO_RT[cur_i]);
            dmso_bd.push_back(DMSO_BD[cur_i]);
            if(use_mask)
                mask.push_back(chr_mask[cur_i]);

            index_array.push_back(cur_i);

            ++cur_i;
        }

        uLONG minBD=10;
        if(low_exp_region(dmso_bd, minBD))
        {
            wstep = wsize - 1;
            scores.assign(nai_rt.size(), null);
        }
        else
        {
            wstep = org_wstep;

            if(use_mask)
            {
                calculate_mask_score(nai_rt, nai_bd, dmso_rt, dmso_bd, scores, mask, param);
            }else{
                calculate_score(nai_rt, nai_bd, dmso_rt, dmso_bd, scores, param);
            }
            //calculate_score(nai_rt, nai_bd, dmso_rt, dmso_bd, scores, param);
        }

        for(uLONG i=0; i<scores.size(); i++)
        {
            if(i < precalculated.size())
                precalculated[i].push_back(scores[i]);
            else{
                FloatArray fa;
                fa.push_back(scores[i]);
                precalculated.push_back(fa);
            }
        }
    }

    if(has_junction)
        if(ji <= junctions.size()-1)
            cerr << RED << "Error: Junction index not reach final " << junctions.at(ji).first << "-" << junctions.at(ji).second << DEF << endl;

    while(precalculated.size() > 0)
    {
        uLONG index = index_array.front();

        unsigned short valid_count;
        float shape_score = mean_score(precalculated.front(), valid_count);
        if(DMSO_BD[index] < min_cov)
            shape_score = null;
        score[ index ] = shape_score;

        string base;
        bool masked = true;
        if(use_mask)
        {
            base = base+chr_seq[index-1]+"\t";
            masked = chr_mask[index];
        }

        if(masked and (NAI_RT[index] >= 1 or DMSO_RT[index] >= 1 or DMSO_BD[index] >= out_min_cov))
        {
            OUT << chr_id << '\t' << strand_char << '\t' << index << '\t' << base
                << NAI_RT[index] << '\t' << NAI_BD[index] << '\t'
                << DMSO_RT[index] << '\t' << DMSO_BD[index] << '\t'
                << shape_score << '\t' << valid_count << '\t';

            for(float v: precalculated.front())
                OUT << v << ',';

            OUT << '\n';
        }

        index_array.pop_front();
        precalculated.pop_front();
    }

    for(uLONG i=0; i<30 and i<=chr_size; i++)
        score[i] = null;
    uLONG s = 0;
    if(chr_size > 30)
        s = chr_size - 30;
    for(uLONG i=s; i<=chr_size; i++)
        score[i] = null;
}

// sliding all junctions in junction regions with NAI + DMSO RT/BD
void sliding_junction( const uIntArray &NAI_RT, const uIntArray &DMSO_RT,
                        const uIntArray &NAI_BD, const uIntArray &DMSO_BD,
                        const JunctionArray &junctions, FloatArray &score, 
                        const string &chr_id, const STRAND &strand, 
                        ofstream &OUT, const icSHAPE_Param &param,
                        const vector<bool> &chr_mask, const string &chr_seq)
{

    ///////////////////////////
    if(param.no_sliding)
        return;
    ///////////////////////////

    const uLONG chr_size = NAI_RT.size() - 1;
    if(DMSO_RT.size() != chr_size+1 or NAI_BD.size() != chr_size+1 or DMSO_BD.size() != chr_size+1)
    {
        cerr << RED << "FATAL Error: RT,BD different length" << DEF << endl;
        exit(-1);
    }

    uLONG junc_num = junctions.size();

    uLONG index = 0;
    for(const Junction &junc: junctions)
    {
        if(++index % 100 == 0)
            clog << "\t junc ratio: " << 100.0*index/junc_num << "%" << endl;
        sliding_single_junction(NAI_RT, DMSO_RT, NAI_BD, DMSO_BD, junc.first, junc.second, score, chr_id, strand, OUT, param, chr_mask, chr_seq);
    }
}


// sliding a single junction with NAI + DMSO RT/BD

void sliding_single_junction(const uIntArray &NAI_RT, const uIntArray &DMSO_RT,
                            const uIntArray &NAI_BD, const uIntArray &DMSO_BD,
                            const uLONG &start, const uLONG &end, FloatArray &score,
                            const string &chr_id, const STRAND &strand, 
                            ofstream &OUT, const icSHAPE_Param &param,
                            const vector<bool> &chr_mask, const string &chr_seq)
{
    bool use_mask = chr_mask.empty() ? false : true;
    if(use_mask and NAI_RT.size() != chr_mask.size())
        throw runtime_error("NAI_RT.size() != chr_mask.size() -- "+to_string(NAI_RT.size())+" != "+to_string(chr_mask.size()));

    const char strand_char = (strand == POSITIVE) ? '+' : '-';
    const uINT out_min_cov = param.out_min_cov;
    const uINT min_cov = param.min_cov;
    const uINT wsize = param.wsize;
    const uINT org_wstep = param.wstep;
    uINT wstep = org_wstep;

    const uLONG junction_len = end-start+1;
    if(junction_len < wsize)
        return;

    deque<bool> mask;
    deque<float> nai_rt;
    deque<float> nai_bd;
    deque<float> dmso_rt;
    deque<float> dmso_bd;

    deque<uLONG> index_array;

    uLONG cur_i = 1;
    if(start > wsize)
        cur_i = start - wsize + 1;

    uLONG de_len = 0;
    while(de_len < wsize and de_len<junction_len and cur_i<=end)
    {
        nai_rt.push_back(NAI_RT[cur_i]);
        nai_bd.push_back(NAI_BD[cur_i]);
        dmso_rt.push_back(DMSO_RT[cur_i]);
        dmso_bd.push_back(DMSO_BD[cur_i]);
        if(use_mask)
            mask.push_back(chr_mask[cur_i]);

        index_array.push_back(cur_i);

        ++cur_i;
        ++de_len;
    }

    //clog << "\tFirst calculate" << endl;
    deque<FloatArray> precalculated;
    FloatArray scores;
    if(use_mask)
    {
        calculate_mask_score(nai_rt, nai_bd, dmso_rt, dmso_bd, scores, mask, param);
    }else{
        calculate_score(nai_rt, nai_bd, dmso_rt, dmso_bd, scores, param);
    }
    //calculate_score(nai_rt, nai_bd, dmso_rt, dmso_bd, scores, param);
    for(const float &cur_score: scores)
    {
        FloatArray fa;
        fa.push_back(cur_score);
        precalculated.push_back(fa);
    }

    uLONG slip_end = end + wsize - 1;
    if(slip_end > NAI_RT.size()-1)
        slip_end = NAI_RT.size() - 1;

    while(cur_i<=end)
    {
        for(long i=0; i<wstep and cur_i<=end; i++)
        {
            uLONG index = index_array.front();
            
            if(index>=start and index<=end)
            {
                unsigned short valid_count;
                float shape_score = mean_score(precalculated.front(), valid_count);
                if(DMSO_BD[index] < min_cov)
                    shape_score = null;
                score[ index ] = shape_score;

                string base;
                bool masked = true;
                if(use_mask)
                {
                    base = base+chr_seq[index-1]+"\t";
                    masked = chr_mask[index];
                }

                if(masked and (NAI_RT[index] >= 1 or DMSO_RT[index] >= 1 or DMSO_BD[index] >= out_min_cov))
                {
                    OUT << chr_id << '\t' << strand_char << '\t' << index << '\t' << base
                        << NAI_RT[index] << '\t' << NAI_BD[index] << '\t'
                        << DMSO_RT[index] << '\t' << DMSO_BD[index] << '\t'
                        << shape_score << '\t' << valid_count << '\t';

                    for(float v: precalculated.front())
                        OUT << v << ',';

                    OUT << '\n';
                }
            }

            index_array.pop_front();
            precalculated.pop_front();
            
            nai_rt.pop_front();
            nai_bd.pop_front();
            dmso_rt.pop_front();
            dmso_bd.pop_front();
            if(use_mask)
                mask.pop_front();
            
            nai_rt.push_back(NAI_RT[cur_i]);
            nai_bd.push_back(NAI_BD[cur_i]);
            dmso_rt.push_back(DMSO_RT[cur_i]);
            dmso_bd.push_back(DMSO_BD[cur_i]);
            if(use_mask)
                mask.push_back(chr_mask[cur_i]);

            index_array.push_back(cur_i);

            ++cur_i;
        }

        uLONG minBD=10;
        if(low_exp_region(dmso_bd, minBD))
        {
            wstep = wsize/2;
            scores.assign(nai_rt.size(), null);
        }
        else
        {
            wstep = org_wstep;
            if(use_mask)
            {
                calculate_mask_score(nai_rt, nai_bd, dmso_rt, dmso_bd, scores, mask, param);
            }else{
                calculate_score(nai_rt, nai_bd, dmso_rt, dmso_bd, scores, param);
            }
            //calculate_score(nai_rt, nai_bd, dmso_rt, dmso_bd, scores, param);
        }

        for(uLONG i=0; i<scores.size(); i++)
        {
            if(i < precalculated.size())
                precalculated[i].push_back(scores[i]);
            else{
                FloatArray fa;
                fa.push_back(scores[i]);
                precalculated.push_back(fa);
            }
        }
    }

    while(index_array.size() > 0)
    {
        uLONG index = index_array.front();

        if(index>=start and index<=end)
        {
            unsigned short valid_count;
            float shape_score = mean_score(precalculated.front(), valid_count);
                if(DMSO_BD[index] < min_cov)
                    shape_score = null;
                score[ index ] = shape_score;
            
            string base;
            bool masked = true;
            if(use_mask)
            {
                base = base+chr_seq[index-1]+"\t";
                masked = chr_mask[index];
            }

            if(masked and (NAI_RT[index] >= 1 or DMSO_RT[index] >= 1 or DMSO_BD[index] >= out_min_cov))
            {
                OUT << chr_id << '\t' << strand_char << '\t' << index << '\t' << base
                    << NAI_RT[index] << '\t' << NAI_BD[index] << '\t'
                    << DMSO_RT[index] << '\t' << DMSO_BD[index] << '\t'
                    << shape_score << '\t' << valid_count << '\t';

                for(float v: precalculated.front())
                    OUT << v << ',';

                OUT << '\n';
            }
        }

        index_array.pop_front();
        precalculated.pop_front();
    }
}

// calculate SHAPE score with NAI + DMSO RT/BD
void calculate_mask_score(const deque<float> &nai_rt, const deque<float> &nai_bd, 
                    const deque<float> &dmso_rt, const deque<float> &dmso_bd, 
                    FloatArray &scores, const deque<bool> &mask, const icSHAPE_Param &param)
{
    scores.clear();

    deque<float> nai_rt_mask, nai_bd_mask, dmso_rt_mask, dmso_bd_mask;
    for(int i=0;i<mask.size();i++)
    {
        if(mask[i])
        {
            nai_rt_mask.push_back(nai_rt[i]);
            nai_bd_mask.push_back(nai_bd[i]);
            dmso_rt_mask.push_back(dmso_rt[i]);
            dmso_bd_mask.push_back(dmso_bd[i]);
        }
    }

    /* 2020-04-13
    if(nai_rt_mask.size() < 20)
    {
        scores.assign(nai_rt.size(), null);
        return;
    }
    */

    //FloatArray
    // Modify in 2020-04-13
    if(not valid_cov(dmso_bd_mask, nai_rt_mask, param.min_cov, 0.0))
    {
        scores.assign(nai_rt.size(), null);
        return;
    }

    float nai_rt_sf = calcScalingFactor(nai_rt_mask, param);
    float nai_bd_sf = calcScalingFactor(nai_bd_mask, param);
    float dmso_rt_sf = calcScalingFactor(dmso_rt_mask, param);
    float dmso_bd_sf = calcScalingFactor(dmso_bd_mask, param);

    /* 2020-04-13
    if(nai_rt_sf == 0 or nai_bd_sf == 0 or dmso_rt_sf == 0 or dmso_bd_sf == 0)
    {
        scores.assign(nai_rt.size(), null);
        return;
    }
    */
    if(nai_rt_sf == 0) nai_rt_sf = 1;
    if(nai_bd_sf == 0) nai_bd_sf = 1;
    if(dmso_rt_sf == 0) dmso_rt_sf = 1;
    if(dmso_bd_sf == 0) dmso_bd_sf = 1;

    FloatArray scores_mask;
    calcEnrich(dmso_bd_mask, dmso_rt_mask, nai_bd_mask, nai_rt_mask, dmso_bd_sf, dmso_rt_sf, 
                nai_bd_sf, nai_rt_sf, scores_mask, param);

    winsorization(scores_mask, param.winsor_factor);

    for(int i=0,j=0;i<mask.size();i++)
    {
        if(mask[i])
        {
            scores.push_back(scores_mask[j++]);
        }else{
            scores.push_back(null);
        }
    }
}

// calculate SHAPE score with NAI + DMSO RT/BD
void calculate_score(const deque<float> &nai_rt, const deque<float> &nai_bd, 
                    const deque<float> &dmso_rt, const deque<float> &dmso_bd, 
                    FloatArray &scores, const icSHAPE_Param &param)
{
    scores.clear();

    //FloatArray
    if(not valid_cov(dmso_bd, nai_rt))
    {
        scores.assign(nai_rt.size(), null);
        return;
    }

    float nai_rt_sf = calcScalingFactor(nai_rt, param);
    float nai_bd_sf = calcScalingFactor(nai_bd, param);
    float dmso_rt_sf = calcScalingFactor(dmso_rt, param);
    float dmso_bd_sf = calcScalingFactor(dmso_bd, param);

    /*
    if(nai_rt_sf == 0 or nai_bd_sf == 0 or dmso_rt_sf == 0 or dmso_bd_sf == 0)
    {
        scores.assign(nai_rt.size(), null);
        return;
    }
    */

    if(nai_rt_sf == 0) nai_rt_sf = 1;
    if(nai_bd_sf == 0) nai_bd_sf = 1;
    if(dmso_rt_sf == 0) dmso_rt_sf = 1;
    if(dmso_bd_sf == 0) dmso_bd_sf = 1;

    calcEnrich(dmso_bd, dmso_rt, nai_bd, nai_rt, dmso_bd_sf, dmso_rt_sf, 
                nai_bd_sf, nai_rt_sf, scores, param);

    winsorization(scores, param.winsor_factor);
}

// calculate enrichment with DMSO and NAI
void calcEnrich(const deque<float> &dmso_bd, const deque<float> &dmso_rt,
                const deque<float> &nai_bd, const deque<float> &nai_rt,
                const float &dmso_bd_sf, const float &dmso_rt_sf, 
                const float &nai_bd_sf, const float &nai_rt_sf, 
                FloatArray &score, const icSHAPE_Param &param)
{
    float sub_fac = param.sub_factor;
    float div_fac = param.div_factor;
    float add_fac = param.add_factor;

    score.clear();
    const uLONG length = dmso_bd.size();

    if(param.enrich_method==SUBSTRACTION)
    {
        for(uLONG i=0; i<length; i++)
            score.push_back( (nai_rt[i]/nai_rt_sf)-sub_fac*(dmso_rt[i]/dmso_rt_sf) );
    }else if(param.enrich_method==DIVIDING)
    {
        for(uLONG i=0; i<length; i++)
            if(dmso_bd[i] > 0)
                score.push_back( div_fac * (nai_bd[i]/nai_bd_sf) / (dmso_bd[i]/dmso_bd_sf) );
            else
                score.push_back(null);

    }else if(param.enrich_method==COMPLEX){
        for(uLONG i=0; i<length; i++)
            if(dmso_bd[i] > 0)
                score.push_back( div_fac * ((nai_rt[i]/nai_rt_sf)-sub_fac*(dmso_rt[i]/dmso_rt_sf)) / (dmso_bd[i]/dmso_bd_sf) );
            else
                score.push_back(null);
    }else if(param.enrich_method==LOG){
        for(uLONG i=0; i<length; i++)
            if(dmso_bd[i] > 0)
            {
                //cout << "add_fac: " << add_fac << endl;
                //cout << "(nai_rt[i]+add_fac)/(dmso_rt[i]+add_fac): " << (nai_rt[i]+add_fac)/(dmso_rt[i]+add_fac) << endl;
                score.push_back( log2( (nai_rt[i]+add_fac)/(dmso_rt[i]+add_fac) ) );
            }
            else
                score.push_back(null);
    }else{
        throw Unexpected_Error("enrich_method unrecognized option");
    }
}

//###############################
//###############################
//    Insert 1 /**** smart-SHAPE fansion (with NAI only) ****/
//###############################
//###############################

// calculate SHAPE score with NAI RT/BD in non-junction regions
void sliding_non_junction(  const uIntArray &NAI_RT, const uIntArray &NAI_BD,
                            const JunctionArray &junctions, FloatArray &score, 
                            const string &chr_id, const STRAND &strand, 
                            ofstream &OUT, const smartSHAPE_Param &param,
                            const vector<bool> &chr_mask, const string &chr_seq)
{
    bool use_mask = chr_mask.empty() ? false : true;
    if(use_mask and NAI_RT.size() != chr_mask.size())
        throw runtime_error("NAI_RT.size() != chr_mask.size() -- "+to_string(NAI_RT.size())+" != "+to_string(chr_mask.size()));

    const char strand_char = (strand == POSITIVE) ? '+' : '-';
    const uINT out_min_cov = param.out_min_cov;
    const uINT min_cov = param.min_cov;
    const uINT wsize = param.wsize;
    const uINT org_wstep = param.wstep;
    uINT wstep = org_wstep;

    deque<bool> mask;     // =========
    deque<float> nai_rt;
    deque<float> nai_bd;

    const uLONG chr_size = NAI_RT.size() - 1;
    if(NAI_BD.size() != NAI_RT.size())
    {
        cerr << RED << "FATAL Error: RT,BD different length" << DEF << endl;
        exit(-1);
    }

    //////////////////////////////////////////Only useful when non-junction///////////////////////////////////////////////////
    if(param.no_sliding)
    {

        if(chr_size < 50)
        {
            cerr << RED << "Warning: Chromosome " << chr_id << "(" << strand_char << ")" << " is too short. " << chr_size;
            cerr << "Skip this one" << DEF << endl;
            return;
        }

        FloatArray scores;
        for(uLONG i=0; i<chr_size-30; i++)
        {
            if(use_mask)
                mask.push_back(chr_mask[i]);
            nai_rt.push_back(NAI_RT[i]);
            nai_bd.push_back(NAI_BD[i]);
        }

        if(use_mask)
            calculate_mask_score(nai_rt, nai_bd, scores, mask, param);
        else
            calculate_score(nai_rt, nai_bd, scores, param);

        for(uLONG i=6; i<chr_size-30; i++)
            if(NAI_BD[i] < min_cov)
                score[i] = null;
            else
                score[i] = scores[i];

        for(uLONG index=1; index<chr_size; index++)
        {
            string base;
            bool masked = true;
            if(use_mask)
            {
                base = base+chr_seq[index-1]+"\t";
                masked = chr_mask[index];
            }

            if(masked and (NAI_RT[index] >= 1 or NAI_BD[index] >= out_min_cov))
            {
                OUT << chr_id << '\t' << strand_char << '\t' << index << '\t' << base
                    << NAI_RT[index] << '\t' << NAI_BD[index] << '\t'
                    << score[index] << '\t' << (score[index]==null ? 0 : 1) << '\t' 
                    //<< ((index>=scores.size()) ? null : scores[index]) << '\n';
                    << score[index] << '\n';
            }

        }

        return;
    }

    if(chr_size < wsize*2)
    {
        cerr << RED << "Warning: Chromosome " << chr_id << "(" << strand_char << ")" << " is too short. " << chr_size;
        cerr << "Not use sliding window (-non-sliding) or use a shorter window size (-wsize) is recommend for this one" << DEF << endl;
        return;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    uINT ji = 0;
    uLONG juncNum = junctions.size();

    bool has_junction = (junctions.size() >= 1);

    deque<uLONG> index_array;

    // Initiazation
    uLONG cur_i = 1;
    // move to start
    while(NAI_RT[cur_i] == 0 and cur_i <= chr_size)
        ++cur_i;
    
    if(has_junction and ji<juncNum)
        while(junctions[ji].first < cur_i)
            ++ji;

    uLONG de_len = 0;
    while(de_len < wsize and cur_i<=chr_size)
    {
        if(has_junction and ji<juncNum)
            if(junctions[ji].first == cur_i)
                cur_i = junctions[ji++].second+1;

        nai_rt.push_back(NAI_RT[cur_i]);
        nai_bd.push_back(NAI_BD[cur_i]);
        if(use_mask)
            mask.push_back(chr_mask[cur_i]);

        index_array.push_back(cur_i);

        ++cur_i;
        ++de_len;
    }

    deque<FloatArray> precalculated;
    FloatArray scores;
    if(use_mask)
    {
        calculate_mask_score(nai_rt, nai_bd, scores, mask, param);
    }else{
        calculate_score(nai_rt, nai_bd, scores, param);
    }
    
    for(const float &score: scores)
    {
        FloatArray fa;
        fa.push_back(score);
        precalculated.push_back(fa);
    }

    //clog << "start to loding ...." << endl;
    while(cur_i<=chr_size)
    {
        for(long i=0; i<wstep and cur_i<=chr_size; i++)
        {
            if(cur_i % 1000000 == 0)
            //if(cur_i % 1000 == 0)
                clog << "\t lines: " << cur_i << "(" << 100.0*cur_i/chr_size << "%) wStep: " << wstep << endl;

            uLONG index = index_array.front();
            unsigned short valid_count;
            float shape_score = mean_score(precalculated.front(), valid_count);
            if(NAI_BD[index] < min_cov)
                shape_score = null;
            score[ index ] = shape_score;

            string base;
            bool masked = true;
            if(use_mask)
            {
                base = base+chr_seq[index-1]+"\t";
                masked = chr_mask[index];
            }

            if(masked and (NAI_RT[index] >= 1 or NAI_BD[index] >= out_min_cov))
            {
                OUT << chr_id << '\t' << strand_char << '\t' << index << '\t' << base
                    << NAI_RT[index] << '\t' << NAI_BD[index] << '\t'
                    << shape_score << '\t' << valid_count << '\t';

                for(float v: precalculated.front())
                    OUT << v << ',';

                OUT << '\n';
            }


            index_array.pop_front();
            precalculated.pop_front();
            
            nai_rt.pop_front();
            nai_bd.pop_front();
            if(use_mask)
                mask.pop_front();

            if(has_junction and ji<juncNum)
                if(junctions[ji].first == cur_i)
                    cur_i = junctions.at(ji++).second+1;

            nai_rt.push_back(NAI_RT[cur_i]);
            nai_bd.push_back(NAI_BD[cur_i]);
            if(use_mask)
                mask.push_back(chr_mask[cur_i]);

            index_array.push_back(cur_i);

            ++cur_i;
        }


        uLONG minBD=10;
        if(low_exp_region(nai_bd, minBD))
        {
            wstep = wsize - 1;
            scores.assign(nai_rt.size(), null);
        }
        else
        {
            wstep = org_wstep;
            if(use_mask)
            {
                calculate_mask_score(nai_rt, nai_bd, scores, mask, param);
            }else{
                calculate_score(nai_rt, nai_bd, scores, param);
            }
        }

        for(uLONG i=0; i<scores.size(); i++)
        {
            if(i < precalculated.size())
                precalculated[i].push_back(scores[i]);
            else{
                FloatArray fa;
                fa.push_back(scores[i]);
                precalculated.push_back(fa);
            }
        }

    }

    if(has_junction)
        if(ji <= junctions.size()-1)
            cerr << RED << "Error: Junction index not reach final " << junctions.at(ji).first << "-" << junctions.at(ji).second << DEF << endl;

    //clog << "finish loding ...." << endl;
    while(precalculated.size() > 0)
    {
        uLONG index = index_array.front();
        unsigned short valid_count;
        float shape_score = mean_score(precalculated.front(), valid_count);
        if(NAI_BD[index] < min_cov)
            shape_score = null;
        score[ index ] = shape_score;
        
        string base;
        bool masked = true;
        if(use_mask)
        {
            base = base+chr_seq[index-1]+"\t";
            masked = chr_mask[index];
        }

        if(masked and (NAI_RT[index] >= 1 or NAI_BD[index] >= out_min_cov))
        {
            OUT << chr_id << '\t' << strand_char << '\t' << index << '\t' << base
                << NAI_RT[index] << '\t' << NAI_BD[index] << '\t'
                << shape_score << '\t' << valid_count << '\t';

            for(float v: precalculated.front())
                OUT << v << ',';

            OUT << '\n';
        }

        index_array.pop_front();
        precalculated.pop_front();
    }

    for(uLONG i=0; i<30 and i<=chr_size; i++)
        score[i] = null;
    uLONG s = 0;
    if(chr_size > 30)
        s = chr_size - 30;
    for(uLONG i=s; i<=chr_size; i++)
        score[i] = null;
}

// sliding all junctions in junction regions with NAI RT/BD
void sliding_junction( const uIntArray &NAI_RT, const uIntArray &NAI_BD,
                        const JunctionArray &junctions, FloatArray &score, 
                        const string &chr_id, const STRAND &strand, 
                        ofstream &OUT, const smartSHAPE_Param &param,
                        const vector<bool> &chr_mask, const string &chr_seq)
{

    ///////////////////////////
    if(param.no_sliding)
        return;
    ///////////////////////////

    if(NAI_RT.size() != NAI_BD.size())
    {
        cerr << RED << "FATAL Error: RT,BD different length" << DEF << endl;
        exit(-1);
    }

    uLONG junc_num = junctions.size();

    uLONG index = 0;
    for(const Junction &junc: junctions)
    {
        if(++index % 100 == 0)
            clog << "\t junc ratio: " << 100.0*index/junc_num << "%" << endl;
        sliding_single_junction(NAI_RT, NAI_BD, junc.first, junc.second, score, chr_id, strand, OUT, param, chr_mask, chr_seq);
    }
}


// sliding a single junction with NAI RT/BD
void sliding_single_junction(const uIntArray &NAI_RT, const uIntArray &NAI_BD,
                            const uLONG &start, const uLONG &end, FloatArray &score,
                            const string &chr_id, const STRAND &strand, 
                            ofstream &OUT, const smartSHAPE_Param &param,
                            const vector<bool> &chr_mask, const string &chr_seq)
{
    bool use_mask = chr_mask.empty() ? false : true;
    if(use_mask and NAI_RT.size() != chr_mask.size())
        throw runtime_error("NAI_RT.size() != chr_mask.size() -- "+to_string(NAI_RT.size())+" != "+to_string(chr_mask.size()));

    const char strand_char = (strand == POSITIVE) ? '+' : '-';
    const uINT out_min_cov = param.out_min_cov;
    const uINT min_cov = param.min_cov;
    const uINT wsize = param.wsize;
    const uINT org_wstep = param.wstep;
    uINT wstep = org_wstep;

    deque<bool> mask;
    deque<float> nai_rt;
    deque<float> nai_bd;

    deque<uLONG> index_array;


    uLONG cur_i = 1;
    if(start > wsize)
        cur_i = start - wsize + 1;

    uLONG de_len = 0;
    while(de_len < wsize and cur_i<=end)
    {
        nai_rt.push_back(NAI_RT[cur_i]);
        nai_bd.push_back(NAI_BD[cur_i]);
        if(use_mask)
            mask.push_back(chr_mask[cur_i]);

        index_array.push_back(cur_i);

        ++cur_i;
        ++de_len;
    }

    deque<FloatArray> precalculated;
    FloatArray scores;
    //calculate_score(nai_rt, nai_bd, scores, param);
    if(use_mask)
    {
        calculate_mask_score(nai_rt, nai_bd, scores, mask, param);
    }else{
        calculate_score(nai_rt, nai_bd, scores, param);
    }
    for(const float &score: scores)
    {
        FloatArray fa;
        fa.push_back(score);
        precalculated.push_back(fa);
    }

    uLONG slip_end = end + wsize - 1;
    if(slip_end > NAI_RT.size()-1)
        slip_end = NAI_RT.size() - 1;

    while(cur_i<=slip_end)
    {
        for(uLONG i=0; i<wstep and cur_i<=slip_end; i++)
        {
            uLONG index = index_array.front();

            if(index>=start and index<=end)
            {
                unsigned short valid_count;
                float shape_score = mean_score(precalculated.front(), valid_count);
                if(NAI_BD[index] < min_cov)
                    shape_score = null;
                
                score[ index ] = shape_score;

                string base;
                bool masked = true;
                if(use_mask)
                {
                    base = base+chr_seq[index-1]+"\t";
                    masked = chr_mask[index];
                }

                if(masked and (NAI_RT[index] >= 1 or NAI_BD[index] >= out_min_cov))
                {
                    OUT << chr_id << '\t' << strand_char << '\t' << index << '\t' << base
                        << NAI_RT[index] << '\t' << NAI_BD[index] << '\t'
                        << shape_score << '\t' << valid_count << '\t';

                    for(float v: precalculated.front())
                        OUT << v << ',';

                    OUT << '\n';
                }
            }
            index_array.pop_front();
            precalculated.pop_front();
            
            nai_rt.pop_front();
            nai_bd.pop_front();
            if(use_mask)
                mask.pop_front();

            nai_rt.push_back(NAI_RT[cur_i]);
            nai_bd.push_back(NAI_BD[cur_i]);
            if(use_mask)
                mask.push_back(chr_mask[cur_i]);

            index_array.push_back(cur_i);

            ++cur_i;
        }

        uLONG minBD=10;
        if(low_exp_region(nai_bd, minBD))
        {
            wstep = wsize/2;
            scores.assign(nai_rt.size(), null);
        }
        else
        {
            wstep = org_wstep;
            if(use_mask)
            {
                calculate_mask_score(nai_rt, nai_bd, scores, mask, param);
            }else{
                calculate_score(nai_rt, nai_bd, scores, param);
            }
            //calculate_score(nai_rt, nai_bd, scores, param);
        }

        for(uLONG i=0; i<scores.size(); i++)
        {
            if(i < precalculated.size())
                precalculated.at(i).push_back(scores[i]);
            else{
                FloatArray fa;
                fa.push_back(scores[i]);
                precalculated.push_back(fa);
            }
        }
    }

    while(precalculated.size() > 0)
    {
        uLONG index = index_array.front();

        if(index>=start and index<=end)
        {
            unsigned short valid_count;
            float shape_score = mean_score(precalculated.front(), valid_count);
            if(NAI_BD[index] < min_cov)
                shape_score = null;
            score[ index ] = shape_score;
            

            string base;
            bool masked = true;
            if(use_mask)
            {
                base = base+chr_seq[index-1]+"\t";
                masked = chr_mask[index];
            }

            if(masked and (NAI_RT[index] >= 1 or NAI_BD[index] >= out_min_cov))
            {
                OUT << chr_id << '\t' << strand_char << '\t' << index << '\t' << base
                    << NAI_RT[index] << '\t' << NAI_BD[index] << '\t'
                    << shape_score << '\t' << valid_count << '\t';

                for(float v: precalculated.front())
                    OUT << v << ',';

                OUT << '\n';
            }
        }

        index_array.pop_front();
        precalculated.pop_front();
    }
}



// calculate smart-SHAPE score with NAI RT and BD
void calculate_mask_score(const deque<float> &nai_rt, const deque<float> &nai_bd, FloatArray &scores, const deque<bool> &mask, const smartSHAPE_Param &param)
{
    scores.clear();

    deque<float> nai_rt_mask, nai_bd_mask;
    for(int i=0;i<mask.size();i++)
    {
        if(mask[i])
        {
            nai_rt_mask.push_back(nai_rt[i]);
            nai_bd_mask.push_back(nai_bd[i]);
        }
    }

    /* 2020-04-13
    if(nai_rt_mask.size() < 20)
    {
        scores.assign(nai_rt.size(), null);
        return;
    }
    */

    //FloatArray
    // Modify in 2020-04-13
    if(not valid_cov(nai_bd_mask, nai_rt_mask, param.min_cov, 0.0))
    {
        scores.assign(nai_rt.size(), null);
        return;
    }

    if(nai_bd_mask.size() != nai_bd_mask.size())
    {
        cerr << RED << "Error: nai_bd_mask.size() != nai_bd_mask.size()" << DEF << endl;
        exit(-1);
    }

    FloatArray scores_mask;
    calcEnrich(nai_bd_mask, nai_rt_mask, scores_mask, param);
    sink_score(scores_mask, param.sink_index);
    winsorization(scores_mask, param.winsor_factor);

    for(int i=0,j=0;i<mask.size();i++)
    {
        if(mask[i])
        {
            scores.push_back(scores_mask[j++]);
        }else{
            scores.push_back(null);
        }
    }
}


// calculate smart-SHAPE score with NAI RT and BD
void calculate_score(const deque<float> &nai_rt, const deque<float> &nai_bd, FloatArray &scores, const smartSHAPE_Param &param)
{
    scores.clear();

    if(not valid_cov(nai_bd, nai_rt))
    {
        scores.assign(nai_rt.size(), null);
        return;
    }

    if(nai_rt.size() != nai_bd.size())
    {
        cerr << RED << "Error: nai_rt.size() != nai_bd.size()" << DEF << endl;
        exit(-1);
    }

    calcEnrich(nai_bd, nai_rt, scores, param);
    sink_score(scores, param.sink_index);
    winsorization(scores, param.winsor_factor);
}


// calculate enrichment with only NAI
void calcEnrich(const deque<float> &nai_bd, const deque<float> &nai_rt, FloatArray &score, const smartSHAPE_Param &param)
{
    //float sub_fac = param.sub_factor;
    //float div_fac = param.div_factor;
    float add_fac = param.add_factor;

    const uLONG length = nai_bd.size();

    if(param.enrich_method==RT)
    {
        for(uLONG i=0; i<length; i++)
            if(nai_bd[i] > 0)
                score.push_back( nai_rt[i] );
            else
                score.push_back( null );
    }else if(param.enrich_method==DIVIDING)
    {
        for(uLONG i=0; i<length; i++)
            if(nai_bd[i] > 0)
                score.push_back( nai_rt[i]/nai_bd[i] );
            else
                score.push_back(null);

    }else if(param.enrich_method==LOG){
        for(uLONG i=0; i<length; i++)
            if(nai_bd[i] > 0)
                score.push_back( log2(nai_rt[i]+add_fac) );
            else
                score.push_back(null);
    }else{
        throw Unexpected_Error("enrich_method unrecognized option");
    }
}

// sink all SHAPE that means xi <- xi - sorted_scores[sink_index*len]
void sink_score(FloatArray &score, const float &sink_index)
{
    if(sink_index == 0)
        return;

    FloatArray sortedTrimed_array = score;

    sort(sortedTrimed_array.begin(), sortedTrimed_array.end());

    uLONG len = sortedTrimed_array.size();
    uLONG winsorLen = sink_index * len;

    float sink_value = sortedTrimed_array.at(winsorLen);

    for(float &v: score)
        v -= sink_value;
}

/**** statistics functions ****/

// winsorization a data array that means all data are assigned to 0-1
void winsorization( FloatArray &score, const float &winsor_factor)
{

    float winsorUpper, winsorLower;
    winsorWindow(score, winsor_factor, winsorLower, winsorUpper);

    if(winsorLower < 0) winsorLower = 0;

    if(winsorLower >= winsorUpper)
    {
        for(uLONG i=0; i<score.size(); i++)
            score[i] = null;

        return;
    }

    for(uLONG i=0; i<score.size(); i++)
    {
        if(score[i] != null)
        {
            if(score[i] > winsorUpper)
                score[i] = 1;
            else if(score[i] < winsorLower)
                score[i] = 0;
            else
                score[i] = (score[i] - winsorLower)/(winsorUpper-winsorLower);
        }
    }
}

// Get winsor upper(U) and lower(L), and normalize each raw xi to (xi-L)/(U-L)
void winsorWindow(  const FloatArray &score, const float &winsor_factor, float &winsorLower, float &winsorUpper)
{

    FloatArray sortedTrimed_array;
    for(uINT i=0; i<=score.size(); i++)
    {
        if(score[i] != null)
            sortedTrimed_array.push_back(score[i]);
    }

    sort(sortedTrimed_array.begin(), sortedTrimed_array.end());

    uLONG len = sortedTrimed_array.size();
    if(len < 20)
    {
        winsorLower = winsorUpper = 0;
        return;
    }

    uLONG winsorLen = winsor_factor * len;

    winsorLower = sortedTrimed_array[winsorLen];
    winsorUpper = sortedTrimed_array[len-winsorLen-1];
}


// calculate calcScalingFactor when input a array
float calcScalingFactor(const deque<float> &data_array, const icSHAPE_Param &param)
{
    FloatArray sorted_data(data_array.cbegin(), data_array.cend());
    sort(sorted_data.rbegin(), sorted_data.rend());

    uINT len = sorted_data.size();
    FloatArray rangeOfSelection;
    if(param.norm_sample_method == SMART)
    {
        for(uLONG idx=0; idx<len; idx++)
        {
            if(sorted_data[idx] <= 0) break;
            rangeOfSelection.push_back(sorted_data[idx]);
        }
    }else{
        uINT selectStart = 0;
        uINT selectEnd = len - 1;

        float f = param.norm_sample_factor;

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
            throw Unexpected_Error("norm_method unrecognized option");
        }

        for(auto it=sorted_data.cbegin()+selectStart; it!=sorted_data.cbegin()+selectEnd+1; it++)
            rangeOfSelection.push_back(*it);
    }

    len = rangeOfSelection.size();

    float scalling_factor = 1.0;
    if(param.norm_calc_method == MEDIAN)
    {
        float median = 1;
        if ( len % 2 == 0 ) {  median = (rangeOfSelection[len/2-1] + rangeOfSelection[len/2]) /2;  }
        else {  median = rangeOfSelection[(len-1)/2];  }
        scalling_factor = median;

    }else if(param.norm_calc_method == MEAN)
    {
        float Sum = 0.0;
        for(const float &v: rangeOfSelection) Sum += v;
        scalling_factor = Sum/rangeOfSelection.size();

    }else if(param.norm_calc_method == PEAK)
    {
        scalling_factor = rangeOfSelection[0];
    }

    return scalling_factor;
}

// Sync tab files (from sam2tab), must be sorted
bool sync_chrs(vector<ifstream *> &i_vec, StringArray &chr_ids, vector< vector<Map_Record> > &record_array)
{
    uLONG file_num = i_vec.size();

    if(file_num != chr_ids.size() or file_num != record_array.size())
        throw Unexpected_Error("FATAL Error: sync_chrs different size");

    for(uLONG i=0; i<file_num; i++)
        if(not read_chr(record_array[i], i_vec[i], chr_ids[i]))
            return false;
    
    string max_chr = chr_ids[0];
    bool sync = false;
    while(not sync)
    {
        sync = true;
        for(uLONG i=0; i<file_num; i++)
        {
            if(chr_ids[i] > max_chr)
            {
                max_chr = chr_ids[i];
                sync = false;
            }
            else if(chr_ids[i] < max_chr)
            {
                if(not read_chr(record_array[i], i_vec[i], chr_ids[i]))
                    return false;
                sync = false;
            }
        }
    }
    return true;
}

// Check if input handle opened
void check_input_handle(ifstream &IN, const string &fn)
{
    if(not IN)
    {
        cerr << RED << "FATAL Error: " << fn << "cannot be read" << endl;
        exit(-1);
    }
}

// Check if output handle opened
void check_output_handle(ofstream &OUT, const string &fn)
{
    if(not OUT)
    {
        cerr << RED << "FATAL Error: " << fn << "cannot be write" << endl;
        exit(-1);
    }
}



