
#ifndef PARIS_H
#define PARIS_H

#include "pan_type.h"
#include "sstructure.h"
#include "exceptions.h"
#include "align.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include "string_split.h"
#include <exception>
#include <string>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>

#include "sam.h"

namespace pan{

using std::string;
using std::vector;
using std::unordered_map;
using std::istream;
using std::ostream;
using std::istringstream;
using std::pair;
using std::ifstream;
using std::runtime_error;
using std::min;
using std::max;
using std::cout;
using std::endl;


struct Duplex_Hang
{
    string read_id;

    string chr_id_1;
    char strand_1 = '+';
    uINT flag_1 = 0;
    string cigar_1;
    uLONG start_1 = 0;
    uLONG end_1 = 0;
    
    string chr_id_2;
    char strand_2 = '+';
    uINT flag_2 = 0;
    string cigar_2;
    uLONG start_2 = 0;
    uLONG end_2 = 0;

    Duplex_Hang(const string &read_id, 
        const string &chr_id_1, const char &strand_1, const uINT &flag_1, const string &cigar_1, const uLONG &start_1, const uLONG &end_1,
        const string &chr_id_2, const char &strand_2, const uINT &flag_2, const string &cigar_2, const uLONG &start_2, const uLONG &end_2):
    read_id(read_id), chr_id_1(chr_id_1), strand_1(strand_1), flag_1(flag_1), cigar_1(cigar_1), start_1(start_1), end_1(end_1),
    chr_id_2(chr_id_2), strand_2(strand_2), flag_2(flag_2), cigar_2(cigar_2), start_2(start_2), end_2(end_2) 
    {  
        this->read_id.shrink_to_fit();
        this->chr_id_1.shrink_to_fit();
        this->chr_id_2.shrink_to_fit();
        this->cigar_1.shrink_to_fit();
        this->cigar_2.shrink_to_fit();
    }

    Duplex_Hang() = default;

};

ostream &operator<<(ostream &OUT, const Duplex_Hang &dh);
ostream &operator<<(ostream &OUT, const vector<Duplex_Hang> &dh_array);
bool operator<(const Duplex_Hang &dh_1, const Duplex_Hang &dh_2);


/*  build a vector<Duplex_Hang> from a sam file
    
*/
uLONG get_chromosome_hang(
            const string &sam_file_name, 
            vector<Duplex_Hang> &dh_array, 
            uINT min_overhang, 
            uINT min_readlen,
            const string &chr_id,
            const char strand='+');

void read_dh_from_sam(const string &sam_file_name, vector<Duplex_Hang> &dh_array);

/*  init and fill a symmetric matrix with vector<Duplex_Hang>

*/
template<typename T>
void fill_sym_matrix(Matrix<T> &matrix, 
                    const vector<Duplex_Hang> &dh_array, 
                    uLONG chr_len);


/*  compute a block feature from a Matrix

*/

enum FEATURE
{
    FEATURE_MEAN, FEATURE_MAX, FEATURE_MIN, FEATURE_SUM
};

template<typename T>
T matrix_block_feature( const Matrix<T> &matrix,
                        typename Matrix<T>::size_type x,
                        typename Matrix<T>::size_type y,
                        typename Matrix<T>::size_type x_len,
                        typename Matrix<T>::size_type y_len,
                        FEATURE feature);

template<typename T1, typename T2>
void compress_matrix(const Matrix<T1> &raw_matrix,
                    Matrix<T2> &target_matrix,
                    uLONG target_size,
                    FEATURE feature);

void compress_regions(const RegionArray &raw_regions,
                    RegionArray &target_regions,
                    const uLONG raw_size,
                    const uLONG target_size);

/*  Interaction Region Border and Support Reads

*/

struct InterRegion {
    InterRegion(pair<Region,Region> region, double max_cov, Point point): region(region),max_cov(max_cov),max_point(point) {}

    pair<Region,Region> region;

    double max_cov = 0;
    Point max_point;
};

ostream& operator<<(ostream& OUT, const InterRegion& iter_reg);
ostream& operator<<(ostream& OUT, const vector<InterRegion>& iter_reg);
bool operator<(const InterRegion &inter_1, const InterRegion &inter_2);
pair<uLONG, uLONG> overlap(const InterRegion &iter_1, const InterRegion &iter_2);

template<typename T>
void scan_interaction(  const Matrix<T> raw_matrix, 
                        vector<InterRegion> &interact_regions,
                        uLONG min_dist, 
                        uLONG min_window_size, 
                        uLONG max_window_size, 
                        T percep_threshold, 
                        T extend_threshold);

template<typename T>
uLONG read_matrix(const string &matrix_file, Matrix<T> &matrix);

vector<Point> search_TT_cross_linking(const SStructure &structure);


template<typename T>
void fill_sym_matrix_with_aligned_dh(   Matrix<T> &matrix,
                                    const std::vector<Duplex_Hang> &dh_array, 
                                    const Multi_Align &alignment, 
                                    const string &chromosome_name, 
                                    const uINT target_size );


void read_domain_file(const string &domain_file_name, RegionArray &regions);

/* remove PARIS background according to GRID-Seq method */
template<typename T>
void remove_paris_background(const Matrix<T> &raw_matrix, Matrix<T> &matrix, const uINT around=5);

/* get a sub-matrix from raw matrix: coordination is 0-based */
template<typename T>
void sub_matrix(const Matrix<T> &raw_matrix, Rect<T> &sub_rect, 
                uLONG x_start=0, uLONG x_end=-1UL, 
                uLONG y_start=0, uLONG y_end=-1UL);



































































/* ================= implemetation ================= */


/*  init and fill a symmetric matrix with vector<Duplex_Hang>

*/
template<typename T>
void fill_sym_matrix(Matrix<T> &matrix, 
                    const vector<Duplex_Hang> &dh_array, 
                    uLONG chr_len)
{
    const string &chr_name = dh_array[0].chr_id_1;
    const char &strand = dh_array[0].strand_1;

    init_matrix(matrix, chr_len);
    for(const Duplex_Hang &dh: dh_array)
    {
        if(dh.chr_id_1 != chr_name or dh.strand_1 != strand or dh.chr_id_2 != chr_name or dh.strand_2 != strand)
            throw Unexpected_Error("fill_sym_matrix input Duplex_Hang should be same Chromosome and same strand", true);
        for(uLONG x=dh.start_1-1;x<dh.end_1;x++)
        {
            if( dh.end_1 > chr_len )
                throw Unexpected_Error("Bad Chromosome Length");
            for(uLONG y=dh.start_2-1;y<dh.end_2;y++)
            {
                if( dh.end_2 > chr_len )
                    throw Unexpected_Error("Bad Chromosome Length");
                matrix[x][y]++;
                matrix[y][x]++;
            }
        }
    }
}



/*  compute a block feature from a Matrix

*/



template<typename T>
T matrix_block_feature( const Matrix<T> &matrix,
                        typename Matrix<T>::size_type x,
                        typename Matrix<T>::size_type y,
                        typename Matrix<T>::size_type x_len,
                        typename Matrix<T>::size_type y_len,
                        FEATURE feature)
{
    T total_v = T(), min_v = T(), max_v = T();
    for(decltype(matrix.size()) idx=x; idx<x+x_len; idx++)
        for(decltype(matrix.size()) idy=y; idy<y+y_len; idy++)
        {
            //std::cerr << "Fetch " << idx << "\t" << idy << std::endl;
            auto cur_v = matrix.at(idx).at(idy);
            total_v += cur_v;
            min_v = min(min_v, cur_v);
            max_v = max(max_v, cur_v);
            /*
            if(min_v > cur_v)
                min_v = cur_v;
            if(max_v < cur_v)
                max_v = cur_v;
            */
        }

    if(feature == FEATURE_MEAN)
        return total_v / (x_len*y_len);
    else if(feature == FEATURE_MAX)
        return max_v;
    else if(feature == FEATURE_MIN)
        return min_v;
    else if(feature == FEATURE_SUM)
        return total_v;
    else
        return 0;
}


/*  compress a Matrix

*/

template<typename T1, typename T2>
void compress_matrix(const Matrix<T1> &raw_matrix,
                    Matrix<T2> &target_matrix,
                    uLONG target_size,
                    FEATURE feature)
{
    using size_type = typename Matrix<T2>::size_type;

    if(target_size*2 >= raw_matrix.size())
    {
        throw Unexpected_Error("Invalid target_size to compress matrix");
    }
    init_matrix(target_matrix, target_size);
    double step = 1.0 * raw_matrix.size() / target_size;
    for(size_type idx=0; idx<target_size; idx++)
        for(size_type idy=0; idy<target_size; idy++)
        {
            size_type x_base = idx * step;
            size_type y_base = idy * step;
            size_type x_len = min( decltype(raw_matrix.size())((idx+1) * step), raw_matrix.size() ) - x_base;
            size_type y_len = min( decltype(raw_matrix.size())((idy+1) * step), raw_matrix.size() ) - y_base;
            T1 value = matrix_block_feature( raw_matrix, x_base, y_base, x_len, y_len, feature );
            target_matrix.at(idx).at(idy) = value;
        }
}



/*

./call_interaction -in 59.sam -chr KU501215.1 -out 59.interaction -min_armlen 15 -min_dist 100 -min_window_size 30 -max_window_size 200 -percep_threshold 30 -extend_threshold 2

*/

template<typename T>
T locate_matrix_max_point(const Matrix<T> &raw_matrix,
        Point &max_pos,
        uLONG x_lower,
        uLONG x_upper,
        uLONG y_lower,
        uLONG y_upper,
        function<bool(uLONG,uLONG)> is_valid_pos)
{
    T block_max = T();
    max_pos.first = x_lower;
    max_pos.second = y_lower;
    if(x_lower > x_upper or y_lower > y_upper)
    {
        throw range_error("matrix_block_max parameter error");
    }else if(x_upper > raw_matrix.size() or y_upper > raw_matrix.size())
    {
        throw range_error("matrix_block_max parameter error");
    }else{
        for(typename Matrix<T>::size_type idx = x_lower; idx < x_upper; idx++)
            for(typename Matrix<T>::size_type idy = y_lower; idy < y_upper; idy++)
                if( is_valid_pos(idx, idy) and block_max < raw_matrix[idx][idy])
                {
                    max_pos.first = idx;
                    max_pos.second = idy;
                    block_max = raw_matrix[idx][idy];
                }
    }
    return block_max;
}

template<typename T>
void scan_interaction(  const Matrix<T> raw_matrix, 
                        vector<InterRegion> &interact_regions,
                        const uLONG min_dist, 
                        const uLONG min_window_size, 
                        const uLONG max_window_size, 
                        const T percep_threshold, 
                        const T extend_threshold)
{
    using size_type = typename Matrix<T>::size_type;

    interact_regions.clear();
    auto matrix = raw_matrix;
    auto matrix_size = matrix.size();

    auto valid_pos_func = [min_window_size, min_dist](uLONG x, uLONG y){ return y+min_window_size+min_dist<=x ? true : false; };
    //bool combine = true;

    for(size_type idx=0; idx<matrix_size; idx++)
        for(size_type idy=0; idy<matrix_size; idy++)
            if(not valid_pos_func(idx, idy))
                matrix[idx][idy] = T();



    for(size_type x_idx=0; x_idx<matrix_size; x_idx+=min_window_size)
    {
        //for(size_type y_idx=x_idx+max(min_dist,min_window_size); y_idx<matrix_size;y_idx+=min_window_size)
        for(size_type y_idx=0; y_idx+min_window_size+min_dist<=x_idx;y_idx+=min_window_size)
        {
           // cout << "Y Step " << y_idx << endl;

            size_type width = min_window_size;
            size_type height = min_window_size;
            size_type proper_x_upper = min(x_idx+width, matrix_size);
            size_type proper_y_upper = min(y_idx+height, matrix_size);
            Point last_max_pos, max_pos;
            T current_max_rc = locate_matrix_max_point<T>(matrix, last_max_pos, x_idx, proper_x_upper, y_idx, proper_y_upper, valid_pos_func);

            if( current_max_rc >= percep_threshold )
            {
               // cout << "Find One: " << x_idx << "\t" << y_idx << endl;
                // posit center
                
                x_idx = last_max_pos.first >= min_window_size / 2 ? last_max_pos.first - min_window_size / 2 : 0;
                y_idx = last_max_pos.second >= min_window_size / 2 ? last_max_pos.second - min_window_size / 2 : 0;

                //if( not valid_pos_func(x_idx, y_idx) )
                //

                proper_x_upper = min(x_idx+width, matrix_size);
                proper_y_upper = min(y_idx+height, matrix_size);
                T max_rc_record = locate_matrix_max_point<T>(matrix, max_pos, x_idx, proper_x_upper, y_idx, proper_y_upper, valid_pos_func);
                if(max_pos.second + min_window_size + min_dist > max_pos.first)
                    max_pos = last_max_pos;

                //const uINT max_iter = 10;
                while(max_pos != last_max_pos and max_pos.second + min_window_size + min_dist <= max_pos.first)
                {
                    //std::cout << max_pos << "\t" << last_max_pos << std::endl;
                    //std::cout << x_idx << "\t" << y_idx << std::endl;
                    x_idx = last_max_pos.first >= min_window_size / 2 ? last_max_pos.first - min_window_size / 2 : 0;
                    y_idx = last_max_pos.second >= min_window_size / 2 ? last_max_pos.second - min_window_size / 2 : 0;

                    proper_x_upper = min(x_idx+width, matrix_size);
                    proper_y_upper = min(y_idx+height, matrix_size);
                    last_max_pos = max_pos;
                    max_rc_record = locate_matrix_max_point<T>(matrix, max_pos, x_idx, proper_x_upper, y_idx, proper_y_upper, valid_pos_func);
                }
                //if(max_pos.second + min_window_size + min_dist > max_pos.first)
                max_pos = last_max_pos;

/*
                if(y_idx+height+min_dist>x_idx)
                {
                    //std::cerr << "score\n";
                    auto bias = y_idx+height+min_dist-x_idx;
                    x_idx = min(bias/2+x_idx, matrix_size-1);
                    y_idx = (y_idx>bias/2) ? (y_idx-bias/2) : 0;

                    combine = true;
                    //std::cerr << x_idx << "\t" << y_idx << std::endl;
                }
*/

               // cout << "Start to extend..." << endl;
                bool extend_up(true), extend_left(true), extend_right(true), extend_down(true);
                uLONGArray tmp_vec; uLONG max_rc;
                while( extend_up or extend_left or extend_right or extend_down )
                {
                    //extend_up = extend_left = extend_right = extend_down = true;
                    // extend left
                    if(extend_left)
                    {
                        tmp_vec.clear();
                        for(size_type idx=0;idx<min(height, matrix_size-y_idx);idx++)
                            tmp_vec.push_back(matrix.at(x_idx).at(y_idx+idx));
                        max_rc = *max_element(tmp_vec.cbegin(), tmp_vec.cend());
                        //if(max_rc>=extend_threshold and y_idx-x_idx>min_dist and x_idx>0)
                        if(max_rc>=extend_threshold and y_idx+height+min_dist<x_idx and x_idx>0 and width < max_window_size)
                        {
                            x_idx--;
                            width++;
                        }else{
                            extend_left = false;
                        }
                    }
                    //if(width >= max_window_size) break;

                    // extend down
                    if(extend_down)
                    {
                        tmp_vec.clear();
                        for(size_type idx=0;idx<min(width, matrix_size-x_idx);idx++)
                            tmp_vec.push_back(matrix.at(x_idx+idx).at(y_idx));
                        max_rc = *max_element(tmp_vec.cbegin(), tmp_vec.cend());
                        //if(max_rc>=extend_threshold and y_idx-x_idx>min_dist and y_idx>0 and y_idx>x_idx+width+min_dist)
                        if(max_rc>=extend_threshold and y_idx+height+min_dist<x_idx and y_idx>0 and height < max_window_size)
                        {
                            y_idx--;
                            height++;
                        }else{
                            extend_down = false;
                        }
                    }
                    //if(height >= max_window_size) break;

                    // extend right
                    if(extend_right)
                    {
                        tmp_vec.clear();
                        for(size_type idx=0;idx<min(height, matrix_size-y_idx);idx++)
                            tmp_vec.push_back(matrix.at( min(x_idx+width-1,matrix_size-1) ).at(y_idx+idx));
                        max_rc = *max_element(tmp_vec.cbegin(), tmp_vec.cend());
                        //if(max_rc>=extend_threshold and x_idx+width<matrix_size and y_idx>x_idx+width+min_dist)
                        if(max_rc>=extend_threshold and y_idx+height+min_dist<x_idx and width<max_window_size and x_idx+width<=matrix_size)
                        {
                            width++;
                        }else{
                            extend_right = false;
                        }
                    }
                    //if(width >= max_window_size) break;

                    // extend upper
                    if(extend_up)
                    {
                        tmp_vec.clear();
                        for(size_type idx=0;idx<min(width, matrix_size-x_idx);idx++)
                            tmp_vec.push_back(matrix.at(x_idx+idx).at( min(y_idx+height-1,matrix_size-1) ));
                        max_rc = *max_element(tmp_vec.cbegin(), tmp_vec.cend());
                        //if(max_rc>=extend_threshold and y_idx+height<matrix_size)
                        if(max_rc>=extend_threshold and y_idx+height+min_dist<x_idx and height < max_window_size)
                        {
                            height++;
                        }else{
                            extend_up = false;
                        }
                    }
                    //if(height >= max_window_size) break;
                }
              //  cout << "extend finish..." << endl;
                proper_x_upper = min(x_idx+width, matrix_size);
                proper_y_upper = min(y_idx+height, matrix_size);

                Region right(x_idx+1, proper_x_upper);
                Region left(y_idx+1, proper_y_upper);
                max_rc_record = locate_matrix_max_point<T>(matrix, max_pos, x_idx, proper_x_upper, y_idx, proper_y_upper, valid_pos_func);
                interact_regions.push_back( InterRegion(std::make_pair(left, right), max_rc_record, Point(max_pos.second, max_pos.first)) );

              //  cout << x_idx << "-" << x_idx+width << "\t" << y_idx << "-" << y_idx+height << endl;

                //std::cout << left << "\t" << right << std::endl;
                for(size_type idx=x_idx; idx<proper_x_upper; idx++)
                    for(size_type idy=y_idx; idy<proper_y_upper; idy++)
                        matrix.at(idx).at(idy) = T();

                /*
                for(auto x_iter=matrix.begin()+x_idx; x_iter<matrix.begin()+x_idx+min(width,matrix_size-x_idx); x_iter++)
                    for(auto y_iter=x_iter->begin()+y_idx; y_iter<x_iter->begin()+y_idx+min(height,matrix_size-y_idx); y_iter++)
                        *y_iter = T();
                */

                // again
                x_idx = 0;
               // cout << "end find..." << endl;
                break;
            }
        }
    }

    sort(interact_regions.begin(), interact_regions.end());

    //remove overlapped regions
    
    /*
    if(combine)
    {
        for(auto iter_1=interact_regions.begin(); iter_1!=interact_regions.end()-1; iter_1++)
        {
            for(auto iter_2=iter_1+1; iter_2!=interact_regions.end(); )
            {
                auto arm_overlap = overlap(*iter_1, *iter_2);
                if( arm_overlap.first > 0 and arm_overlap.second > 0 )
                {
                    auto left_start = min(iter_1->region.first.first, iter_2->region.first.first);
                    auto left_end = max(iter_1->region.first.second, iter_2->region.first.second);

                    auto right_start = min(iter_1->region.second.first, iter_2->region.second.first);
                    auto right_end = max(iter_1->region.second.second, iter_2->region.second.second);

                    iter_1->region.first.first = left_start;
                    iter_1->region.first.second = left_end;
                    iter_1->region.second.first = right_start;
                    iter_1->region.second.second = right_end;
                    
                    if( iter_2->max_cov >  iter_1->max_cov)
                    {
                        iter_1->max_cov = iter_2->max_cov;
                        iter_1->max_point = iter_2->max_point;
                    }
                    
                    iter_2 = interact_regions.erase(iter_2);

                }else{
                    iter_2++;
                }
            }
        }
    }
    */
}


template<typename T>
uLONG read_matrix(const string &matrix_file, Matrix<T> &matrix)
{
    matrix.clear();
    typename Matrix<T>::size_type matrix_size = 0;

    auto StringArray2DoubleArray = [](const StringArray &array)->DoubleArray
    {
         DoubleArray d_array;
         for(string num: array)
            d_array.push_back( stod(num) );
        return d_array;
    };

    ifstream IN(matrix_file, ifstream::in);
    if(not IN)
    {
        throw Bad_IO(matrix_file+" is unreadable", true);
    }

    string this_line;
    while(this_line.empty())
    {
        if(not getline(IN, this_line))
            throw Unexpected_Error(matrix_file+" is empty", true);
    }
    StringArray items;
    trim(this_line, '\t');
    split(this_line, items);
    matrix_size = items.size();
    matrix.push_back( StringArray2DoubleArray(items) );

    while(getline(IN, this_line))
    {
        trim(this_line, '\t');
        split(this_line, items);
        if(matrix_size != items.size())
            throw Unexpected_Error(matrix_file+" has different line length", true);
        matrix.push_back( StringArray2DoubleArray(items) );
    }

    if(matrix_size != matrix.size())
        throw Unexpected_Error(matrix_file+" column number - " + to_string(matrix_size) + " != row number - " + to_string(matrix.size()), true);

    return matrix_size;
}

template<typename T>
void fill_sym_matrix_with_aligned_dh(   Matrix<T> &matrix,
                                    const std::vector<Duplex_Hang> &dh_array, 
                                    const Multi_Align &alignment, 
                                    const string &chromosome_name, 
                                    const uINT target_size )
{
    using namespace std::placeholders;

    if(not alignment.has(chromosome_name))
        throw Unexpected_Error("Chromosome "+chromosome_name+" is not in alignment file");

    matrix.clear();
    init_matrix(matrix, target_size);
    auto coor_covert = std::bind(&Multi_Align::raw_coor_to_align_coor, _1, chromosome_name, _2);

    uLONG align_length = alignment.length();
    double step_size = 1.0 * align_length / target_size;
    for(const Duplex_Hang &dh: dh_array)
    {
        /*  Method 1 : sparse
        double left_center = 1.0*( coor_covert(alignment, dh.start_1-1) + coor_covert(alignment, dh.end_1-1) )/2 - 1;
        double right_center = 1.0*( coor_covert(alignment, dh.start_2-1) + coor_covert(alignment, dh.end_2-1) )/2 - 1;
        uLONG x_idx = floor(left_center / step_size);
        uLONG y_idx = floor(right_center / step_size);
        matrix.at(x_idx).at(y_idx) += 1;
        */

        /* Method 2 */
        for(uLONG left_start=dh.start_1; left_start <= dh.end_1; left_start++)
            for(uLONG right_start = dh.start_2; right_start <= dh.end_2; right_start++)
            {
                uLONG x_covert = coor_covert(alignment, left_start-1)-1;
                uLONG y_covert = coor_covert(alignment, right_start-1)-1;
                x_covert /= step_size;
                y_covert /= step_size;
                matrix.at(x_covert).at(y_covert) += 1;
            }
    }

    for(size_t x=0; x<target_size; x++)
        matrix.at(x).at(x) = T();
}


template<typename T>
void remove_paris_background(const Matrix<T> &raw_matrix, Matrix<T> &matrix, const uINT around)
{
    using size_type = typename Matrix<T>::size_type;
    auto matrix_size = raw_matrix.size();
    init_matrix(matrix, matrix_size);

    // build background
    vector<double> background;
    for(size_type idx=0; idx<matrix_size; idx++)
    {
        // arround
        size_type lower = idx>around ? idx-around : 0;
        size_type upper = idx+around<matrix_size ? idx+around : matrix_size-1;
        T around_total = T();
        T line_total = max( accumulate(raw_matrix.at(idx).cbegin(), raw_matrix.at(idx).cend(), T()), static_cast<T>(1) );
        for(size_type row=lower; row<=upper; row++)
        {
            around_total += accumulate(raw_matrix.at(row).cbegin(), raw_matrix.at(row).cend(), T());
        }
        background.push_back( max(1.0*around_total/line_total, 1.0) );
    }

    // calculate average
    vector<double> average;
    for(size_type idy=0; idy<matrix_size; idy++)
    {
        double Sum = 0;
        for(size_type idx=0; idx<matrix_size; idx++)
        {
            Sum += raw_matrix.at(idx).at(idy);
        }
        average.push_back( max(Sum/matrix_size, 1.0) );
    }

    // recalculate matrix
    for(size_type idx=0; idx<matrix_size; idx++)
        for(size_type idy=0; idy<matrix_size; idy++)
        {
            auto sqr_bg = sqrt(background.at(idx) * background.at(idy));
            //auto ave_ave = (average.at(idx)+average.at(idy))/2;
            auto sqr_ave = sqrt(average.at(idx) * average.at(idy)); //(average.at(idx)+average.at(idy))/2;
            //matrix.at(idx).at(idy) = 1.0 * raw_matrix.at(idx).at(idy) / sqr_bg  / ave_ave;
            matrix.at(idx).at(idy) = 1.0 * raw_matrix.at(idx).at(idy) / (sqr_bg  * sqr_ave);
        }
}

/* get a sub-matrix from raw matrix: coordination is 0-based, and left-close; right-open */
template<typename T>
void sub_matrix(const Matrix<T> &raw_matrix, Rect<T> &sub_rect, uLONG x_start, uLONG x_end, uLONG y_start, uLONG y_end)
{
    if(x_start==0 and y_start==0 and x_end==-1UL and y_end==-1UL)
    {
        sub_rect = raw_matrix;
        return;
    }

    sub_rect.clear();
    const uLONG raw_size = raw_matrix.size();
    if(x_end>raw_size) x_end = raw_size;
    if(y_end>raw_size) y_end = raw_size;
    if(x_start >= x_end) throw Unexpected_Error("sub_matrix: range is illegal");
    if(y_start >= y_end) throw Unexpected_Error("sub_matrix: range is illegal");

    const uLONG new_x_width = x_end - x_start;
    const uLONG new_y_width = y_end - y_start;

    init_rect(sub_rect, new_x_width, new_y_width);
    for(uLONG x=x_start; x<x_end; x++)
        for(uLONG y=y_start; y<y_end; y++)
            sub_rect[x-x_start][y-y_start] = raw_matrix[x][y];
}



}
#endif


