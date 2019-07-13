
#ifndef PAN_TYPE_H
#define PAN_TYPE_H

#include <unordered_map>
#include <vector>
#include <iostream>
#include <string>
#include <functional>
#include <memory>
#include <cstring>

#define PsBL_LIB_VERSION "PsBL 1.1.0(2018-11-2)"

using std::unordered_map;
using std::vector;
using std::string;
using std::pair;
using std::ostream;
using std::shared_ptr;

namespace pan{

using uINT = unsigned int;
using uLONG = unsigned long;
using uLONGLONG = unsigned long long;

using StringArray = vector<string>;
using StringMatrix = vector<vector<string>>;
using IntArray = vector<int>;
using uIntArray = vector<unsigned int>;
using uLONGArray = vector<unsigned long>;
using uLONGLONGArray = vector<unsigned long long>;
using FloatArray = vector<float>;
using DoubleArray = vector<double>;

using MapStringString = unordered_map<string, string>;
using MapStringuLONG = unordered_map<string, uLONG>;
using MapStringDouble = unordered_map<string, double>;


//using Region = pair<uLONG, uLONG>;
//using Point = pair<uLONG, uLONG>;
//using PointF = pair<double, double>;

struct Region
{
    uLONG first;
    uLONG second;

    //bool is_null(){ return first<=second?false:true; };
    Region(pair<uLONG, uLONG>init_pair):first(init_pair.first), second(init_pair.second){}
    Region(uLONG first, uLONG second):first(first), second(second){}
    Region()=default;
};

inline bool operator==(const Region &r_1, const Region &r_2){ return (r_1.first==r_2.first and r_1.second==r_2.second) ? true : false; }
inline bool operator!=(const Region &r_1, const Region &r_2){ return not(r_1==r_2); }
inline ostream& operator<<(ostream& OUT, const Region& r){ OUT << r.first << "-" << r.second; return OUT; }
bool operator<(const Region &r_1, const Region &r_2);
uLONG overlap(const Region &r_1, const Region &r_2);

vector<Region> operator+(const Region &r_1, const Region &r_2);
vector<Region> operator-(const Region &r_1, const Region &r_2);



struct Point
{
    uLONG first;
    uLONG second;

    Point(pair<uLONG, uLONG>init_pair):first(init_pair.first), second(init_pair.second){}
    Point(uLONG first, uLONG second):first(first), second(second){}
    Point()=default;
};

inline bool operator==(const Point &p_1, const Point &p_2){ return (p_1.first==p_2.first and p_1.second==p_2.second) ? true : false; }
inline bool operator!=(const Point &p_1, const Point &p_2){ return not(p_1==p_2); }
inline Point operator+(const Point &p_1, const Point &p_2){ return Point( p_1.first+p_2.first, p_1.second+p_2.second ); }
inline Point operator-(const Point &p_1, const Point &p_2){ return Point( p_1.first-p_2.first, p_1.second-p_2.second ); }
inline ostream& operator<<(ostream& OUT, const Point& p){ OUT << "(" << p.first << "," << p.second << ")"; return OUT; }

struct PointF
{
    double first;
    double second;

    PointF(pair<double, double>init_pair):first(init_pair.first), second(init_pair.second){}
    PointF(uLONG first, uLONG second):first(first), second(second){}
    PointF()=default;
};

inline bool operator==(const PointF &p_1, const PointF &p_2){ return (p_1.first==p_2.first and p_1.second==p_2.second) ? true : false; }
inline bool operator!=(const PointF &p_1, const PointF &p_2){ return not(p_1==p_2); }
inline PointF operator+(const PointF &p_1, const PointF &p_2){ return PointF( p_1.first+p_2.first, p_1.second+p_2.second ); }
inline PointF operator-(const PointF &p_1, const PointF &p_2){ return PointF( p_1.first-p_2.first, p_1.second-p_2.second ); }
inline ostream& operator<<(ostream& OUT, const PointF& p){ OUT << "(" << p.first << ", " << p.second << ")"; return OUT; }

using RegionArray = vector<Region>;
bool sorted(const RegionArray &ra);
bool sorted_no_overlap(const RegionArray &ra);

//typedef template unordered_map<string, T> MapStringT<T>;
template<typename T> 
using MapStringT = unordered_map<string, T>;


// =========== MATRIX ===========
template<typename T>
using Matrix = vector<vector<T>>;

template<typename F_T, typename S_T>
ostream& operator<<(ostream& OUT, const pair<F_T, S_T> &cur_pair);
template<typename T>
ostream& operator<<(ostream& OUT, const vector<T> &cur_vec);
template<typename T>
ostream& operator<<(ostream& OUT, const vector<vector<T>> &cur_matrix);
template<typename K_T, typename V_T>
ostream& operator<<(ostream& OUT, const unordered_map<K_T, V_T> &cur_map);



/* ========================

         Pointer

   ======================== */
template<typename T>
using sp = shared_ptr<T>;










template<typename T>
ostream& operator<<(ostream& OUT, const vector<T> &cur_vec)
{
    for(auto iter=cur_vec.cbegin(); iter!=cur_vec.cend(); iter++)
    {
        OUT << *iter;
        if(iter!=cur_vec.cend()-1)
            OUT << "\t";
    }
    return OUT;
}

template<typename F_T, typename S_T>
ostream& operator<<(ostream& OUT, const pair<F_T, S_T> &cur_pair)
{
    OUT << cur_pair.first << "-" << cur_pair.second;
    return OUT;
}

template<typename T>
ostream& operator<<(ostream& OUT, const vector<vector<T>> &cur_matrix)
{
    for(auto iter_1=cur_matrix.cbegin(); iter_1!=cur_matrix.cend(); iter_1++)
    {
        for(auto iter_2=iter_1->cbegin(); iter_2!=iter_1->cend(); iter_2++)
        {
            OUT << *iter_2;
            if(iter_2!=iter_1->cend()-1)
                OUT << "\t";
        }
        OUT << "\n";
    }

    return OUT;
}

template<typename K_T, typename V_T>
ostream& operator<<(ostream& OUT, const unordered_map<K_T, V_T> &cur_map)
{
    for(auto iter=cur_map.cbegin(); iter!=cur_map.cend(); iter++)
    {
        OUT << iter->first << "\t" << iter->second << "\n";
    }

    return OUT;
}

template<typename T>
void init_matrix(Matrix<T> &matrix, typename Matrix<T>::size_type size, T init_value=T())
{
    matrix.clear();
    vector<T> row;

    T t(init_value);
    row.resize(size, t);
    matrix.resize(size, row);
    for_each(matrix.begin(), matrix.end(), [](vector<T> &cur_row)->void{ cur_row.shrink_to_fit(); });
}

template<typename T>
using Rect = vector<vector<T>>;

template<typename T>
void init_rect(Rect<T> &matrix, typename Rect<T>::size_type row_num, typename Rect<T>::size_type col_num, T init_value=T())
{
    matrix.clear();
    vector<T> row;

    T t(init_value);
    row.resize(col_num, t);
    matrix.resize(row_num, row);
    for_each(matrix.begin(), matrix.end(), [](vector<T> &cur_row)->void{ cur_row.shrink_to_fit(); });
}




}
#endif