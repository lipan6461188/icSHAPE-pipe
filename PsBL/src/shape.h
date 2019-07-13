

#ifndef ICSHAPE_H
#define ICSHAPE_H

#include "pan_type.h"

#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <unordered_map>
#include <cstring>
#include <memory>
#include <iostream>

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
using std::shared_ptr;
using std::cout;
using std::endl;


#define NULL_Score -1

enum Shape_NULL_Type
{
    minus_999, minus_555, minus_1, zero, plus_1, plus_555, plus_999, remove
};

struct Chr_Shape
{
    string chr_id;
    double rpkm = 0.0;
    FloatArray shape_score;

    friend istream &operator >>(istream &, Chr_Shape &);
};

class icSHAPE
{
public:

    using sp_Chr_Shape = shared_ptr<Chr_Shape>;

    icSHAPE(const string &shape_file_name){ read_shape(shape_file_name); }
    icSHAPE(const icSHAPE &);
    icSHAPE& operator=(const icSHAPE &);
    icSHAPE(icSHAPE &&);
    icSHAPE& operator=(icSHAPE &&);

    uINT size() const { return chr_id_list.size(); }
    const StringArray& keys() const { return chr_id_list; }
    // not find: out_of_range
    shared_ptr<Chr_Shape> get_shape(const string &transID) const;
    bool has(const string &transID)const { return shape_list.find(transID) != shape_list.end(); }

    // output icSHAPE score to ostream
    void to_bedGraph(ostream &, Shape_NULL_Type flag=remove) const; //all transcript
    void to_bedGraph(ostream &, const string &, Shape_NULL_Type flag=remove) const; //single transcript
    void to_bedGraph(string &, Shape_NULL_Type flag=remove) const; //all transcript
    void to_bedGraph(string &, const string &, Shape_NULL_Type flag=remove) const; //single transcript

    MapStringT<sp_Chr_Shape>::const_iterator cbegin() const { return shape_list.cbegin(); }
    MapStringT<sp_Chr_Shape>::const_iterator cend() const { return shape_list.cend(); }
    MapStringT<sp_Chr_Shape>::iterator begin(){ return shape_list.begin(); }
    MapStringT<sp_Chr_Shape>::iterator end(){ return shape_list.end(); }

private:
    uINT capacity = 0;

    MapStringT<sp_Chr_Shape> shape_list;
    StringArray chr_id_list;

    void read_shape( const string &shape_file_name );
    shared_ptr<string> to_bedGraph(const string &, Shape_NULL_Type flag=remove) const;
};

// Exception
class Invalid_Shape_Line: public runtime_error 
{ 
public:
    explicit Invalid_Shape_Line(const std::string& what_arg, bool fatal=false): runtime_error(what_arg), fatal(fatal) {};
    explicit Invalid_Shape_Line(const char* what_arg, bool fatal=false): runtime_error(what_arg), fatal(fatal) {};

    bool fatal;

};

class End_Of_File: public runtime_error
{
public:
    End_Of_File(): runtime_error("End of File") {}
};

}
#endif // ICSHAPE_H















