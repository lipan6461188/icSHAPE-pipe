 

#include "shape.h"
#include "string_split.h"

using namespace std;

namespace pan{

istream &operator >>(istream &IN, Chr_Shape &shape)
{
    string this_line, this_word;
    if(getline(IN, this_line))
    {
        trim(this_line);
        if(this_line.empty())
        {
            bool fatal = false;
            throw Invalid_Shape_Line("Empty Line", fatal);
        }else if(this_line[0] == '#')
        {
            bool fatal = false;
            throw Invalid_Shape_Line("Comments Line", fatal);
        }

        uLONG length;

        istringstream input_stream(this_line);
        input_stream >> shape.chr_id >> length >> shape.rpkm;
        while(input_stream >> this_word)
        {
            if(this_word == "NULL")
                shape.shape_score.push_back(NULL_Score);
            else
                shape.shape_score.push_back( stof(this_word) );
        }
        if(shape.shape_score.size() != length or length == 0)
        {
            bool fatal = true;
            throw Invalid_Shape_Line(shape.chr_id+" Bad Shape Line -- different length or length of 0", fatal);
        }
    }else{
        throw End_Of_File();
    }
    return IN;
}


void icSHAPE::read_shape(const string &shape_file_name)
{
    ifstream IN(shape_file_name, ifstream::in);
    if(IN)
    {
        while(1)
        {
            shared_ptr<Chr_Shape> new_shape(new Chr_Shape);
            try{
                IN >> *new_shape;
            }catch(Invalid_Shape_Line error)
            {
                if(error.fatal)
                {
                    throw error;
                }
                if(IN.fail())
                    break;
                continue;
            }catch(End_Of_File)
            {
                break;
            }
            chr_id_list.push_back(new_shape->chr_id);
            shape_list[new_shape->chr_id] = new_shape;
        }
        sort(chr_id_list.begin(), chr_id_list.end());
    }else{
        throw runtime_error(shape_file_name+" file is unreadable");
    }
    IN.close();
}

shared_ptr<Chr_Shape> icSHAPE::get_shape(const string &chr_id) const
{
    // not find: out_of_range

    return shape_list.at(chr_id);
}

icSHAPE::icSHAPE(const icSHAPE &other):
    capacity(other.capacity),
    shape_list(other.shape_list),
    chr_id_list(other.chr_id_list)
{
    // null
}

icSHAPE& icSHAPE::operator=(const icSHAPE &other)
{
    //capacity = other.capacity;
    shape_list = other.shape_list;
    chr_id_list = other.chr_id_list;

    return *this;
}

icSHAPE::icSHAPE(icSHAPE &&other):
    //capacity(std::move(other.capacity)),
    shape_list(std::move(other.shape_list)),
    chr_id_list(std::move(other.chr_id_list))
{
    other.shape_list.clear();
    other.chr_id_list.clear();
}

icSHAPE& icSHAPE::operator=(icSHAPE &&other)
{
    //capacity = std::move(other.capacity);
    shape_list = std::move(other.shape_list);
    chr_id_list = std::move(other.chr_id_list);

    other.shape_list.clear();
    other.chr_id_list.clear();

    return *this;
}

shared_ptr<string> icSHAPE::to_bedGraph(const string &chr_id, Shape_NULL_Type flag) const
{
    double null_type = 0;
    switch(flag){
        case minus_999:
            null_type = -999;
            break;
        case minus_555:
            null_type = -555;
            break;
        case minus_1:
            null_type = -1;
            break;
        case zero:
            null_type = 0;
            break;
        case plus_1:
            null_type = 1;
            break;
        case plus_555:
            null_type = 555;
            break;
        case plus_999:
            null_type = 999;
            break;
        default:
            break;
    }

    ostringstream out_sstream;
    out_sstream.precision(3);
    out_sstream << std::fixed << showpoint;
    unsigned start_loc = 0;
    for(float shape_score: shape_list.at(chr_id)->shape_score)
    {
        if(shape_score != NULL_Score)
            out_sstream << chr_id << "\t" << start_loc << "\t" << start_loc+1 << "\t" << shape_score << "\n";
        else if( flag == remove )
        {
            // skip, do_nothing
        }else{
            out_sstream << chr_id << "\t" << start_loc << "\t" << start_loc+1 << "\t" << null_type << "\n";
        }
        ++start_loc;
    }

    return make_shared<string>(out_sstream.str());
}

// output icSHAPE score to ostream
void icSHAPE::to_bedGraph(ostream &OUT, Shape_NULL_Type flag) const //all transcript
{
    for(const string &chr_id: chr_id_list)
        OUT << *to_bedGraph(chr_id, flag);
}

void icSHAPE::to_bedGraph(ostream &OUT, const string &chr_id, Shape_NULL_Type flag) const //single transcript
{
    OUT << *to_bedGraph(chr_id, flag);
}

void icSHAPE::to_bedGraph(string & bigString, Shape_NULL_Type flag) const //all transcript
{
    for(const string &chr_id: chr_id_list)
        bigString += *to_bedGraph(chr_id, flag);
}

void icSHAPE::to_bedGraph(string & bigString, const string &chr_id, Shape_NULL_Type flag) const //single transcript
{
    bigString += *to_bedGraph(chr_id, flag);
}


}
