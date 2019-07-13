#include "string_split.h"

#include <sstream>
#include <ctype.h>

using namespace std;

namespace pan{


string::size_type findany(const string &String, const char *find_char, string::size_type begin)
{
    for(;begin<String.size();++begin)
        for(const char *c=find_char; *c != '\0'; ++c)
            if(*c == String.at(begin))
                return begin;
    return string::npos;
}

void split(const string &String, const char *seperators, StringArray &Vec)
{
    using size_t = string::size_type;
    size_t last = 0;
    Vec.clear();
    size_t begin = findany(String, seperators, 0);  //String.find(c, 0);
    while(begin != string::npos)
    {
        Vec.emplace_back(String.substr(last, begin-last));
        last = begin + 1;
        begin = findany(String, seperators, begin+1); //String.find(c, begin+1);
    }
    Vec.emplace_back(String.substr(last, begin-last));
}

// split a string with a char
StringArray split(const string &String, const char &c)
{
    StringArray Vec;
    size_t last = 0;
    size_t begin = String.find(c, 0);
    while(begin != string::npos)
    {
        Vec.emplace_back(String.substr(last, begin-last));
        last = begin + 1;
        begin = String.find(c, begin+1);
    }
    Vec.emplace_back(String.substr(last, begin-last));
    return Vec;
}


// a higher efficient method to split s a string
void split(const string &String, const char &c, StringArray &Vec)
{
    string::size_type last = 0;
    Vec.clear();
    string::size_type begin = String.find(c, 0);
    while(begin != string::npos)
    {
        Vec.emplace_back(String.substr(last, begin-last));
        last = begin + 1;
        begin = String.find(c, begin+1);
    }
    Vec.emplace_back(String.substr(last, begin-last));
}

// split a string with blank spaces
void split(const string &String, StringArray &Vec)
{
    Vec.clear();
    istringstream ss(String);
    string cur_item;
    while( ss >> cur_item ) Vec.push_back(cur_item);
}

// trimming a string
void trim(string &String, const char &c)
{

    /*
    if(String.empty())
        return;

    string::size_type start_point = 0;
    for(;start_point<String.size();start_point++)
        if(String.at(start_point) != c)
            break;

    string::size_type end_point = String.size() - 1;
    for(;end_point>=start_point;end_point--)
        if(String.at(end_point) != c)
            break;

    if(end_point < String.size())
        String.erase(end_point+1, String.size());

    if(start_point > 0)
        String.erase(0, start_point);

    */

    trim_right(String, c);
    trim_left(String, c);

}

void trim(string &String)
{
    trim_right(String);
    trim_left(String);
}

void trim_left(string &String)
{
    if(String.empty())
        return;

    if( not isspace(String.front()) )
        return;

    string::size_type start_point = 0;
    for(;start_point<String.size()-1;start_point++)
        if( not isspace(String.at(start_point+1)) )
            break;

    //if(start_point > 0)
    String.erase(0, start_point+1);
}

void trim_right(string &String)
{
    if(String.empty())
        return;

    string::size_type end_point = String.size();
    for(;end_point!=0;end_point--)
        if( not isspace(String.at(end_point-1)) )
            break;

    //if(end_point < String.size()-1)
    String.erase(end_point, String.size());
}

void trim_left(string &String, const char &c)
{
    if(String.empty())
        return;

    if( String.front() != c )
        return;

    string::size_type start_point = 0;
    for(;start_point<String.size()-1;start_point++)
        if( String.at(start_point+1) != c )
            break;

    //if(start_point > 0)
    String.erase(0, start_point+1);
}

void trim_right(string &String, const char &c)
{
    if(String.empty())
        return;

    string::size_type end_point = String.size();
    for(;end_point!=0;end_point--)
        if( String.at(end_point-1) != c )
            break;

    //if(end_point < String.size()-1)
    String.erase(end_point, String.size());
}


bool startswith(const char *ref_string, const char *head_string)
{
    size_t ref_len = strlen(ref_string);
    size_t query_len = strlen(head_string);

    if(ref_len < query_len)
        return false;
    for(size_t i=0; i<query_len; i++)
        if(ref_string[i] != head_string[i])
            return false;

    return true;
}

bool endswith(const char *ref_string, const char *tail_string)
{
    size_t ref_len = strlen(ref_string);
    size_t query_len = strlen(tail_string);

    if(ref_len < query_len)
        return false;
    int query_i = 0;
    for(size_t i=ref_len-query_len; i<ref_len; i++, query_i++)
        if(ref_string[i] != tail_string[query_i])
            return false;
    return true;
}

string reverse_string(const string &input_string)
{
    string new_string;
    for(auto iter=input_string.crbegin(); iter!=input_string.crend(); iter++)
        new_string += *iter;
    return new_string;
}


/*

void trim2(string &String, const char &c)
{
    if(String.empty())
        return;
    string::size_type start_point = 0;
    for(;start_point<String.size();start_point++)
        if(String.at(start_point) != c)
            break;

    string::size_type end_point = String.size() - 1;
    for(;end_point>=start_point;end_point--)
        if(String.at(end_point) != c)
            break;

    String = String.substr(start_point,end_point-start_point+1);
}

*/


}