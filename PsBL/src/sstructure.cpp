
#include "sstructure.h"
#include "string_split.h"

#include <sstream>
#include <iostream>
#include <iomanip>

using namespace std;

namespace pan{

void SStructure::parse_dot_bracket(const string &dot_bracket)
{
    uINT index = 1;
    stack<uINT> base_bucket;
    for(auto iter=dot_bracket.cbegin(); iter!=dot_bracket.cend(); iter++)
    {
        if(*iter == '(' or *iter == '{' or *iter == '<' or *iter == '[')
            base_bucket.push( index );
        else if(*iter == ')' or *iter == '}' or *iter == '>' or *iter == ']')
        {
            Base_Pair bp(base_bucket.top(), index);
            base_bucket.pop();
            base_pairs.push_back( bp );
        }else{
            if(*iter != '.' and *iter != '_' and *iter != '-' and *iter != '*' and *iter != '+' and *iter != '=')
                throw runtime_error("Bad dot-bracket");
        }
        index++;
    }
    if(base_bucket.size() != 0)
    {
        throw runtime_error("Bad dot-bracket");
    }

    length = dot_bracket.size();
    sort(base_pairs.begin(), base_pairs.end(), [](const Base_Pair &bp_1, const Base_Pair &bp_2){return bp_1.left < bp_2.left;});
}

string SStructure::to_dot_bracket() const
{
    string dot_bracket(length, '.');
    for(auto iter=base_pairs.cbegin(); iter!=base_pairs.cend(); iter++)
    {
        uLONG left_base_index = iter->left - 1;
        uLONG right_base_index = iter->right - 1;

        dot_bracket[left_base_index] = '(';
        dot_bracket[right_base_index] = ')';
    }

    return dot_bracket;
}

void  SStructure::to_ct(ostream &OUT) const
{
    OUT << setw(6) << length << "\t" << title << "\n";
    
    unordered_map<uLONG, uLONG> bp_map;
    for(auto iter=base_pairs.cbegin(); iter!=base_pairs.cend(); iter++)
        bp_map[iter->left] = iter->right;

    for(uLONG index=0; index<length; index++)
    {
        char base('.');
        if(not sequence.empty())
            base = sequence[index];
        if(bp_map.find(index+1) != bp_map.end())
            OUT << setw(6) << index + 1 << ' ' << base << setw(6) 
                << index << setw(6) 
                << ((index+2>length)?0:(index+2)) << setw(6) 
                << bp_map[index+1] << setw(6)
                << index+1 << "\n";
        else
            OUT << setw(6) << index + 1 << ' ' << base << setw(6) 
                << index << setw(6) 
                << ((index+2>length)?0:(index+2)) << setw(6) 
                << 0 << setw(6)
                << index+1 << "\n";
    }
}

void SStructure::to_ct(const string &file_name) const
{
    ofstream OUT(file_name, ofstream::out);
    if(not OUT)
    {
        throw runtime_error(file_name+" cannot writable");
    }
    to_ct(OUT);
    OUT.close();
}

vector<SStructure> SStructure::read_ct_file(const string &file_name)
{
    vector<SStructure> structures;
    ifstream IN(file_name, ifstream::in);
    if(not IN)
    {
        throw runtime_error(file_name+" cannot readable");
    }

    string cur_line;

    BpArray bp_array;
    string title;
    uLONG length;
    string sequence;
    // read head
    while(getline(IN, cur_line))
    {
        trim(cur_line);
        if(cur_line.empty())
            continue;
        StringArray items;
        split(cur_line, items);
        length = stoul(items[0]);
        title = items.back();
        break;
    }

    while(getline(IN, cur_line))
    {
        trim(cur_line);
        if(cur_line.empty())
            continue;
        StringArray items;
        split(cur_line, items);
        if(items.size() != 6)
            throw runtime_error("Bad Ct Line");
        uLONG index = stoul(items[0]);
        sequence.push_back(items[1][0]);
        uLONG pair_base = stoul(items[4]);
        if(pair_base != 0)
        {
            Base_Pair bp(index, pair_base);
            bp_array.push_back(bp);
        }
        if(index==length)
        {
            SStructure structure(bp_array, length);
            structure.setTitle(title);
            structure.setSequence(sequence);
            structures.push_back(structure);

            bp_array.clear();
            length = 0;
            sequence.clear();
            title.clear();

            // read next title
            while(getline(IN, cur_line))
            {
                trim(cur_line);
                if(cur_line.empty())
                    continue;
                StringArray items;
                split(cur_line, items);
                length = stoul(items[0]);
                title = items.back();
                break;
            }
        }
    }

    return structures;
}










}