#include "pan_type.h"
#include <cmath>
#include <algorithm>

using namespace std;

namespace pan{


bool operator<(const Region &r_1, const Region &r_2)
{ 
    if(r_1.first<r_2.first) 
        return true; 
    else if(r_1.first==r_2.first and r_1.second < r_2.second)
        return true;
    else
        return false;
}

vector<Region> operator+(const Region &r_1, const Region &r_2)
{
    vector<Region> regions;
    if(r_1.first <= r_2.second and r_2.first <= r_1.second)
    {
        regions.push_back( Region(min(r_1.first, r_2.first), max(r_1.second, r_2.second)) );
    }else{
        regions.push_back( r_1 );
        regions.push_back( r_2 );
        sort(regions.begin(), regions.end());
    }
    return regions;
}

vector<Region> operator-(const Region &r_1, const Region &r_2)
{
    vector<Region> regions;
    if(r_1.first <= r_2.second and r_2.first <= r_1.second)
    {
        if(r_1.first < r_2.first)
            regions.push_back(Region(r_1.first, r_2.first-1));
        if(r_1.second > r_2.second)
            regions.push_back(Region(r_2.second+1, r_1.second));
    }
    return regions;
}

uLONG overlap(const Region &r_1, const Region &r_2)
{
    if(r_1.first <= r_2.second and r_2.first <= r_1.second)
    {
        return min(r_1.second, r_2.second) - max(r_1.first, r_2.first) + 1;
    }else{
        return 0;
    }
}

bool sorted(const RegionArray &ra)
{
    if(ra.empty())
        return true;

    for(auto iter=ra.cbegin(); iter!=ra.cend()-1; iter++)
        if(iter->first > (iter+1)->first)
        {
            //cout << iter->first << "\t" << (iter+1)->first << endl;
            return false;
        }

    return true;
}

bool sorted_no_overlap(const RegionArray &ra)
{
    if(ra.empty())
        return true;

    if(sorted(ra))
    {
        for(auto iter=ra.cbegin(); iter != ra.cend()-1; iter++)
            if(iter->second >= (iter+1)->first)
            {
                //cout << iter->second << "\t" << (iter+1)->first << endl;
                return false;
            }
        return true;
    }else{
        return false;
    }
}

}