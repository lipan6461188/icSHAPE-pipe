#include <iostream>
#include "fold.h"

using namespace std;

namespace pan{

/* number starts from 1, returned pairs is 1-base coorinates */
vector< pair<uLONG, uLONG> > get_pairs(RNA &ct, const uLONG number)
{
    if(ct.GetStructureNumber() < static_cast<int>(number))
    {
        //cout << "\nGetStructureNumber: " << ct.GetStructureNumber() << "\tnumber:" << number << endl;
        throw Unexpected_Error( "Structure number out of range: " + to_string(ct.GetStructureNumber()) + "<" + to_string(number) );
    }

    int seq_len = ct.GetSequenceLength();
    vector< pair<uLONG, uLONG> > pairs;


    for(int j=1; j<=seq_len; j++)
    {
        int paired_index = ct.GetPair(j, number);

        if(j < paired_index)
            pairs.push_back(make_pair(j, paired_index));
    }


    sort(pairs.begin(), pairs.end(), [](const pair<uLONG, uLONG> &pair_1, const pair<uLONG, uLONG> &pair_2){ return pair_1.first < pair_2.first; });
    return pairs;
}

sp<StringArray> get_dot_bracket(RNA &ct)
{
    int number = ct.GetStructureNumber();
    int seq_len = ct.GetSequenceLength();

    sp<StringArray> structures = make_shared<StringArray>();

    //char GetNucleotide(const int i);
    //int GetSequenceLength() const;
    //const char* GetSequence();

    for(int i=1; i<=number; i++)
    {
        string cur_structure;
        for(int j=1; j<=seq_len; j++)
        {
            int paired_index = ct.GetPair(j, i);

            if(j < paired_index)
                cur_structure += '(';
            else if(paired_index == 0)
                cur_structure += '.';
            else
                cur_structure += ')';
        }
        structures->push_back(cur_structure);
    }

    return structures;
}

sp<StringArray> fold_two_seq(const string &seq_1, const string &seq_2, const bool forbidIntramolecular)
{
    const double temperature = 310.15;
    const bool isRNA = true;

    HybridRNA strand(seq_1.c_str(), seq_2.c_str(), isRNA);
    strand.SetForbidIntramolecular( forbidIntramolecular );
    strand.SetTemperature( temperature );

    // Initialize the maximum percent energy difference.
    const double percent = 50.0;
    // Initialize the maximum number of structures.
    const int maxStructures = 20;
    // Initialize the folding window size.
    const int windowSize = 0;
    //Initialize the maximum internal bulge loop size.
    const int maxLoop = 30;

    strand.FoldBimolecular( percent, maxStructures, windowSize, "", maxLoop );

    sp<StringArray> structures = get_dot_bracket(strand);
    string sequence = strand.GetSequence();

    structures->insert( structures->begin(), sequence.substr(0, strand.GetSequenceLength()) );

    return structures;
}











}