

#include "../../src/fasta.h"
#include <iostream>

using namespace std;
using namespace pan;

int main(int argc, char *argv[])
{

    if(argc < 2)
    {
        cerr << "Usage: test_Fasta_class input_fasta.fa" << endl;
        return 0;
    }

    // Condition 1
    Fasta seq1(argv[1]);

    cout << "Chr nums: " << seq1.get_chr_num() << endl;
    cout << "test get_chr_lens..." << endl;
    auto chr_lens1 = seq1.get_chr_lens();
    for(auto it=chr_lens1.cbegin(); it!=chr_lens1.cend(); it++)
    {
        cout << it->first << "\t" << it->second << "\n";
    }
    cout << endl << endl;

    // Condition 2
    Fasta seq2;
    seq2.add_fasta(argv[1]);

    cout << "test get_chr_lens..." << endl;
    auto chr_lens2 = seq2.get_chr_lens();
    for(auto it=chr_lens2.cbegin(); it!=chr_lens2.cend(); it++)
    {
        cout << it->first << "\t" << it->second << "\n";
    }
    cout << endl << endl;

    cout << "test get_chr_seq..." << endl;
    cout << "ENSMUST00000174899.2 sequence" << endl;
    cout << seq2.get_chr_seq("ENSMUST00000174899.2") << endl << endl;
    cout << "Flatten sequence ENSMUST00000215494.1" << endl;
    cout << flat_seq(seq2.get_chr_seq("ENSMUST00000215494.1")) << endl;

    cout << "test get_chr_subbseq..." << endl;
    cout << "ENSMUST00000146358.7 39, 76 +" << endl;
    cout << seq2.get_chr_subbseq("ENSMUST00000146358.7", 39, 76, POSITIVE) << endl;
    cout << "ENSMUST00000146358.7 39, 76 -" << endl;
    cout << seq2.get_chr_subbseq("ENSMUST00000146358.7", 39, 76, NEGATIVE) << endl << endl;

    return 0;
}
