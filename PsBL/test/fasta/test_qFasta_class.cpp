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

    qFasta seq(argv[1]);
    //seq.load_fasta_file("mm10_transcriptome.fa");

    cout << "test get_chr_lens..." << endl;
    auto chr_lens = seq.get_chr_lens();
    for(auto it=chr_lens.cbegin(); it!=chr_lens.cend(); it++)
    {
        cout << it->first << "\t" << it->second << "\n";
    }
    cout << endl << endl;

    cout << "test get_chr_seq..." << endl;
    cout << "ENSMUST00000174899.2 sequence" << endl;
    cout << seq.get_chr_seq("ENSMUST00000174899.2") << endl << endl;
    cout << "Flatten sequence ENSMUST00000215494.1" << endl;
    cout << flat_seq(seq.get_chr_seq("ENSMUST00000215494.1")) << endl;

    cout << "test get_chr_subbseq..." << endl;
    cout << "ENSMUST00000146358.7 39, 76 +" << endl;
    cout << seq.get_chr_subbseq("ENSMUST00000146358.7", 39, 76, POSITIVE) << endl;
    cout << "ENSMUST00000146358.7 39, 76 -" << endl;
    cout << seq.get_chr_subbseq("ENSMUST00000146358.7", 39, 76, NEGATIVE) << endl << endl;

    return 0;
}
