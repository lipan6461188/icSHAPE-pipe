
#include "../../src/sam.h"
#include "../../src/htslib.h"
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <fstream>

using namespace std;
using namespace pan;

int main(int argc, char *argv[])
{

    if(argc < 2)
    {
        cerr << "Usage: test_read_bamorsam input.bam/input.sam" << endl;
        return 0;
    }

    if( endswith(argv[1], ".bam") )
    {

        BGZF* fn_hd = bgzf_open(argv[1], "r");
        if(not fn_hd)
        {
            cerr << "open failed" << endl;
            return -1;
        }

        // read the bam head
        bam_hdr_t *hdr = bam_hdr_read(fn_hd);

        vector<Sam_Record> records;     // reads
        Sam_Record *p = nullptr;        // cache
        while(read_a_read_record(fn_hd, hdr, records, p))
        {
            cout << "=======" << endl;
            cout << records;
        }

        bgzf_close(fn_hd);

    }else if( endswith(argv[1], ".sam") )
    {
        ifstream IN;
        IN.open(argv[1]);

        Sam_Head sam_head;
        read_sam_head(IN, sam_head);
        
        vector<Sam_Record> records;
        while( read_a_read_record(IN, records) )
        {
            cout << "=======" << endl;
            cout << records;
        }

        IN.close();
    }
}



