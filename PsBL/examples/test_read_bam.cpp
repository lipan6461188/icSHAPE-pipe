
#include <sam.h>
#include <htslib.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>

using namespace std;
using namespace pan;

int main(int argc, char *argv[])
{
    BGZF* fn_hd = bgzf_open(argv[1], "r");

    if(not fn_hd)
    {
        cerr << "open failed" << endl;
        return -1;
    }

    // read the sam head
    bam_hdr_t *hdr = bam_hdr_read(fn_hd);

    vector<Sam_Record> records;
    
    //Sam_Record records;
    Sam_Record *p = nullptr;
    while(read_a_read_record(fn_hd, hdr, records, p))
    {
        //cout << p << endl;
        cout << "=======" << endl;
        cout << records;
    }

    bgzf_close(fn_hd);
}



