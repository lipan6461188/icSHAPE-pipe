
#include "../../src/sam.h"
#include "../../src/htslib.h"
#include <htslib/sam.h>
#include <htslib/bgzf.h>

using namespace std;
using namespace pan;

int main(int argc, char *argv[])
{

    if(argc < 3)
    {
        cerr << "Usage: test_htslib input_bam.bam num" << endl;
        return 0;
    }

    BGZF* fn_hd = bgzf_open(argv[1], "r");

    if(not fn_hd)
    {
        cerr << "open failed" << endl;
        return -1;
    }

    // read the sam head
    bam_hdr_t *hdr = bam_hdr_read(fn_hd);

    for(int i=0; i<stoi(argv[2]); i++)
    {
        bam1_t *record = bam_init1();
        int ret = bam_read1(fn_hd, record);
        if(ret < 0)
        {
            cerr << "end of file" << endl;
            break;
        }
        
        Sam_Record read_record;
        read_record.read_id = getBamQName(record);
        read_record.flag = getBamFlag(record);
        if(read_record.flag & 4)
        {
           // unmap
            read_record.chr_id = "*";
            read_record.pos = 0;
            read_record.map_quanlity = 0;
        }else{
            read_record.chr_id = getBamRef(record, hdr);
            read_record.pos = getBamRefPos(record);
            read_record.map_quanlity = getBamMapQuanlity(record);
            read_record.cigar = getBamCigar(record);
            read_record.read_id_next = getBamMateRef(record, hdr);
            read_record.pos_next = getBamMateRefPos(record);
            read_record.temp_len = 0;
            read_record.read_seq = getBamSeq(record);
            read_record.read_quality = getBamQuanlity(record);

            string cur_attributes;
            stringstream string_in(getBamTag(record));
            while( string_in.good() )
            {
                string_in >> cur_attributes;
                read_record.attributes.push_back(cur_attributes);
            }
        }
        cout << read_record;
        bam_destroy1(record);
    }

    bgzf_close(fn_hd);

    return 0;
}



