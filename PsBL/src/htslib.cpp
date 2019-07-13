#include "htslib.h"

namespace pan{

string getBamQName(bam1_t *record)
{
    return bam_get_qname(record);
}

string getBamSeq(bam1_t *record)
{
    auto seq_l = record->core.l_qseq;
    auto codedSeq = bam_get_seq(record);

    char *seq = new char[seq_l+1];
    seq[seq_l] = '\0';

    string reads;
    for(int i=0; i<seq_l; i++)
    {
        auto base = bam_seqi(codedSeq, i);
        if(base == 1)
            seq[i] = 'A';
        else if(base == 2)
            seq[i] = 'C';
        else if(base == 4)
            seq[i] = 'G';
        else if(base == 8)
            seq[i] = 'T';
        else
            seq[i] = 'N';
    }
    
    string seq_str(seq);
    delete[] seq;
    
    return seq_str;
}

string getBamQuanlity(bam1_t *record)
{
    auto seq_l = record->core.l_qseq;
    auto coded_quanlity = bam_get_qual(record);

    char *quan = new char[seq_l+1];
    quan[seq_l] = '\0';

    
    for(int i=0; i<seq_l; i++)
    {
        quan[i] = coded_quanlity[i]+33;
    }

    string quanlity(quan);
    delete[] quan;

    return quanlity;
}

uint16_t getBamFlag(bam1_t *record)
{
    return record->core.flag;
}

string getBamRef(bam1_t *record, bam_hdr_t *hdr)
{
    int32_t tid = record->core.tid;
    return hdr->target_name[tid];
}

int32_t getBamRefPos(bam1_t *record)
{
    return record->core.pos + 1;
}

int32_t getBamMapQuanlity(bam1_t *record)
{
    return record->core.qual;
}

string getBamTag(bam1_t *record)
{
    // TODO: add more tags

    uint8_t *NH = bam_aux_get(record, "NH");
    uint8_t *HI = bam_aux_get(record, "HI");
    uint8_t *AS = bam_aux_get(record, "AS");
    uint8_t *MD = bam_aux_get(record, "MD");
    uint8_t *NM = bam_aux_get(record, "NM");

    string tag;

    if(NH)
    {
        if(not tag.empty()) tag += "\t";
        tag += "NH:i:" + to_string( bam_aux2i(NH) );
    }
    if(HI)
    {
        if(not tag.empty()) tag += "\t";
        tag += "HI:i:" + to_string( bam_aux2i(HI) );
    }
    if(AS)
    {
        if(not tag.empty()) tag += "\t";
        tag += "AS:i:" + to_string( bam_aux2i(AS) );
    }
    if(MD)
    {
        if(not tag.empty()) tag += "\t";
        tag += "MD:i:" + string( bam_aux2Z(MD) );
    }
    if(NM)
    {
        if(not tag.empty()) tag += "\t";
        tag += "NM:i:" + to_string( bam_aux2i(NM) );
    }

    return tag;
}

string getBamMateRef(bam1_t *record, bam_hdr_t *hdr)
{
    int32_t mtid = record->core.mtid;
    if(mtid == -1)
        return "*";
    else
        return hdr->target_name[mtid];
}

int32_t getBamMateRefPos(bam1_t *record)
{
    return record->core.mpos + 1;
}

string getBamCigar(bam1_t *record)
{
    uint32_t n_cigar = record->core.n_cigar;
    auto codedCigar = bam_get_cigar(record);

    string cigar;
    for(int i=0; i<n_cigar; i++)
    {
        auto op = bam_cigar_op(codedCigar[i]);
        auto op_len = bam_cigar_oplen(codedCigar[i]);

        if(op == BAM_CMATCH)
            cigar += to_string(op_len)+'M';
        else if(op == BAM_CINS)
            cigar += to_string(op_len)+'I';
        else if(op == BAM_CDEL)
            cigar += to_string(op_len)+'D';
        else if(op == BAM_CREF_SKIP)
            cigar += to_string(op_len)+'N';
        else if(op == BAM_CSOFT_CLIP)
            cigar += to_string(op_len)+'S';
        else if(op == BAM_CHARD_CLIP)
            cigar += to_string(op_len)+'H';
        else if(op == BAM_CPAD)
            cigar += to_string(op_len)+'P';
        else if(op == BAM_CEQUAL)
            cigar += to_string(op_len)+'=';
        else if(op == BAM_CDIFF)
            cigar += to_string(op_len)+'X';
        else
        { cigar += to_string(op_len)+"?"; }
    }   

    return cigar;
}


bool getBamHead(bam_hdr_t *hdr, MapStringT<uLONG> &chr_len)
{
    chr_len.clear();
    for(int i=0; i<hdr->n_targets; i++)
        chr_len[hdr->target_name[i]] = hdr->target_len[i];

    return true;
}




};