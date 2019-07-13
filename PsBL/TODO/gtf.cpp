
#include "gtf.h"

using namespace std;

namespace pan{





sp<GTF_Line> read_a_gtf_line(ifstream &IN)
{
    sp<GTF_Line> sp_gtf_line = make_shared<GTF_Line>();

    string cur_line;
    if(getline(IN, cur_line))
    {
        trim(cur_line);
        if(not cur_line.empty() and cur_line[0] != '#')
        {
            istringstream IN_Stream(cur_line);
            string tmp_score, tmp_frame;
            IN_Stream >> sp_gtf_line->chr_id >> sp_gtf_line->source >> sp_gtf_line->feature >> sp_gtf_line->start >> sp_gtf_line->end
                      >> tmp_score >> sp_gtf_line->strand >> tmp_frame;
            if(not IN_Stream)
                return nullptr;

            if(tmp_score != ".")
                sp_gtf_line->score = stod(tmp_score);

            if(tmp_frame == "0")
                sp_gtf_line->frame = GTF_Line::ZERO;
            else if(tmp_frame == "1")
                sp_gtf_line->frame = GTF_Line::ONE;
            else if(tmp_frame == "2")
                sp_gtf_line->frame = GTF_Line::TWO;
            else
                sp_gtf_line->frame = GTF_Line::NO;

            string attr_line;
            getline(IN_Stream, attr_line);
            istringstream attr_Stream(attr_line);
            while(attr_Stream)
            {
                string attr_key;
                string attr_value;
                attr_Stream >> attr_key >> attr_value;
                trim(attr_value, ';');
                trim(attr_value, '\"');
                sp_gtf_line->attributes[attr_key] = attr_value;
            }

            sp_gtf_line->chr_id.shrink_to_fit();
            sp_gtf_line->source.shrink_to_fit();
            sp_gtf_line->feature.shrink_to_fit();
            //for(pair<string, string> &item: sp_gtf_line->attributes)
            //  item.first.shrink_to_fit(), item.second.shrink_to_fit();

            return sp_gtf_line;
        }else
            return nullptr;
    }else
        return nullptr;
}


bool operator< (const Base_GTF_Exon &exon_1, const Base_GTF_Exon &exon_2)
{
    if(exon_1.strand() == exon_2.strand())
    {
        if( exon_1.strand() == '+')
            return exon_1.genome_start() < exon_2.genome_start();
        else
            return exon_1.genome_start() > exon_2.genome_start();
    }else{
        return exon_1.strand() < exon_2.strand();
    }   
}

bool operator< (const Base_GTF_CDS &cds_1, const Base_GTF_CDS &cds_2)
{
    if(cds_1.strand() == cds_2.strand())
    {
        if( cds_1.strand() == '+')
            return cds_1.genome_start() < cds_2.genome_start();
        else
            return cds_1.genome_start() > cds_2.genome_start();
    }else{
        return cds_1.strand() < cds_2.strand();
    }   
}

void GTF2_Gencode::read_file(const string &gtf2_file_name, const GTF2_Parse_Key * const gtf_parser)
{
    ifstream IN(gtf2_file_name);
    if(not IN)
        throw Bad_IO("Bad GTF File", false);
    
    MapStringT< vector<sp<GTF2_CDS>> > trans_cds;
    MapStringT< vector<sp<GTF2_Exon>> > trans_exon;
    uLONG gene_count(0), trans_count(0), exon_count(0), cds_count(0), other_count(0), line_count(0);

    clog << "line_count\tgene_count\ttrans_count\texon_count\tcds_count\tother_count\n";
    while(IN)
    {
        sp<GTF_Line> sp_gtf_line = read_a_gtf_line(IN);
        if(sp_gtf_line)
        {
            line_count++;
            if(line_count % 100000 == 0)
                clog << line_count << "\t" << gene_count << "\t" << trans_count << "\t" << exon_count << "\t" << cds_count << "\t" << other_count << "\n";

            if(sp_gtf_line->feature == gtf_parser->gene_feature_key)
            {
                sp<GTF2_Gene> gene = make_shared<GTF2_Gene>(sp_gtf_line, gtf_parser);
                //this->genes.emplace(gene.gene_id(), sp_gtf_line, gtf_parser);
                this->genes[ gene->gene_id() ] = gene;
                gene_count++;
            }
            else if(sp_gtf_line->feature == gtf_parser->transcript_feature_key)
            {
                sp<GTF2_Transcript> transcript = make_shared<GTF2_Transcript>(sp_gtf_line, gtf_parser);
                //this->transcripts.emplace(transcript.trans_id(), sp_gtf_line, gtf_parser);
                this->transcripts[ transcript->trans_id() ] = transcript;
                trans_count++;
            }
            else if(sp_gtf_line->feature == gtf_parser->exon_feature_key)
            {
                sp<GTF2_Exon> exon = make_shared<GTF2_Exon>(sp_gtf_line, gtf_parser);
                //this->transcripts.emplace(transcript.gene_id(), sp_gtf_line, gtf_parser);
                trans_exon[ exon->trans_id() ].push_back(exon);
                exon_count++;
            }
            else if(sp_gtf_line->feature == gtf_parser->cds_feature_key)
            {   
                sp<GTF2_CDS> cds = make_shared<GTF2_CDS>(sp_gtf_line, gtf_parser);
                trans_cds[ cds->trans_id() ].push_back(cds);
                cds_count++;
            }else{
                // other features are ignored
                other_count++;
            }
        }
    }

    for(const pair<string, vector<sp<GTF2_Exon>>> &item: trans_exon)
    {
        this->transcripts.at(item.first)->exon_list = item.second;
        // sort exons
        sort( this->transcripts.at(item.first)->exon_list.begin(), this->transcripts.at(item.first)->exon_list.end(), [](const sp<GTF2_Exon> &exon_1, const sp<GTF2_Exon> &exon_2){ return *exon_1 < *exon_2; } );
    }

    for(const pair<string, vector<sp<GTF2_CDS>>> &item: trans_cds)
    {
        this->transcripts.at(item.first)->cds_list = item.second;
        // sort cds
        sort( this->transcripts.at(item.first)->cds_list.begin(), this->transcripts.at(item.first)->cds_list.end(), [](const sp<GTF2_CDS> &cds_1, const sp<GTF2_CDS> &cds_2){ return *cds_1 < *cds_2; } );
    }

    sp<MapStringT<vector<sp<GTF2_Transcript>>>> gene_trans_map = get_gene_trans_map();
    for(const pair<string, vector<sp<GTF2_Transcript>>> &item: *gene_trans_map)
    {
        this->genes.at(item.first)->trans_list = item.second;
        for(const sp<GTF2_Transcript> &trans_item: item.second)
            trans_item->gene = this->genes.at(item.first);
    }
}

sp<StringArray> GTF2_Gencode::get_gene_list() const
{
    sp<StringArray> sp_gene_list = make_shared<StringArray>();
    for_each(genes.cbegin(), genes.cend(), [&](const pair<string, sp<GTF2_Gene>> &gene_item){ sp_gene_list->push_back(gene_item.first); });
    return sp_gene_list;
}

sp<StringArray> GTF2_Gencode::get_trans_list() const
{
    sp<StringArray> sp_trans_list = make_shared<StringArray>();
    for_each(transcripts.cbegin(), transcripts.cend(), [&](const pair<string, sp<GTF2_Transcript>> &trans_item){ sp_trans_list->push_back(trans_item.first); });
    return sp_trans_list;
}

sp<MapStringT<vector<sp<GTF2_Transcript>>>> GTF2_Gencode::get_gene_trans_map() const
{
    sp<MapStringT< vector<sp<GTF2_Transcript> >>> sp_gene_trans_map = make_shared< MapStringT<vector<sp<GTF2_Transcript>>> >();

    for(const pair<string, sp<GTF2_Transcript>> &trans_item: transcripts)
    {
        (*sp_gene_trans_map)[ trans_item.second->uniq_gene_id() ].push_back( trans_item.second );
    }

    return sp_gene_trans_map;
}

bool GTF2_Gencode::has_gene(const string &gene_id) const
{
    if(genes.find(gene_id) == genes.cend())
        return false;
    else
        return true;
}

bool GTF2_Gencode::has_trans(const string &trans_id) const
{
    if(transcripts.find(trans_id) == transcripts.cend())
        return false;
    else
        return true;
}


const sp<GTF2_Gene> GTF2_Gencode::get_gene_handle(const string &gene_id) const
{
    return genes.at(gene_id);
}

const sp<GTF2_Transcript> GTF2_Gencode::get_trans_handle(const string &trans_id) const
{
    return transcripts.at(trans_id);
}

void GTF2_Gencode::write_genome_coor_annotation(ostream &OUT) const
{
    for( const pair<string, sp<GTF2_Transcript>> &trans_item: transcripts )
    {

        OUT << trans_item.second->chr_id() << "\t" << trans_item.second->genome_start() << "\t" 
            << trans_item.second->genome_end() << "\t" << trans_item.second->strand() << "\t"
            << trans_item.second->gene->gene_name() << "=" << trans_item.second->gene->gene_id() << "\t"
            << trans_item.second->trans_id() << "\t" << trans_item.second->trans_type() << "\t";

        for(auto iter=trans_item.second->exon_list.cbegin(); iter!=trans_item.second->exon_list.cend(); iter++)
        {
            OUT << (*iter)->genome_start() << "-" << (*iter)->genome_end();
            if(iter!=trans_item.second->exon_list.cend()-1)
                OUT << ",";
        }

        OUT << "\t";
        for(auto iter=trans_item.second->cds_list.cbegin(); iter!=trans_item.second->cds_list.cend(); iter++)
        {
            OUT << (*iter)->genome_start() << "-" << (*iter)->genome_end();
            if(iter!=trans_item.second->cds_list.cend()-1)
                OUT << ",";
        }
        OUT << "\n";
    }
}






















































}