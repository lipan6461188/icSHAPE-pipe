
#include "fasta.h"

using namespace std;

namespace pan{

// **************************
//  Fasta class
// **************************

void Fasta::clear()
{
    this->chr_ids.clear();
    this->sequence.clear();
    this->chr_annotation.clear();
}

void Fasta::read_fasta_file(const string &fastaFn, bool strict)
{
    ifstream SEQ(fastaFn, ifstream::in);
    if(not SEQ){
        throw runtime_error( "Bad_Input_File: "+fastaFn );
    }

    string this_line;  // Current line content
    string cur_chrID;  // Current chrID
    while(getline(SEQ, this_line))
    {
        if(this_line.empty()){ continue; }
        if(this_line[0] == '>')
        {
            istringstream head_line(this_line);
            
            // >chr1 The first chromosome
            // >chr1
            head_line >> cur_chrID;
            cur_chrID = cur_chrID.substr(1);
            auto pos = this->sequence.find(cur_chrID);
            
            // A duplicate sequence
            if( pos != this->sequence.end() )
            {
                if(strict)
                    throw Unexpected_Error("Error: duplicate Fatsa Sequence: "+cur_chrID);
                else
                {
                    cerr << "Warning: " << "duplicate fatsa sequence: " << cur_chrID << "; only preserve the last one" << endl;
                    this->sequence[cur_chrID] = "";
                    this->chr_annotation[cur_chrID] = "";
                }
            }else{
                // A new sequence, init it
                this->sequence[cur_chrID];

                string annotation;
                while( head_line.good() )
                {
                    string cur_annotation;
                    head_line >> cur_annotation;
                    annotation.append( (annotation.empty() ? "" : " ")+cur_annotation );
                }
                this->chr_annotation[cur_chrID] = annotation;
                this->chr_ids.push_back(cur_chrID);
            }
        }else{
            this->sequence[cur_chrID].append(this_line);
        }
    }
    SEQ.close();
}

Fasta::Fasta(const string &fastaFn, bool strict)
{
    this->read_fasta_file(fastaFn, strict);
}

void Fasta::add_fasta(const string &fastaFn, bool strict)
{
    this->read_fasta_file(fastaFn, strict);
}

StringArray Fasta::get_chr_ids() const
{
    return this->chr_ids;
}

uLONG Fasta::get_chr_len(const string &chrID) const
{
    return this->sequence.at(chrID).size();
}

MapStringuLONG Fasta::get_chr_lens() const
{
    MapStringuLONG len;
    for(auto iter=this->sequence.cbegin(); iter!=this->sequence.cend(); iter++)
        len[iter->first] = iter->second.size();
    return len;
}

uLONG Fasta::get_chr_num() const
{
    return this->sequence.size();
}

const string &Fasta::get_chr_seq(const string &chrID) const
{
    return this->sequence.at(chrID);
}

const string &Fasta::get_chr_anno(const string &chrID) const
{
    return this->chr_annotation.at(chrID);
}

string Fasta::get_chr_subbseq(const string &chrID, uLONG start, uLONG len, STRAND strand) const
{
    if(strand == POSITIVE)
        return this->sequence.at(chrID).substr(start, len);
    else{
        return reverse_comp(this->sequence.at(chrID).substr(start, len));
    }
}

// **************************
//  Quick fasta class
// **************************

void qFasta::build_fasta_index(const string &fastaFn)
{
    // index file
    string fastaIndFn = fastaFn + ".fai";

    // input handle and output handle
    ifstream SEQ(fastaFn, ifstream::in);
    if(not SEQ){
        throw runtime_error( "Bad_Input_File: "+fastaFn );
    }

    ofstream INDEX(fastaIndFn, ofstream::out);
    if(not INDEX){
        throw runtime_error( "Bad_Output_File: "+fastaFn );
    }

    // chromosome list
    set<string> chr_set;

    // current file position
    uLONGLONG cur_pos = 0;

    // current line count
    uLONG cur_line_count = 0;
    
    // current line cotent
    string this_line;

    // current chromosome id
    string cur_chrID;

    // current chromosome length
    uLONG cur_chrLen = 0;

    // current chromosome start position
    uLONG cur_chrStart = 0;

    // current line length for chromosome
    uLONG chr_line_len = 0;

    // flag
    bool first_line = false;
    bool last_line = false;

    while(getline(SEQ, this_line))
    {
        cur_line_count++;

        auto line_size = this_line.size();

        if(line_size == 0)
        {
            // do nothing
        }else if(this_line[0] == '>')
        {
            if(cur_chrLen)
            {
                // Output current chromosome
                INDEX << cur_chrID << "\t" 
                    << cur_chrLen << "\t" 
                    << cur_chrStart << "\t" 
                    << chr_line_len << "\t" 
                    << chr_line_len+1 << "\n";
                cur_chrLen = 0;
            }

            // Read a new chromosome
            istringstream head_line(this_line);
            head_line >> cur_chrID;
            cur_chrID = cur_chrID.substr(1);

            if(chr_set.find(cur_chrID) != chr_set.end())
                cerr << "Warning: duplicate chrID: " << cur_chrID << " found" << endl;
            
            chr_set.insert(cur_chrID);
            first_line = true;
            last_line = false;
        }else{
            if(last_line == true)
                throw Unexpected_Error( "FATAL Error: different line count for "+cur_chrID+" in line "+to_string(cur_line_count-1) );

            if(first_line)
            {
                cur_chrStart = cur_pos;
                chr_line_len = line_size;
                first_line = false;
            }

            if(line_size != chr_line_len)
                last_line = true;

            cur_chrLen += line_size;
        }

        cur_pos += line_size + 1;
    }

    if(cur_chrLen)
        INDEX << cur_chrID << "\t" << cur_chrLen << "\t" << cur_chrStart << "\t" << chr_line_len << "\t" << chr_line_len+1 << "\n";

    INDEX.close();
    SEQ.close();
}

void qFasta::load_fasta_file(const string &fastaFn)
{
    this->FASTA.open(fastaFn, ifstream::binary);
    if(not this->FASTA)
        throw runtime_error( "Bad_Input_File: "+fastaFn );

    ifstream IN(fastaFn+".fai", ifstream::in);
    if(not IN)
    {
        cerr << "Warning: fai index not found, now build " << fastaFn+".fai" << "..." << endl;
        this->build_fasta_index(fastaFn);
        IN.open(fastaFn+".fai");
    }

    string this_line;
    while(getline(IN, this_line))
    {
        fai_item fi;
        StringArray vec;
        split(this_line, vec);
        fi.chr_len = std::stol(vec.at(1));
        fi.chr_start = std::stol(vec.at(2));
        fi.chr_line_len = std::stol(vec.at(3));
        this->fa_index[vec.at(0)] = fi;
    }

    IN.close();
}

qFasta::qFasta(const string &fastaFn)
{
    this->load_fasta_file(fastaFn);
}

StringArray qFasta::get_chr_ids() const
{
    StringArray chr_ids;
    for(auto it=this->fa_index.cbegin();it!=this->fa_index.cend();it++)
        chr_ids.push_back(it->first);
    return chr_ids;
}

uLONG qFasta::get_chr_len(const string &chrID) const
{
    return this->fa_index.at(chrID).chr_len;
}

MapStringuLONG qFasta::get_chr_lens() const
{
    MapStringuLONG chr_lens;
    for(auto it=this->fa_index.cbegin();it!=this->fa_index.cend();it++)
    {
        chr_lens[ it->first ] =  it->second.chr_len;
    }
    return chr_lens;
}

uLONG qFasta::get_chr_num() const
{
    return this->fa_index.size();
}

string qFasta::get_chr_seq(const string &chrID)
{
    return get_chr_subbseq(chrID, 0, string::npos/2, POSITIVE);
}

string qFasta::get_chr_subbseq(const string &chrID, uLONG start, uLONG len, STRAND strand)
{
    // The length shouldn't exceed the chr length
    uLONG chr_len = this->get_chr_len(chrID);
    len = (start+len > chr_len) ? chr_len-start : len;

    const fai_item index = this->fa_index.at(chrID);
    
    // Locate to the position
    uLONG offset = index.chr_start + start + start / index.chr_line_len;
    this->FASTA.seekg(offset, ios_base::beg);

    // Read the block
    uLONG read_len = len+len/index.chr_line_len;
    char *buffer = new char[read_len+10];
    // !! No \0 append to buffer in read function
    this->FASTA.read(buffer, read_len+5);
    if(this->FASTA.eof())
        this->FASTA.clear();

    string read_seq;
    uLONG seq_len = 0;

    for(uLONG i=0; seq_len!=len; i++)
    {
        if(buffer[i] != '\n')
        {
            read_seq += buffer[i];
            ++seq_len;
        }
    }
    delete[] buffer;

    if(read_seq.size() != len)
    {
        cerr << "FATAL Error: different size " << read_seq.size() << "\t" << len << endl;
    }

    if(strand == POSITIVE)
        return read_seq;
    else
        return reverse_comp(read_seq);
}

// **************************
//  Other common functions
// **************************


string reverse_comp(const string &raw_seq, NUCTYPE ntype)
{
    string rev_comp_seq;
    for(auto iter=raw_seq.crbegin(); iter!=raw_seq.crend(); iter++)
    {
        switch(*iter)
        {
            case 'A':
                if(ntype==DNA)
                    rev_comp_seq += 'T';
                else if(ntype==RNA)
                    rev_comp_seq += 'U';
                break;
            case 'a':
                if(ntype==DNA)
                    rev_comp_seq += 't';
                else if(ntype==RNA)
                    rev_comp_seq += 'u';
                break;

            case 'T': case 'U':
                rev_comp_seq += 'A';
                break;
            case 't': case 'u':
                rev_comp_seq += 'a';
                break;

            case 'C':
                rev_comp_seq += 'G';
                break;
            case 'c':
                rev_comp_seq += 'g';
                break;

            case 'G':
                rev_comp_seq += 'C';
                break;
            case 'g':
                rev_comp_seq += 'c';
                break;

            default:
                rev_comp_seq += 'N';
        }
    }
    return rev_comp_seq;
}



string flat_seq(const string &raw_seq, size_t line_width)
{
    if(line_width >= raw_seq.size())
        return raw_seq+"\n";

    string flatted_seq;
    decltype(raw_seq.size()) idx=0;
    for(; idx<=raw_seq.size()-line_width; idx+=line_width)
    {
        flatted_seq += raw_seq.substr(idx, line_width);
        flatted_seq += "\n";
    }
    if(idx != raw_seq.size())
    {
        flatted_seq += raw_seq.substr(idx);
        flatted_seq += "\n";
    }
    return flatted_seq;
}



int global_align_2_seq(const string &a,
        const string &b,
        string &a_aligned, 
        string &b_aligned,
        const int gap_penalty,
        const int mismatch_penalty)
{
    a_aligned.clear();
    b_aligned.clear();

    struct DP_node
    {
        enum DIRECT{ left, top, left_top };
        int penalty;
        DIRECT direct;
        
        DP_node(int p=0, DIRECT d=left_top): penalty(p), direct(d){}
    };

    size_t n = a.size();
    size_t m = b.size();

    vector<vector<DP_node>> A(n + 1, vector<DP_node>(m + 1));

    for (size_t i = 0; i <= m; ++i)
        A[0][i] = DP_node(gap_penalty * i, DP_node::top);
    for (size_t i = 0; i <= n; ++i)
        A[i][0] = DP_node(gap_penalty * i, DP_node::left);

    for (size_t i = 1; i <= n; ++i)
    {
        for (size_t j = 1; j <= m; ++j)
        {
            char x_i = a[i-1];
            char y_j = b[j-1];

            int penalty_lt = A[i-1][j-1].penalty + ( (x_i == y_j) ? 0 : mismatch_penalty );
            int penalty_l = A[i-1][j].penalty + gap_penalty;
            int penalty_t = A[i][j-1].penalty + gap_penalty;

            if(penalty_lt <= penalty_l and penalty_lt <= penalty_t)
            {
                A[i][j] = DP_node(penalty_lt, DP_node::left_top);
            }else if(penalty_l <= penalty_lt and penalty_l <= penalty_t)
            {
                A[i][j] = DP_node(penalty_l, DP_node::left);
            }else if(penalty_t <= penalty_lt and penalty_t <= penalty_l)
            {
                A[i][j] = DP_node(penalty_t, DP_node::top);
            }else{
                // no thing
            }
        }
    }

    size_t j = m;
    size_t i = n;
    int a_count = 0, b_count = 0;
    for (; i >= 1 and j >= 1; )
    {
        char x_i = a[i-1];
        char y_j = b[j-1];
        if (A[i][j].direct == DP_node::left_top)
        {
            a_aligned.push_back(x_i);
            b_aligned.push_back(y_j);
            --j; --i;
            a_count++; b_count++;
        }
        else if (A[i][j].direct == DP_node::left)
        {
            a_aligned.push_back(x_i);
            b_aligned.push_back('-');
            --i;
            a_count++;
        }
        else if(A[i][j].direct == DP_node::top) 
        {
            a_aligned.push_back('-');
            b_aligned.push_back(y_j);
            --j;
            b_count++;
        }else{
            cerr << "Error" << endl;
        }
    }

    while (i >= 1 && j < 1)
    {
        a_aligned.push_back(a[i-1]);
        b_aligned.push_back('-');
        --i;
        a_count++;
    }
    while (j >= 1 && i < 1)
    {
        a_aligned.push_back('-');
        b_aligned.push_back(b[j-1]);
        --j;
        b_count++;
    }

    reverse(a_aligned.begin(), a_aligned.end());
    reverse(b_aligned.begin(), b_aligned.end());

    return A[n][m].penalty;
}

}
