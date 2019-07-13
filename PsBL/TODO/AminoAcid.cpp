
struct Amino_Acid
{
    string animo_acid;
    string three_letter;
    string one_letter;
    
    string side_chain_name;
    string side_chain_polarity;
    string side_chain_charge;
    double hydropathy_index;
};

ostream& operator<<(ostream& OUT, const Amino_Acid &aa);

void load_amino_acid(MapStringT<Amino_Acid> &aa_vec, 
        const string &file_name="conf_data/amino_acid.txt");

void load_codon_map(MapStringString &codon_map, 
        const string &file_name="conf_data/codon_table.txt");


void load_amino_acid(MapStringT<Amino_Acid> &aa_vec, const string &file_name)
{
    ifstream IN(file_name, ifstream::in);
    if(not IN)
    {
        throw runtime_error( "Bad_Input_File: "+file_name );
    }

    aa_vec.clear();
    string cur_line;
    while(getline(IN, cur_line))
    {
        if(cur_line.empty())
            continue;
        istringstream line_in(cur_line);
        Amino_Acid aa;
        line_in >>  aa.animo_acid >> aa.three_letter >> aa.one_letter >> 
                    aa.side_chain_name >> aa.side_chain_polarity >> aa.side_chain_charge >>
                    aa.hydropathy_index;
        aa_vec[ aa.three_letter ] = aa;
    }
    IN.close();
}

void load_codon_map(MapStringString &codon_map, const string &file_name)
{
    ifstream IN(file_name, ifstream::in);
    if(not IN)
    {
        throw runtime_error( "Bad_Input_File: "+file_name );
    }

    codon_map.clear();
    string cur_line;
    while(getline(IN, cur_line))
    {
        if(cur_line.empty())
            continue;

        auto aa_codon = split(cur_line, '=');
        if(aa_codon.size() != 2)
        {
            throw runtime_error("Bad Codon_table " + file_name);
        }
        auto codons = split(aa_codon.at(1), '|');
        for(string codon: codons)
        {
            codon_map[ codon ] = aa_codon.at(0);
            std::replace(codon.begin(), codon.end(), 'U', 'T');
            codon_map[ codon ] = aa_codon.at(0);
        }
        
    }
    IN.close();
}

ostream& operator<<(ostream& OUT, const Amino_Acid &aa)
{
    OUT << aa.animo_acid << "\t" << aa.three_letter << "\t" << aa.one_letter
        << "\t" << aa.side_chain_name << "\t" << aa.side_chain_polarity 
        << "\t" << aa.side_chain_charge << "\t" << aa.hydropathy_index;

    return OUT;
}






StringArray translation(const string &raw_seq, const MapStringString &codon_map);
StringArray translation(const string &raw_seq, const MapStringString &codon_map)
{
    StringArray protein;
    for(size_t idx=0; idx<raw_seq.size()-2; idx+=3)
    {
        protein.push_back( codon_map.at(raw_seq.substr(idx, 3)) );
    }
    return protein;
}






