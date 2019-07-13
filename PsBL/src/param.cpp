
#include "param.h"
#include "string_split.h"
#include "exceptions.h"

namespace pan{

void push_params(int &argc, char *argv[], const string &param)
{
    argc++;
    argv[argc-1] = new char[20];
    strcpy(argv[argc-1], param.c_str());
}

std::string param_string(const int argc, char *argv[])
{
    string param;
    for(int i=0; i<argc; i++)
    {
        param += string( argv[i] );
        if(i != argc-1)
            param += " ";
    }
    return param;
}


namespace FILE_FORMAT
{

file_format_type guess_file_type(const string &file_name)
{
    StringArray file_items;
    split(file_name, '.', file_items);
    string postfix = file_items.back();

    if(postfix == "sam")
        return SAM_FILE;
    else if(postfix == "bam")
        return BAM_FILE;
    else if(postfix == "fasta" or postfix == "fa")
        return FASTA_FILE;
    else if(postfix == "fastq" or postfix == "fq")
        return FASTQ_FILE;
    else if(postfix == "stockholm" or postfix == "sto")
        return STOCKHOLM_FILE;
    else if(postfix == "matrix")
        return MATRIX_FILE;
    else if(postfix == "cm")
        return CM_FILE;
    else if(postfix == "dg")
        return DG_FILE;
    else
        return UNKNOWM;
}

string fft_to_string(const file_format_type &fft)
{
    switch(fft)
    {
        case SAM_FILE:
            return "sam";
            break;
        case BAM_FILE:
            return "bam";
            break;
        case FASTA_FILE:
            return "fasta";
            break;
        case FASTQ_FILE:
            return "fastq";
            break;
        case STOCKHOLM_FILE:
            return "stockholm";
            break;
        case CM_FILE:
            return "cm";
            break;
        case MATRIX_FILE:
            return "matrix";
            break;
        case DG_FILE:
            return "dg";
            break;
        case UNKNOWM:
            return "unknown";
            break;
        
        default:
            throw Unexpected_Error("Unkown file format");
    }
}

}

const string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return buf;
}


/*
    enum file_format_type
    {
        SAM_FILE,
        BAM_FILE,
        FASTA_FILE,
        FASTQ_FILE,
        STOCKHOLM_FILE,
        CM_FILE,

        // I defined type
        MATRIX_FILE,
        DG_FILE,

        // UNKNOWM
        UNKNOWM
    };

*/







}