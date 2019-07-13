#include "RsampleData.h"
#include "src/rna_library.h"
#include <cmath>
#include <cstring>
#include <stdarg.h>
#include <fstream>
#include <iostream>
using namespace std;

#define ERR_BAD_FILE_UU 1
#define ERR_BAD_FILE_PE 2
#define ERR_BAD_FILE_PM 3

//constuctor
//! check ErrorCode after construction. (2 == bad file)
RsampleData::RsampleData(const char* const uu_file, const char* const pend_file, const char* const  pmid_file) {
	double reactivity;
	ErrorCode = 0;
	
  string s_uu_file =   (uu_file   == NULL ? "" : uu_file);
  string s_pend_file = (pend_file == NULL ? "" : pend_file);
  string s_pmid_file = (pmid_file == NULL ? "" : pmid_file);

// unpaired reference data
	if (s_uu_file.empty()) {
	  s_uu_file = getDataPath();
	  s_uu_file += "/rsample/unpaired.txt";
	} 

// paired_end reference data
	if (s_pend_file.empty()) {
	  s_pend_file = getDataPath();
	  s_pend_file += "/rsample/paired_end.txt";
	} 

// paired_mid reference data
	if (s_pmid_file.empty()) {
	  s_pmid_file = getDataPath();
	  s_pmid_file += "/rsample/paired_mid.txt";
	} 


ifstream infile;

 // unpaired uu
  infile.open( s_uu_file.c_str() );
 if (!infile.good()) {ErrorCode = ERR_BAD_FILE_UU; return;}
  while(infile >> reactivity ) {
    react_uu.push_back(reactivity);
  }
  infile.close();
  
  
  // paired-end i.e. pend
  infile.open(s_pend_file.c_str());
  if (!infile.good()) {ErrorCode = ERR_BAD_FILE_PE; return;}
  while(infile >> reactivity ) {
    react_pend.push_back(reactivity);
  }
  infile.close();
  
  // paired-mid i.e. pmid
  infile.open(s_pmid_file.c_str());
  if (!infile.good()) {ErrorCode = ERR_BAD_FILE_PM; return;}
  while(infile >> reactivity ) {
    react_pmid.push_back(reactivity);
  }
  infile.close();	
}

const char* RsampleData::GetErrorMessage(const int code) {
  #define BAD_RSAMPLE_FILE_MSG(NUC_TYPE) "The Rsample reactivity data for " NUC_TYPE " could not be opened or read."
  switch(code) {
    case ERR_BAD_FILE_UU: return BAD_RSAMPLE_FILE_MSG("unpaired nucleotides");
    case ERR_BAD_FILE_PE: return BAD_RSAMPLE_FILE_MSG("paired-end nucleotides");
    case ERR_BAD_FILE_PM: return BAD_RSAMPLE_FILE_MSG("paired-middle nucleotides");
    default: return "Unknown error code.";
  }
  #undef BAD_RSAMPLE_FILE_MSG
}
