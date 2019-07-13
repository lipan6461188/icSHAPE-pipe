
#ifndef HTSLIB_H
#define HTSLIB_H

#include <htslib/bgzf.h>
#include <htslib/sam.h>
#include <iostream>
#include <string>
#include "pan_type.h"

namespace pan{

using std::string;
using std::to_string;


/* 
A read:
SRR3194440.15842008 32  ENST00000527779.1   495 1   5S10M100N10M6S  *   0   0   TTGTGAAGACATTGGTGTTGAACCTGAAAAT JJIHIJJJGJJJJJJHHHIIHIIJJJIJJJJ AS:i:0  XS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:31 YT:Z:UU
*/

// **************************
//  Some functions to parse information from Bam file
// **************************

/*

#include <sam.h>
#include <htslib.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>

## Common framework
BGZF* fn_hd = bgzf_open(argv[1], "r");
			=>
bam_hdr_t *hdr = bam_hdr_read(fn_hd);
			=>
bool read_a_read_record(BGZF* fn_hd, bam_hdr_t *hdr, vector<Sam_Record> &read_records, Sam_Record* &p_cache);
			=>
bgzf_close(fn_hd);

## Build a record
bam1_t *record = bam_init1();
int ret = bam_read1(fn_hd, record);

*/

/* 
record             -- One record from bam file
*/
string getBamQName(bam1_t *record);  					// Get read ID
uint16_t getBamFlag(bam1_t *record);					// Get read flag
string getBamRef(bam1_t *record, bam_hdr_t *hdr);		// Get reference ID
int32_t getBamRefPos(bam1_t *record);					// Get reference map position
int32_t getBamMapQuanlity(bam1_t *record);				// Get mapping quality
string getBamCigar(bam1_t *record);						// Get ciger
string getBamMateRef(bam1_t *record, bam_hdr_t *hdr);	// Get read mate reference
int32_t getBamMateRefPos(bam1_t *record);				// Get read mate reference map position
/*** TLEN skip ***/
string getBamSeq(bam1_t *record);						// Get read sequence
string getBamQuanlity(bam1_t *record);					// Get read quanlity
string getBamTag(bam1_t *record);						// Get bam tags -- TODO: get full attribute list

bool getBamHead(bam_hdr_t *hdr, MapStringT<uLONG> &chr_len);	// Get bam head

};
#endif