#if !defined(INTERMOLECULAR_H)
#define INTERMOLECULAR_H


/*=======================================================================
intermolecular.h and intermolecular.cpp inculde funcions calculating and report
different free energy for binding in OligoWalk.
They are revised based on Mathews' code from RNAStructure.
olig() generate the energy data;
report save the generated data;
siprefileter and filterbysirna() were used to prefilter and postfileter the 
functional siRNA respectively, using criterias other than thermodynamics


															----Nov., 2005
															John(Zhi Lu)

=======================================================================*/

#include <cmath>
//#include "stdafx.h"
#include "rna_library.h"
#include "thermo.h"
#include "pclass.h"
#include "algorithm.h"
#include "siRNAfilter.h"
#include "alltrace.h"
#include "alltrace_intermolecular.h"
#include "stochastic.h"


#define FILTER_PASS 6	//The criteria for prefiltering of functional siRNA in OligoWalk
#define maxfile 75		//maximum length of file names

//=======================================================================
//this structure rddata contains all the thermodynamic parameters needed
//for determining the stability of an RNA-DNA duplex
struct rddata {

   short int stack[5][5][5][5];
   short int init;
};


//=======================================================================
//olig is the backend function for oligo walk
void olig(bool isdna, int option, structure *ct, int length,double c, int **table,int **numofsubstructures,datatable& data,datatable& ddata, rddata *hybriddata, int Usesub,ProgressHandler *update,
		  thermo* helixstack,int start, int stop, siPREFILTER *prefilter,int foldsize,int distance,char *shapefile,int *TEST,bool WRITE=false); 

//=======================================================================
//this function reads data into structure rddata
int readrd (rddata* data,char* dnarna);

//=======================================================================
//this function writes a tab delimited (or HTML) file with the oligowalk data
	//Filename is the target file.
	//ct contains the sequnece
	//table is table filled by OligoWalk
	//numofsubstructure
	//length is the length of the oligonucleotides
	//isdna is whether the oligos are dna (true = DNA, false = RNA)
	//conc is the oligonucleotide concentration
	//Usesub is whether suboptimal structures were used
	//start and stop are indexes for the start and stop of the walk, respectively
	//prefilter is John Lu's siRNA prefilter
	//foldsize is the size of the folding region that is centered around the oligo
	//mask
	//asuf
	//tofe
	//fnnfe
	//isHTML is whether this is HTML or tab delimited text (false = tab delimited, true = HTML)
void report(const char* filename, structure *ct, int **table,int **numofsubstructures, const int length, const bool isdna,
			const double conc, const int Usesub,const int start,const int stop,const siPREFILTER *prefilter,const int foldsize,
			const bool *mask=NULL, const double asuf=0, const double tofe=0, const double fnnfe=0, const bool isHTML=true);
//Alternative format that is tab-delimited
void report(const char* filename, structure *ct, int **table, const int length, const bool isdna,
		const double conc, const int Usesub,const bool *mask, const double asuf, const double tofe, const double fnnfe);
//=======================================================================
//returns the base complementary to base i
//with 1 == A, 2 == C, 3 == G, 4 == U
int complement(int i, structure *ct);

//=======================================================================
//post-filter of functional siRNA written by Dave.
void filterbysirna(structure *ct, int **table, int length, datatable *data, 
				   bool *mask, double asuf, double tofe, double fnnfe);

//=======================================================================
//copy the arrays to be reused with different index when the folded region move one nucleotide to the right
inline void scancopy(OligoPclass *region, OligoPclass *copyregion);
//copy every arrays to be reused when the folded region did not move right yet 
inline void scancopyend(OligoPclass *region, OligoPclass *copyregion) ;

#endif

