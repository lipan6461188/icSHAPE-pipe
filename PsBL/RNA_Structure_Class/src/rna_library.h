
#if !defined (RNA_LIBRARY)
#define RNA_LIBRARY



#include <cmath>
#include <cstring>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <math.h>

#include "stackstruct.h"
#include "defines.h"
#include "common_utils.h"

// names of datatable files
#define DT_RNA "rna" // alphabet prefix for RNA -- i.e.  rna.specification.dat
#define DT_DNA "dna" // alphabet prefix for DNA -- i.e.  dna.specification.dat
#define EXT_ENERGY ".dg"  // suffix for energy files
#define EXT_ENTHALPY ".dh" // suffix for enthalpy files
#define EXT_SPEC ".dat" // suffix for specification files
#define F_spec "specification"
#define F_loop "loop"
#define F_stack "stack"
#define F_tstackh "tstackh"
#define F_tstacki "tstacki"
#define F_tloop "tloop"
#define F_miscloop "miscloop"
#define F_dangle "dangle"
#define F_int22 "int22"
#define F_int21 "int21"
#define F_triloop "triloop"
#define F_coax "coaxial"
#define F_tstackcoax "tstackcoax"
#define F_coaxstack "coaxstack"
#define F_tstack "tstack"
#define F_tstackm "tstackm"
#define F_int11 "int11"
#define F_hexaloop "hexaloop"
#define F_tstacki23 "tstacki23"
#define F_tstacki1n "tstacki1n"

#define DATATABLE_ERRMSG "The thermodynamic parameters could not be read. Please make sure the DATAPATH environment variable is set (e.g. to the RNAstructure/data_tables directory)."

using namespace std;


struct datatable {
    //this structure contains all the info read from thermodynamic data files


    //define types for vectors.
    typedef std::vector<std::vector<int> > intVector2D;
    typedef std::vector<std::vector<short int> > shortVector2D;
    typedef std::vector<std::vector<std::vector<std::vector<short int> > > > shortVector4D;
    typedef std::vector<std::vector<std::vector<std::vector<std::vector
            <std::vector<short int> > > > > > shortVector6D;
    typedef std::vector<std::vector<std::vector<std::vector<std::vector
            <std::vector<std::vector<short int> > > > > > > shortVector7D;
    typedef std::vector<std::vector<std::vector<std::vector<std::vector
            <std::vector<std::vector<std::vector<short int> > > > > > > > shortVector8D;

    // a constructor to get RT.
    datatable();
	
	////! Set whether or not to allow unknown nucleobases. If true, basetonum will return 0 instead of -1 for unknown bases.
	//void setAllowUnknownBases(const bool allow);
	////! Get whether or the datatable is in permissive mode. See setPermissive.
	//bool getAllowUnknownBases();
	//! Get the name of the alphabet that is currently loaded. This will be "" if no alphabet is loaded or if the name is unknown.
	string GetAlphabetName();

    float RT;
	string alphabetName; // the name of the currently loaded alphabet, e.g. "rna" or "dna" etc.
	bool allowUnknownBases;   // whether or not to allow unknown nucleobases. If true, basetonum will return 0 instead of -1 for unknown bases.
	bool loadedAlphabet;      // whether the alphabet was loaded (by a call to opendat)
	bool loadedTables;        // whether the thermodynamic parameters were loaded (by a call to opendat with the loadThermo parameter set to true)
    string data_directory;    // the directory from which the current parameters have been loaded. This is re-used in ScaleToTemperature etc.
	double temperature;       // the temperature corresponding to the currently loaded free energy parameters.

    // declare variables for specification.dat data. alphabet, paring, not_pairing.
    // declare alphabet vector, storing alphabets, 1st dimension size is the ALPHABET_SIZE.
    // pairing info 2D bool vector paring    
    // not_paring info 1D char vector single
    std::vector<std::vector<char> > alphabet; //alphabets (all bases)
    std::vector<std::vector<bool> > pairing;  //pairing information 
    std::vector<char> not_pairing; // Single bases
    std::vector<char> non_interacting; // non-interacting bases
    std::vector<char> linker; // linker
	std::vector<bool> LinkerInts; //indicates where an int representation of a nuc is the linker
    
    //variables were used in old tloop, triloop and hexaloop. New function does not need them.
    //short int numofhexaloops, numoftriloops, numoftloops;

    //declare variables for data read from miscloop.dat. 
    short int maxpen, auend, gubonus, cint, cslope, c3, efn2a, efn2b,
              efn2c, init, mlasym, strain, singlecbulge;
    float prelog;

	//Some variables to track the number of triloops, tetraloops, and hexaloops, used for backwards compatability
	short numofhexaloops, numoftloops, numoftriloops;

    //poppen and eparam will be signed with news size vector in the the constructor.
    std::vector<short int> poppen; 
    std::vector<short int> eparam;


    //declare variables for data read from loop.dat. 
    std::vector<short int> inter;
    std::vector<short int> bulge;
    std::vector<short int> hairpin;


    //declare variables for data read from (t,tri,hexa)loop.dg. 
    //tloop and triloop and hexaloop are built as int vector for larger number.
    intVector2D tloop; 
    intVector2D triloop; 
    intVector2D hexaloop;




    // //declare 4D matrices of short for dangle, stack, tstack, tstkh, tstki, coax, tstackcoax
    // //coaxstack, tstkm, tstki23, tstki1n
    shortVector4D dangle, stack, tstack, tstkh, tstki, coax, tstackcoax, 
                coaxstack, tstkm, tstki23, tstki1n;               


    // //declare vectors of short for iloop11, iloop21 and iloop22 data 
    //read from int11.dg, int21.dg, and int22.dg.           
    shortVector6D iloop11;
    shortVector7D iloop21;
    shortVector8D iloop22;
   
    //functions
	//given a sequence (i.e. the numseq member of the structure class), check if two nucleotides
	//are allowed to pair
	bool can_pair(int i, int j, short* sequence);
    //convert base to a int number, based on read in alphabet vector from specification.dat file.
    int basetonum(char base);
	//convert a number to a base
	char numtobase(int i);
    //convert a sequence to a int number, by an algorithme using aphabet size as the power factor.
    int seqtonum(std::string seq);
    //read in specification.dat to get alphabet vector info, paring bool vector info and not_paring char 1D vector.
    bool read_spec_file(char* spec_fname);
    //read in data table file to a string vector.
    bool read_datatable_lines(char* ifname, std::vector<std::string> &v_datatable_lines);
    //read in loop.dat type file, get inter, bulge and hairpin info
    bool read_loop(char* ifname, std::vector<short int> &inter,
                   std::vector<short int> &bulge, std::vector<short int> &hairpin);




    // //read in tloop.dat type files, get a 2D matrix info
    bool read_xloop(char* ifname, intVector2D &matrix);
    // //read in dangle.dat type files, get a 4D matrix info
    bool read_dangle(char* ifname, shortVector4D &matrix);



    //read in stack.dat type files, get a 4D matrix info
    bool read_4D_table(char* ifname, shortVector4D &matrix);



    // read in int11.dat type file, get a 6D matrix info
    bool read_6D_table(char* ifname, shortVector6D &matrix);
    //read in int21.dat type file, get a 7D matrix info
    bool read_7D_table(char* ifname, shortVector7D &matrix);
    //read in int22.dat type file, get a 8D matrix info
    bool read_8D_table(char* ifname, shortVector8D &matrix);



    // //read in miscloop.dat
    bool read_miscloop(char* ifname, float &prelog, short int &maxpen,
                   short int &efn2a, short int &efn2b, short int &efn2c, 
                   short int &strain, short int &auend, short int &gubonus, 
                   short int &cslope, short int &cint, short int &c3, short int &init,
                   short int &singlecbulge, std::vector<short int> &poppen,
                   std::vector<short int> &eparam);

	//When reading data tables from save files, the tables need to be allocated without reading the individiual data files.
	//This function is able to allocate all the arrays.  It requires that the alphabet size be determined.
	void allocate_data_tables();

	//Resize a 4D table
	void size4D(shortVector4D &matrix);
    
    //! Reset properties in preparation for (re)loading data tables.
    void reset();

	//Function opens data files to read thermodynamic data
	//The names of the files are generated from the directory (e.g. $DATAPATH) and the alphabet name,
	//along with the name of the data parameters and the appropriate extension (.dg for energies, .dh for enthalpies).
	// If skipThermoTables is true, only the alphabet is loaded. The enthalpy/entropy tables are NOT.
	// If skipThermoTables is false both the alphabet and thermo tables are fully loaded.
	int opendat(const char *directory, const char *alphabet,  const bool useEnthalpies=false, const bool skipThermoTables = false);

    //Function opens data files to read thermodynamic data
    //name for Getdat are changed, since they will have same name with data parameters.
    //rename rule(for files output only one type of parameters, file name is 'f'+ parameter name, 
    //exceptions are miscloop and loop, they generate multiple data variables.)
    //***********************
    // int opendat (char *fspec, char *floop, char *fdangle, char *fstack,  
    //         char *fcoax, char *fcoaxstack, char *ftstackcoax, char *ftstki,  
    //         char *ftstack, char *ftstkh, char *ftstkm, char *ftstki23,    
    //         char *ftstki1n, char *filoop11, char *filoop21, char *filoop22, 
    //         char *ftloop, char *ftriloop, char *fhexaloop, char *fmiscloop);
    int opendat(char *fspec, char *floop,  char *fstack, char *ftstkh, char *ftstki, 
		char *ftloop, char *fmiscloop, char *fdangle, char *filoop22, char *filoop21,
            char *fcoax, char *ftstackcoax, char *fcoaxstack,  
            char *ftstack,  char *ftstkm, char *ftriloop, char *filoop11, char *fhexaloop,  
            char *ftstki23, char *ftstki1n);
    //***********************

    //! Reads in the *.dh enthalpy tables and calls dG_T to calculate the free energy tables at the specified temperature.
    //! \param temperature The new temperature in degrees Kelvin
    //! \return Returns 0 on success. On failure, returns an error code compatible with RNA::GetErrorMessage
    int ScaleToTemperature(const double temperature);

    //! Returns the temperature at which the current free energy tables are valid.
    //! By default this is 37 C (310.15 K)
    //! But calling ScaleToTemperature(new_temp) causes the free energy parameters to be recalculated for a different temperature.
    inline double GetLoadedTemperature() const {
        return temperature;
    }

	//Return true if the base number or character is the intermolecular linker
	bool isLinker(const int baseNumber);
	bool isLinker(const char baseChar);

	private:
		void processAlphabet();   // perform final tasks after reading specification file
		void processDat(); // perform final tasks after reading data files
};

#define ERR_BAD_RESTRAINT_NUC_POS 2004  // corresponds to RNA::GetErrorMessage for "Invalid Nucleotide Number in Restraint file"
#define ERR_BAD_RESTRAINT_FILE 2002 // corresponds to RNA::GetErrorMessage for "Could not open or read Restraint file"
#define DEFAULT_RESTRAINT_VALUE -999.0 // default restraint value for nucleotides that are not listed explicitly in a restraints file (e.g. SHAPE)

//! ReadRestraints reads a SHAPE/DMS file
//!
//! /param valuesRead is a vector of doubles that is filled with the restraint data. 
//!                   This vector is 0-indexed, i.e. the value for base 1 is stored at valuesRead[0].
//!                   The vector be empty or have its size set before the call (if the number of nucleotides are known).
//! /return Returns ERR_BAD_RESTRAINT_NUC_POS if a bad position is listed in the file or ERR_BAD_RESTRAINT_FILE if the file cannot be read. Otherwise returns 0 on success.
int ReadRestraints(vector<double> &valuesRead, const char* SHAPEFile);
// Write restraints (e.g. SHAPE, DMS etc) to a file
//! /param values is a vector of doubles that is filled with the restraint data. 
//!               This vector must be 0-indexed, i.e. the value for base 1 is stored at values[0].
//! /param outfile A string specifying the path to the output file.
//! /param append Whether to append data to a file if it exists (the default is to truncate existing files before writing).
//! /return Returns ERR_BAD_RESTRAINT_FILE if the file cannot be opened for writing. Otherwise returns 0 on success.
int WriteRestraints(const vector<double> &values, const string& outfile, const bool append = false);

//write is used to write data to a save file
void write(std::ofstream *out,short *i);
void write(std::ofstream *out,bool *i);
void write(std::ofstream *out,int *i);
void write(std::ofstream *out,char *i);
void write(std::ofstream *out,string *i);
void write(std::ofstream *out,float *i);
void write(std::ofstream *out,double *i);
void writesinglechar(std::ofstream *out,char *i);
template<typename T>
void write(ofstream* out, vector<T>* vec){
	int i = vec->size();
	write(out, &i);
	for(typename vector<T>::iterator it=vec->begin(); it!=vec->end(); ++it){
		T val = *it;
		write(out, &val);
	}
}
template<>
inline void write(ofstream* out, vector<char>* vec){
	int i = vec->size();
	write(out, &i);
	for(vector<char>::iterator it=vec->begin(); it!=vec->end(); ++it){
		char val = *it;
		writesinglechar(out, &val);
	}
}
void write(std::ofstream* out, datatable* data);


//read is used to read data from a save file
void read(std::ifstream *out,short *i);
void read(std::ifstream *out,bool *i);
void read(std::ifstream *out,int *i);
void read(std::ifstream *out,char *i);
void read(std::ifstream *out,string *i);
void read(std::ifstream *out,float *i);
void read(std::ifstream *out,double *i);
void readsinglechar(std::ifstream *out,char *i);
template<typename T> T readSingle(std::ifstream *out);

template<typename T>
void read(ifstream* out, vector<T>* v){
	int sz;
	read(out, &sz);
	v->resize(sz);
	for(typename vector<T>::iterator it=v->begin(); it!=v->end(); ++it){
		T val;
		read(out, &val);
		*it = val;
	}
}
template<>
inline void read(ifstream* out, vector<char>* v){
	int sz;
	read(out, &sz);
	v->resize(sz);
	for(vector<char>::iterator it=v->begin(); it!=v->end(); ++it){
		char val;
		readsinglechar(out, &val);
		*it = val;
	}
}
void read(std::ifstream* out, datatable* data);

void push(stackstruct *stack,int a,int b,int c,int d);//push info onto the stack
void pull(stackstruct *stack,int *i,int *j,int *open,int *null,int *stz);
														//pull info from the stack




void dG_T (float T, datatable &data, datatable &dhdata, datatable &dg);
	//Determine the free energy parameters at temperature T and store in dg 
		//from the free energy parameters at 37 degrees C (data) and the
		//temperature independent enathalpies (dhdata).

//Given the absolute temperature, T, return the free energy for a paramater as determined from the 
//free energy at 37 degrees C (dG) and the enthalpy change (assumed to be temperature independent).
inline short Tscale(float T,short dG, short dH)
{
	if (dG==INFINITE_ENERGY) return INFINITE_ENERGY; //keep infinity==infinity
	return dH - (short)floor((float)(dH-dG)*T/T37inK+0.5);
}



//short int tonumi(char *base); //converts base to a numeric

integersize penalty2(int i, int j, datatable *data);

//!Get the directory location of thermodynamic parameter files.
//! In order of priority, this will be:
//!  1. Any non-NULL value set by a prior call to setDataPath
//!  2. The location pointed to by the DATAPATH environment variable, (if it is valid)
//!  3. Any of several standard locations (e.g. ./data_tables or ../data_tables etc, if any are valid)
//!  4. The value of the defaultPath argument, which defaults to "." if not specified.
//! /param defaultPath The path to return if DATAPATH is empty and the directory cannot be auto-detected.
//! /param alphabetName If this is not NULL, it will be used to validate any auto-detected directory. 
//!                     If it is NULL, the "rna" or "dna" alphabets will be used.
const char* getDataPath(const char* const defaultPath = DATAPATH_DEFAULT, const char* const alphabetName=NULL);
//Set the location where programs should look for the thermodynamic parameter files.
//Programs that call getDataPath will receive the new location.
//(This also sets the DATAPATH environment variable for backwards compatibility.)
//This stores a copy of the path and returns a pointer to that copy.
const char* setDataPath(const char* const path);
//Show a warning that DATAPATH was not set and the program is using an auto-detected or default directory.
void warnDataPath(const char* const path, const bool autoDetected=true);
#endif
