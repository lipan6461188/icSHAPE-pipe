/*=======================================================================
intermolecular.h and intermolecular.cpp inculde funcions calculating and report
different free energy for binding in OligoWalk.
intermolecular_test.cpp fold the whole sequence, save partition function in a file 
and reuse it.
stochastic sampling method were also embedded.

They are revised based on Mathews' code from RNAStructure.
olig() generate the energy data;
report save the generated data;
siprefileter and filterbysirna() were used to prefilter and postfileter the 
functional siRNA respectively, using criterias other than thermodynamics


															----Aug. 6, 2006
															John(Zhi Lu)



Changes made by DHM: 7/12/09:
Removed folding when break local structure is used.

=======================================================================*/
//#define debugmode true
#undef debugmode
#define TESTMODE false   //read the save file instead of calculating if in TEST mode
#include "intermolecular.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>

//define the parameters for alltrace()
#define percent 100	//100% of the suboptimal structure windows
#define delta 6	//only the suboptimal structures having energy less than 2.8kcal/mol higher than optimal one
                    //not store structures in ct file for alltrace() to save memory
// #define MaxStructure 10000000 (10M) is defined in alltrace.h



void dgetdat(char *loop, char *stackf, char *tstackh, char *tstacki,
		char *tloop, char *miscloop, char *danglef, char *int22,
      char *int21,char *coax, char *tstackcoax,
      char *coaxstack, char *tstack, char *tstackm, char *triloop,
      char *int11, char *path);

int readrd (rddata* data,char* dnarna);

inline void scancopy(OligoPclass *region, OligoPclass *copyregion);
inline void scancopyend(OligoPclass *region, OligoPclass *copyregion) ;


/*=======================================================================
oligo fills the array table with thermodynamic data for each oligo (2nd dimension)
in the second dimension:
table[i][0] = overall DG
table[i][1] = duplex DG
table[i][2] = free energy of breaking target structure
table[i][3] = intramolecular oligo free energy
table[i][4] = intermolecular oligo free energy
table[i][5] = Tm of duplex (x10)

option 1 - break local target structure to bind oligo 
option 2 - refold target RNA after oligo binding
option 3 - no target structure considered

Usesub 0 - only consider optimal structure
Usesub 1 - like choice 3,using suboptimal structures,but the whole set from alltarce() function prediction
Usesub 2 - using partition funcion considering every possible structure  of target
			   suboptimal 2 can only used with mold 2
Usesub 3 - using suboptimal structures (heuristic method) for both oligo-free and oligo-bound target RNA
Usesub 4 - using stochasit sampling method to sample 1000 structures

prefileter 1 - using creteria to prefill functional siRNA
foldsize >0  - only folding a fragment with size=foldsize+binding length, 
			   which is centered on the siRNA binding region
			   when foldsize>1, only option 2 plus Usesub 2 is the availabe option

distance >0 - limit the maximum distance between nucleotides that can pair
shapefile - specify a SHAPE datafile, set shapefile[0] to '\0' not use SHAPE data

TESTS - run tests, set to -1 for no tests.
write - write sav files to save time in test mode

=======================================================================*/
void olig(bool isdna, int option, structure *ct, int length,double c, int **table,int **numofsubstructures,
		  datatable& data,datatable& ddata, rddata *hybriddata, int Usesub,ProgressHandler *update,
		  thermo* helixstack,int start, int stop, siPREFILTER *prefilter,int foldsize,int distance,char *shapefile,int *TEST,bool WRITE) {
	
	int i,j,k;
	int ip,jp;
	int foldstart, foldstop;
	int dh,ds,dgeff;
	int *energyarray, *ctenergy;//ctenergy is the energy array to store structure energies from alltrace()
	int numofstructures;
	long double k1b,k1u;
	long double fn,sn,sum;
	char savefile[250],pfsfile[250],pos[250];
	int energy=INFINITE_ENERGY;//to store the free energy of lowest free energy structure
	int *temp,**temp2;//temp will store basepairing info in option == 1
	

	PFPRECISION Q,Qc,Qolig,pftemp=310.15;//store partion function without and with constrains
	PFPRECISION rescaleinrefill;//rescaling facor for Qc relative to Q
	pfdatatable *pfdata,*dpfdata;  //store the data tables for partionfunction
	OligoPclass *interoligo,*intraoligo;
	OligoPclass *target,*targetcopy,*targettemp;
	structure *oligo,*oligo1;//to store the structural info for the oligo
						   //oligo is intermolecular structure; oligo1 is intramolecular
	structure *fracct;// to store folded fraction of target centered at siRNA binding site
	

	//Infrastructure for 
	//define shape information for the whole sequnence (ct)
	if ( shapefile[0] != '\0'  ) {
		ct->SHAPEslope = 32;
		ct->SHAPEintercept = -10;
		ct->ReadSHAPE(shapefile) ;
		if (Usesub==2) pfshape(ct,pftemp);
	}


	
	//----------------------------------------------------------------------------------------------
	//allocate space in oligo for the sequence
	//difine for inter oligo sequence, oligo structure can be used both by intra and inter molecular
	//oligo:  intermolecular structure   aagucXXXggcaa
	//oligo1: intramolecular structure   aaguc
	oligo = new structure();
	oligo1 = new structure();

	oligo->allocate(2*length+3);
	oligo1->allocate(length);

	string label;
	label = "Oligo_inter";

	oligo->SetSequenceLabel(label);

	label = "Oligo_intra";

	oligo1->SetSequenceLabel(label);

	
	for (j=1;j<=3;j++) {		oligo->inter[j-1] = length + j;	}

	foldsize= foldsize/2*2; // a trick to make foldsize even
	//----------------------------------------------------------------------------------------------
	//define the base pair distance constraint
	if (distance > 0) {
		ct->SetPairingDistance(distance);
		
	}
	//define the size of fracct, which is the region to be folded
	if (foldsize >0) {
		fracct = new structure();
		fracct->allocate(foldsize+length);
		label="fraction_of_target";
		fracct->SetSequenceLabel(label);
		
		//define the base pair distance constraint
		if (distance > 0) {
			fracct->SetPairingDistance(distance);
		}
	}	
	
	//----------------------------------------------------------------------------------------------
	//some variables declared for each option 
	if (option ==1) {
		if (Usesub==0)	{
			temp = new int [length];//allocate an array for storing base pairing information
			//if (foldsize == 0) {	
			ct->RemoveConstraints();
				//dynamic(ct,&data,1,10,0,0);
			efn2(&data,ct,1);
		   	energy = ct->GetEnergy(1);
			//}
		}
		else if (Usesub==3) {
			//if (foldsize == 0) {	
			ct->RemoveConstraints();
				//dynamic(ct,&data,1000,10,0);
			efn2(&data,ct);
			energyarray = new int [ct->GetNumberofStructures()+1];
			for (k=1;k<=ct->GetNumberofStructures();k++) energyarray[k] = ct->GetEnergy(k);
			temp2 = new int *[ct->GetNumberofStructures()+1];
			for (i=1;i<=ct->GetNumberofStructures();i++)   		temp2[i] = new int [ct->GetSequenceLength()];
			numofstructures=ct->GetNumberofStructures();
			//}
      	}
		
	}
	else if (option==2) {
   		temp = new int [ct->GetSequenceLength()+1];
		for (j=1;j<=ct->GetSequenceLength();++j) {
      		temp[j] = ct->GetPair(j);
		}
		if (Usesub ==0 ) {
			//ct->GetNumberofStructures() = 1;//only interested in lowest free energy structure
			for (int last = ct->GetNumberofStructures();last>1;--last) ct->RemoveLastStructure();

			if (foldsize == 0) {
				ct->RemoveConstraints();
				//ct->nnopair=0;
				strcpy(savefile, ct->GetSequenceLabel().c_str());
				if (savefile[strlen(savefile)-1]=='\n') savefile[strlen(savefile)-1]='\0'; //get off the new line charactor in the string
				strcat(savefile,"_s0.sav");
				//check if a sav file of partition function result exist with a C i/o function
				std::ifstream sav(savefile,std::ios::binary);
				if(TESTMODE && sav) {
					//read information from file if it was found
					read(&sav,&energy);
					sav.close();
				}
				else {
					//close the readonly file sav
					sav.close();
					//write the save file information so that the fold need not to be done again
					dynamic(ct,&data,1000,10,0,0,true);
				
					energy = ct->GetEnergy(1);
					if (WRITE) {
						std::ofstream sav(savefile,std::ios::binary);
						write(&sav,&energy);
						sav.close();
					}
				}
			}
		}
		else if	(Usesub==1)	 {	
			ctenergy=new int [MaxStructure+1];
			if (foldsize ==0) {	
				ct->RemoveConstraints();
				//ct->nnopair=0;
				strcpy(savefile, ct->GetSequenceLabel().c_str());
				if (savefile[strlen(savefile)-1]=='\n') savefile[strlen(savefile)-1]='\0'; //get off the new line charactor in the string
				strcat(savefile,"_s1.sav");
				//check if a sav file of partition function result exist with a C i/o function
				std::ifstream sav(savefile,std::ios::binary);
				if(TESTMODE && sav) {
					//read information from file if it was found
					read(&sav,&numofstructures);
					energyarray = new int [numofstructures+1];
					for (j=1;j<=numofstructures;j++)	read(&sav, energyarray+j);
					sav.close();
				}
				else {
					//close the readonly file sav
					sav.close();
					//calculate the whole set of suboptimal structures
					alltrace(false,ctenergy,ct,&data, percent, delta, NULL,NULL);
					//write the save file information so that the fold need not to be done again
					numofstructures=ct->GetNumberofStructures();
					energyarray = new int [numofstructures+1];
					for (j=1;j<=numofstructures;j++)	{energyarray[j]=ctenergy[j];}
					if (WRITE) {
 						std::ofstream sav(savefile,std::ios::binary);
						write(&sav,&numofstructures);
						for (j=1;j<=numofstructures;j++)		write(&sav, &energyarray[j]);
						sav.close();
					}
				}

			}
		}
		else if	(Usesub==3)	 {	
			if (foldsize ==0) {	
				//ct->nnopair=0;
				ct->RemoveConstraints();
				strcpy(savefile, ct->GetSequenceLabel().c_str());
				if (savefile[strlen(savefile)-1]=='\n') savefile[strlen(savefile)-1]='\0'; //get off the new line charactor in the string
				strcat(savefile,"_s3.sav");
				//check if a sav file of partition function result exist with a C i/o function
				std::ifstream sav(savefile,std::ios::binary);
				if(TESTMODE && sav ) {
					//read information from file if it was found
					read(&sav,&numofstructures);
					energyarray = new int [numofstructures+1];
					for (j=1;j<=numofstructures;j++)	read(&sav, energyarray+j);
					sav.close();
				}
				else {
					//close the readonly file sav
					sav.close();
					dynamic(ct,&data,1000,10,0);
					//write the save file information so that the fold need not to be done again
					numofstructures=ct->GetNumberofStructures();
					energyarray = new int [numofstructures+1];
					for (j=1;j<=numofstructures;j++)	{energyarray[j]=ct->GetEnergy(j);}
					if (WRITE) {
 						std::ofstream sav(savefile,std::ios::binary);
						write(&sav,&numofstructures);
						for (j=1;j<=numofstructures;j++)		write(&sav, &energyarray[j]);
						sav.close();
					}
				}

			}
		}
		else if (Usesub==2) {
			if (isdna) {//oligo is a DNA
				dpfdata = new pfdatatable (&ddata,scalingdefinition,pftemp);
				pfdata = new pfdatatable (&data,scalingdefinition,pftemp);
				intraoligo = new OligoPclass(oligo1,dpfdata);
				interoligo = new OligoPclass(oligo,dpfdata);
			}
			else {//oligo is a RNA
				pfdata = new pfdatatable (&data,scalingdefinition,pftemp);
				intraoligo = new OligoPclass(oligo1,pfdata);
				interoligo = new OligoPclass(oligo,pfdata);
			}
			//calculate partion function for the whole target without constrain
			//ct->nnopair=0;
			ct->RemoveConstraints();
			if (foldsize==0) {//folding the whole sequence at one time
				target = new OligoPclass(ct,pfdata);
				strcpy(savefile, ct->GetSequenceLabel().c_str());
				if (savefile[strlen(savefile)-1]=='\n') savefile[strlen(savefile)-1]='\0'; //get off the new line charactor in the string
				//strcpy(savefile, ct->ctlabel[1]);
				//savefile[strlen(savefile)-1]='\0'; //get off the new line charactor in the string
				strcat(savefile,"_s2.sav");
				strcpy(pfsfile, ct->GetSequenceLabel().c_str());
				if (pfsfile[strlen(savefile)-1]=='\n') pfsfile[strlen(savefile)-1]='\0'; //get off the new line charactor in the string
				strcat(pfsfile,"_s2.pfs");
				//check if a sav file of partition function result exist with a C i/o function
				std::ifstream sav(savefile,std::ios::binary);
				std::ifstream pfs(savefile,std::ios::binary);
				if (TESTMODE && sav && pfs  ){
					pfs.close();
					//read information from file if it was found
					read(&sav,&Q);
					for (i=0;i<=ct->GetSequenceLength();++i) {
						read(&sav,&(target->copyw5[i]));
						for (j=0;j<=ct->GetSequenceLength();++j) {
							read(&sav,&(target->copyv->f(i,j)));
							read(&sav,&(target->copyw->f(i,j)));
							read(&sav,&(target->copywca[i][j]));
							read(&sav,&(target->copywmb->f(i,j)));
							read(&sav,&(target->copywl->f(i,j)));
							read(&sav,&(target->copywmbl->f(i,j)));
							read(&sav,&(target->copywcoax->f(i,j)));
						}	
					}
					sav.close();
				}
				else {
					//close the readonly file sav
					sav.close();
					pfs.close();
					//calculate the partition function if no file found
					if (WRITE) {
						target->partition4refill(&Q,pfsfile);
					}
					else target->partition4refill(&Q);

				
					//write the save file information so that the partition function can be re-folded,
					if (WRITE) {
						std::ofstream sav(savefile,std::ios::binary);
						write(&sav,&Q);
						for (i=0;i<=ct->GetSequenceLength();++i) {
							write(&sav,&(target->copyw5[i]));
							for (j=0;j<=ct->GetSequenceLength();++j) {
								write(&sav,&(target->copyv->f(i,j)));
								write(&sav,&(target->copyw->f(i,j)));
								write(&sav,&(target->copywca[i][j]));
								write(&sav,&(target->copywmb->f(i,j)));
								write(&sav,&(target->copywl->f(i,j)));
								write(&sav,&(target->copywmbl->f(i,j)));
								write(&sav,&(target->copywcoax->f(i,j)));
							}		
						}
						sav.close();
					}
						
				}
			}
			else if(prefilter->useit!=0) {//using prefilter and fold different region each time
				target = new OligoPclass(fracct,pfdata);

			}
		}
		else if (Usesub==4) {
			if (isdna) {//oligo is a DNA
				dpfdata = new pfdatatable (&ddata,scalingdefinition,pftemp);
				pfdata = new pfdatatable (&data,scalingdefinition,pftemp);
				intraoligo = new OligoPclass(oligo1,dpfdata);
				interoligo = new OligoPclass(oligo,dpfdata);
			}
			else {//oligo is a RNA
				pfdata = new pfdatatable (&data,scalingdefinition,pftemp);
				intraoligo = new OligoPclass(oligo1,pfdata);
				interoligo = new OligoPclass(oligo,pfdata);
			}
			//calculate partion function for the whole target without constrain
			ct->RemoveConstraints();
			if (foldsize==0) {//folding the whole sequence at one time
				target = new OligoPclass(ct,pfdata);
				strcpy(pfsfile, ct->GetSequenceLabel().c_str());
				if (pfsfile[strlen(savefile)-1]=='\n') pfsfile[strlen(savefile)-1]='\0'; //get off the new line charactor in the string
				strcat(pfsfile,"_s2.pfs");
				//check if a sav file of partition function result exist with a C i/o function
				std::ifstream pfs(pfsfile,std::ios::binary);
				if(TESTMODE && pfs) {//read information from file if it was found
					pfs.close();				
				}
				else{
					//close the readonly file sav
					pfs.close();
					//calculate the partition function if no file found
					target->partition4refill(&Q,pfsfile); 
					
				}
				stochastic(ct,pfsfile,1000,5);
				efn2(&data,ct);
				//ofstream sout("RL_sample_oligo->out");
				//for (j=1;j<=1000;j++)	{ sout << ct->energy[j]<<"\n";	}
				//sout.close();
				numofstructures=ct->GetNumberofStructures();
				energyarray = new int [numofstructures+1];
				for (j=1;j<=numofstructures;j++)	{energyarray[j]=ct->GetEnergy(j);}
			}
			else if(prefilter->useit!=0) {//using prefilter and fold different region each time
				//allocate target only, if prefilter==0, allocate both target and targetcopy later
				target = new OligoPclass(fracct,pfdata);

			}
		}
		
	}
	//having single-strand constrained for the calculation of constrained energy
	//ct->nnopair=length;
	//ct->checknopair();	

	//------------------------------------------------------------------------------------------
	//------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------
	//begin to scan the target from start to stop
	for (i=start; i<=stop; ++i) {

		
		//communicate progress
		if (update!=NULL) {
   			update->update (int((double (i-start))*100/(double (stop-start))));
		}


		//Remove previous constraints that were imposed:
		ct->RemoveConstraints();

		//define the oligo sequence
		//oligo1->intermolecular = false;
		//oligo1->GetSequenceLength() = length;
		//set sequence
		 for (j=1;j<=length;j++) {
			oligo1->numseq[j] = complement(i+length-j,ct);
			oligo->numseq[j] = complement(i+length-j,ct);
		}
   		//-------------------------------------------------------------------------------------------
		//prefiltering the functional siRNA
		if (prefilter->useit != 0) {		
			prefilter->count(oligo1,i,TEST[i]);
			if (prefilter->score[i] < FILTER_PASS) {
				continue;
			}
		}
	
		//-------------------------------------------------------------------------------------------
		//-------------------------------------------------------------------------------------------
		//calculate the free energy of breaking target structure

		//Option 1: break local structure only
	    
		if (option==1) {
			if (foldsize ==0) {
				//-------------------------------------------------------------------------------------------
				//option 1 + consider only the first structure:
				   if (Usesub ==0){
					//store basepairing info
		     		for (j=0;j<length;j++) {
						temp[j] = ct->GetPair(i+j);
		        		if (ct->GetPair(i+j) > 0) {
							ct->RemovePair(i+j);
		        			
			       		}
		     		}
	         	
					efn2(&data,ct);
		    		table[i][2] = ct->GetEnergy(1) - energy;
	         	
					//restore basepairing
		    		for (j=0;j<length;j++) {
		      			if (temp[j]>0) {
							ct->SetPair(i+j,temp[j]);
		  
		        		}
		    		}
				}	
				//-------------------------------------------------------------------------------------------
				//option 1 + consider all subotimal structures:
				else if (Usesub==3) { 
					sn = 0;
					sum = 0;
				  //store the basepairing info and force the region of hybridization
				 //to have no structure
				  for (k=1;k<=ct->GetNumberofStructures();++k) {
		          		for (j=0;j<length;j++) {//store basepairing info
		     				temp2[k][j] = ct->GetPair(i+j,k);
		           			if (ct->GetPair(i+j,k) != 0) {
								ct->RemovePair(i+j,k);
		         				
		          			}
		      			}

					}
					efn2(&data,ct);
				  for (k=1;k<=ct->GetNumberofStructures();k++) {

		          		fn = exp(-(((long double)energyarray[k]-energyarray[1])/(RT_37C*conversionfactor)));
					   sn = sn + fn;
		   				sum = sum+fn*(((long double)(ct->GetEnergy(k)- energyarray[k])));
				   }
					//restore the basepairing
					for (k=1;k<=ct->GetNumberofStructures();k++) {
			       		for (j=0;j<length;j++) {
							if (temp2[k][j]>0) {
								ct->SetPair(i+j,temp2[k][j],k);
           						
           					}
        				}
					}
	         		
					table[i][2] = (int)(sum/sn);
					numofsubstructures[i][0]=numofstructures; //report in report() function
					numofsubstructures[i][1]=ct->GetNumberofStructures(); 

				}
			}
			if (foldsize > 0 ) {
				//define the region for refolding 
				if ( (i-foldsize/2) > 1 && (i+length-1+foldsize/2) < ct->GetSequenceLength() ) {
					for (j=1;j<=foldsize+length;j++) {
      					fracct->numseq[j] = ct->numseq[i-foldsize/2+j-1];
					}
				}
				else if ( (i-foldsize/2)<=1 ) 
				for (j=1;j<=foldsize+length;j++) {
      				fracct->numseq[j] = ct->numseq[j];
				}
				else if( (i-1+length+foldsize/2)>=ct->GetSequenceLength() )
				for (j=1;j<=foldsize+length;j++) {
      				fracct->numseq[j] = ct->numseq[(ct->GetSequenceLength())+j-foldsize-length];
				}
				fracct->RemoveConstraints();
				//fracct->nnopair=0;
				//calculate the strucutre without constrains
				if (Usesub==0)	{
					//Remove structures, after the first because only need 1 structure
					//for (int structurenum=fracct->GetNumberofStructures();structurenum>1;--structurenum) {
					//	fracct->RemoveLastStructure();

					//}
					
					dynamic(fracct,&data,1,10,0,0);
					efn2(&data,fracct);
		   			energy = fracct->GetEnergy(1);
				}
				else if (Usesub==3) {
					dynamic(fracct,&data,1000,10,0);
					efn2(&data,fracct);
					energyarray = new int [fracct->GetNumberofStructures()+1];
					for (k=1;k<=fracct->GetNumberofStructures();k++) energyarray[k] = fracct->GetEnergy(k);
					temp2 = new int *[fracct->GetNumberofStructures()+1];
					for (k=1;k<=fracct->GetNumberofStructures();k++)   		temp2[k] = new int [fracct->GetSequenceLength()];
					numofstructures=fracct->GetNumberofStructures();
				
      			}
				//set the constrained position
				if (i-foldsize/2<=1) {//constrained positions are different at two ends, as the folded region did not move
					foldstart=i;
					foldstop=i+length-1;
				} 
				else if(i+length-1+foldsize/2>= (ct->GetSequenceLength()) ) {
					foldstart=foldsize+length-(ct->GetSequenceLength()) +i;
					foldstop=foldsize+length+length-1-(ct->GetSequenceLength()) +i;
				}
				else {//folded region begin to move, constrained position will be always in the middle of this region
					foldstart=foldsize/2+1;
					foldstop=foldsize/2+length;
				}
				
				//recalc the energy again
				if (Usesub==0) {
					//store basepairing info
		     		for (j=foldstart;j<=foldstop;j++) {
						temp[j-foldstart] = fracct->GetPair(j);
		        		if (fracct->GetPair(j) > 0) {
							fracct->RemovePair(j);
			       		}
		     		}
	         	
					efn2(&data,fracct);
		    		table[i][2] = fracct->GetEnergy(1) - energy;
	         	
					//restore basepairing
		    		for (j=foldstart;j<=foldstop;j++) {
		      			if (temp[j-foldstart]>0) {
							fracct->SetPair(j,temp[j-foldstart]);
		        		}
		    		}
				}
				else if (Usesub==3) {
					sn = 0;
					sum = 0;
					//store the basepairing info and force the region of hybridization
					//to have no structure
					for (k=1;k<=fracct->GetNumberofStructures();k++) {
		          		for (j=foldstart;j<=foldstop;j++) {//store basepairing info
		     				temp2[k][j] = fracct->GetPair(j,k);
		           			if (fracct->GetPair(j,k) != 0) {
								fracct->RemovePair(j,k);
		         				
		          			}
		      			}

					}
					efn2(&data,fracct);
					for (k=1;k<=fracct->GetNumberofStructures();k++) {

		          		fn = exp(-(((long double)energyarray[k]-energyarray[1])/(RT_37C*conversionfactor)));
						sn = sn + fn;
		   				sum = sum+fn*(((long double)(fracct->GetEnergy(k)- energyarray[k])));
					}
					//restore the basepairing
					for (k=1;k<=fracct->GetNumberofStructures();k++) {
			       		for (j=foldstart;j<=foldstop;j++) {
							if (temp2[k][j]>0) {
								fracct->SetPair(j,temp2[k][j],k);
           						
           					}
        				}
					}
				
					table[i][2] = (int)(sum/sn);
					numofsubstructures[i][0]=numofstructures; //report in report() function
					numofsubstructures[i][1]=fracct->GetNumberofStructures(); 

					delete[] energyarray;
					for (k=1;k<=fracct->GetNumberofStructures();k++) {	delete[] temp2[k];	}
					delete[] temp2;

				}
			}
		}//end option 1


		//-------------------------------------------------------------------------------------------
		//option 2: refold  rna
		else if (option == 2) {
			//not scan: refolding the whole sequence
			if (foldsize ==0) {
				//reset the nopair constrains
				ct->RemoveConstraints();
				for (j=0;j<length;j++)		{\
					ct->AddSingle(i+j);
				}		
				
				//option 2 +notrefold+ consider only the first suboptimal structure
				if (Usesub ==0) {
					dynamic(ct,&data,1000,10,0,0,true);
				
					table[i][2] = ct->GetEnergy(1) - energy;
				}
				//option 2 +notrefold+ consider all suboptimal structures with ensemble energy
				else if (Usesub==1) {
					//sum of unconstrained energy from energyarray[]
					Q = (PFPRECISION) 0;
					for (k=1;k<=numofstructures;k++) {
	          			fn = exp(-(((long double)energyarray[k]-energyarray[1])/(RT_37C*conversionfactor)));
						Q = Q + fn;
         			}
									
					//calculate the whole set of suboptimal structures with constrain
					alltrace(false,ctenergy,ct,&data, percent, delta, NULL,NULL);
					//sum of constrained energy now
					Qc = (PFPRECISION) 0;
					for (k=1;k<=ct->GetNumberofStructures();k++) {
	          			fn = exp(-(((long double)ctenergy[k]-energyarray[1])/(RT_37C*conversionfactor)));
						Qc = Qc + fn;
         			}
					//final free energy difference					
					table[i][2] = (int)( conversionfactor*RT_37C*log( (long double)Q/(long double)Qc ) );

					numofsubstructures[i][0]=numofstructures; //report in report() function
					numofsubstructures[i][1]=ct->GetNumberofStructures(); 
				}
				//option 2 +notrefold+ consider heuristic suboptimal structures with average energy
				else if (Usesub==3) {
					sum = 0;
					sn = 0;
					for (k=1;k<=numofstructures;k++) {
	          			fn = exp(-(((long double)energyarray[k]-energyarray[1])/(RT_37C*conversionfactor)));
						sn = sn + fn;
         				sum = sum + fn*((long double)energyarray[k]);
					}
					table[i][2] = (int)(sum/sn);
								
					//save the energies in a save file
					strcpy(savefile, ct->GetCtLabel(1).c_str());
					savefile[strlen(savefile)-1]='\0'; 
					sprintf(pos,"_%d",i);
					strcat(savefile,pos);
					strcat(savefile,"_s3_0_constrain.sav");
					//fold the constrained structure
					dynamic(ct,&data,1000,10,0);
					if (WRITE) {
						std::ofstream sav(savefile,std::ios::binary);
						int localint = ct->GetNumberofStructures();
						write(&sav,&(localint));
						for (j=1;j<=ct->GetNumberofStructures();j++)		{
							localint = ct->GetEnergy(j);
							write(&sav, & (localint));
						}
						sav.close();
					}
					sum = 0;
					sn = 0;
					for (k=1;k<=ct->GetNumberofStructures();k++) {
	          			fn = exp(-(((long double)ct->GetEnergy(k)-ct->GetEnergy(1))/(RT_37C*conversionfactor)));
						sn = sn + fn;
         				sum = sum + fn*((long double)ct->GetEnergy(k));
						//cout<< ct->energy[k]<<"\n";
					}
					table[i][2] = (int)(sum/sn) - table[i][2] ;
					numofsubstructures[i][0]=numofstructures; //report in report() function
					numofsubstructures[i][1]=ct->GetNumberofStructures(); 
		   		}
				//option 2 +notrefold+ partionfunction
				else if (Usesub ==2) {
					//save the energies in a save file
					strcpy(pfsfile, ct->GetCtLabel(1).c_str());
					pfsfile[strlen(pfsfile)-1]='\0'; 
					sprintf(pos,"_%d",i);
					strcat(pfsfile,pos);
					strcat(pfsfile,"_s2_0_constrain.pfs");
					if (WRITE) {
						target->refill(ct,&Qc,i, i+length-1,rescaleinrefill,pfsfile);
					}
					else target->refill(ct,&Qc,i, i+length-1,rescaleinrefill);


					table[i][2] = (int)(conversionfactor*RT_37C* 
						( (long double)(ct->GetSequenceLength())*log(rescaleinrefill) 
						+ log((long double)Q/(long double)Qc) ) );
					
				}
				//option 2 +notrefold+ stochastic
				else if (Usesub ==4) {
					sum = 0;
					for (k=1;k<=numofstructures;k++) {
	          			sum = sum + energyarray[k];
					}
					table[i][2] = (int)(sum/numofstructures);
				
					strcpy(pfsfile, ct->GetCtLabel(1).c_str());
					pfsfile[strlen(pfsfile)-1]='\0'; 
					sprintf(pos,"_%d",i);
					strcat(pfsfile,pos);
					strcat(pfsfile,"_s2_0_constrain.pfs");
					std::ifstream pfs(pfsfile,std::ios::binary);
					if(TESTMODE && pfs) {
						pfs.close();
					}
					else{
						pfs.close();
						target->refill(ct,&Qc,i, i+length-1,rescaleinrefill,pfsfile);
						
					}
					stochastic(ct,pfsfile,1000,5);
					efn2(&data,ct);
					sum = 0;
					for (k=1;k<=ct->GetNumberofStructures();k++) {
	          			sum = sum + ct->GetEnergy(k);
					}
					table[i][2] = (int)(sum/ct->GetNumberofStructures()) - table[i][2] ;
					numofsubstructures[i][0]=numofstructures; //report in report() function
					numofsubstructures[i][1]=ct->GetNumberofStructures(); 			
				}
				
			}
			//----------------------------------------------------------
			//----------------------------------------------------------
			//scan: only folding the region close to the binding site
			else if (foldsize >0) {
				
				//define the region for refolding 
				if ( (i-foldsize/2) > 1 && (i+length-1+foldsize/2) < ct->GetSequenceLength() ) {
					for (j=1;j<=foldsize+length;j++) {
      					fracct->numseq[j] = ct->numseq[i-foldsize/2+j-1];
					}
				}
				else if ( (i-foldsize/2)<=1 ) 
				for (j=1;j<=foldsize+length;j++) {
      				fracct->numseq[j] = ct->numseq[j];
				}
				else if( (i-1+length+foldsize/2)>=ct->GetSequenceLength() )
				for (j=1;j<=foldsize+length;j++) {
      				fracct->numseq[j] = ct->numseq[(ct->GetSequenceLength())+j-foldsize-length];
				}
				fracct->RemoveConstraints();

				//----------------------------------------------------------
				//fold the scanned region without any constrain:
				//refold for different Usesub options
				if (Usesub ==0) {
					dynamic(fracct,&data,1,10,0,0,true);
					//fracct->GetNumberofStructures() = 1;//only interested in lowest free energy structure
   					energy = fracct->GetEnergy(1);
					
				}
				else if(Usesub==1){
					strcpy(savefile, ct->GetCtLabel(1).c_str());
					savefile[strlen(savefile)-1]='\0'; 
					sprintf(pos,"_%d_s1_%d.sav",i,foldsize+length);
					strcat(savefile,pos);
					std::ifstream sav(savefile,std::ios::binary);
					if (TESTMODE && sav) {
						//read information from file if it was found
						read(&sav,&numofstructures);
						energyarray = new int [numofstructures+1];
						for (j=1;j<=numofstructures;j++)	read(&sav, energyarray+j);
						sav.close();
					}
					else {
						sav.close();
						alltrace(false,ctenergy,fracct,&data, percent, delta, NULL,NULL);
						numofstructures=fracct->GetNumberofStructures();
						energyarray = new int [numofstructures+1];
						for (j=1;j<=numofstructures;j++)		energyarray[j]=ctenergy[j];
					}

					
   				}
				else if(Usesub==3){
					dynamic(fracct,&data,1000,10,0);
					//save the energies in a save file
					strcpy(savefile, ct->GetCtLabel(1).c_str());
					savefile[strlen(savefile)-1]='\0'; 
					sprintf(pos,"_%d_s3_%d.sav",i,foldsize+length);
					strcat(savefile,pos);
					if (WRITE) {
						std::ofstream sav(savefile,std::ios::binary);
						int localint = fracct->GetNumberofStructures();
						write(&sav,&(localint));
						for (j=1;j<=fracct->GetNumberofStructures();j++) {
							int localint = fracct->GetEnergy(j);
							write(&sav, & (localint));

						}
						sav.close();
					}
					numofstructures=fracct->GetNumberofStructures();
					energyarray = new int [numofstructures+1];
					for (j=1;j<=numofstructures;j++)		energyarray[j]=fracct->GetEnergy(j);
   				}
				else if(Usesub==2) {
					strcpy(pfsfile, ct->GetCtLabel(1).c_str());
					pfsfile[strlen(pfsfile)-1]='\0'; 
					sprintf(pos,"_%d_s2_%d.pfs",i,foldsize+length);
					strcat(pfsfile,pos);
			
					//not using prefilter, so arrays can be reused when region move to the right
					//fold the first region, the folded region begin to move to right in the middle of target
					if (prefilter->useit == 0) {
						if (i==start) {
							target=new OligoPclass(fracct,pfdata);
							targetcopy=new OligoPclass(fracct,pfdata);
							target->partition(true,&Q,NULL);
							//char *report1="report1.out";				
							//target->partition(true,&Q,NULL,report1);
						}
						//when folded region begin moving,reuse some arrays overlapped expect for those on the edges
						else if((i-foldsize/2)>1 && (i+length-1+foldsize/2)<= (ct->GetSequenceLength()) ){
				   		//char *report2="report2.out";		
						target->scanfill(fracct,&Q,0);
						}
						//copy the arrays to be reused for next scan region
						//folded region is not moving at two ends, so copy the array outside for next folding without constrain
						if ( (i-foldsize/2)<1 || (i+length-1+foldsize/2)>=(ct->GetSequenceLength()) ) {
							scancopyend(target,targetcopy);
						}
						//folded region begin moving next, copy the overlapped region outside with different index
						else	scancopy(target,targetcopy);
				
					}
					//using prefilter, not reuse arrays,refold the new region every time
					else {
						
						//refold the new region without using any information from previous folding 
						target->reset4oligo(fracct);
						if (WRITE) {
							target->partition(true,&Q,NULL,pfsfile);
							
						}
						else target->partition(true,&Q,NULL);
					}
				}
				else if(Usesub==4) {
					strcpy(pfsfile, ct->GetCtLabel(1).c_str());
					pfsfile[strlen(pfsfile)-1]='\0'; 
					sprintf(pos,"_%d_s2_%d.pfs",i,foldsize+length);
					strcat(pfsfile,pos);
					//not using prefilter, so arrays can be reused when region move to the right
					//fold the first region, the folded region begin to move to right in the middle of target
					if (prefilter->useit == 0) {
						std::cout <<"Must use prefilter!\n";
						return;
					}
					//using prefilter, not reuse arrays,refold the new region every time
					else {
						std::ifstream pfs(pfsfile,std::ios::binary);
						if(TESTMODE && pfs) {
							pfs.close();
						}
						else{
							pfs.close();
							//refold the new region without using any information from previous folding 
							target->reset4oligo(fracct);
							target->partition(true,&Q,NULL,pfsfile);
							
						}
						stochastic(fracct,pfsfile,1000,5);
						efn2(&data,fracct);
					/*	ofstream sout("RL_sample_oligo_scan.out");
						for (j=1;j<=1000;j++)	{ sout << fracct->energy[j]<<"\n";	}
						sout.close();
					*/
						numofstructures=fracct->GetNumberofStructures();
						energyarray = new int [numofstructures+1];
						for (j=1;j<=numofstructures;j++)		energyarray[j]=fracct->GetEnergy(j);
					}
				}

				
				
				
				//set the single-stranded constrain 
				//fracct->nnopair=length;
				//fracct->checknopair();
				fracct->RemoveConstraints();
				if (i-foldsize/2<=1) {//constrained positions are different at two ends, as the folded region did not move
					for (j=0;j<length;j++) {
						fracct->AddSingle(i+j); 
					}
					foldstart=i;
					foldstop=i+length-1;
				} 
				else if(i+length-1+foldsize/2>= (ct->GetSequenceLength()) ) {
					for (j=0;j<length;j++) {
						fracct->AddSingle(foldsize+length-(ct->GetSequenceLength())+i+j);
						//fracct->nopair[j+1] = foldsize+length-(ct->GetSequenceLength())+i+j;
					}
					foldstart=foldsize+length-(ct->GetSequenceLength()) +i;
					foldstop=foldsize+length+length-1-(ct->GetSequenceLength()) +i;
				}
				else {//folded region begin to move, constrained position will be always in the middle of this region
					for (j=0;j<length;j++) {
						fracct->AddSingle(foldsize/2+1+j);
						//fracct->nopair[j+1] = foldsize/2+1+j;
					}
					foldstart=foldsize/2+1;
					foldstop=foldsize/2+length;
				}

				//---------------------------------------------------------------------
				//refold with constrain:
				//option 2 + refold + consider only the first suboptimal structure
				if (Usesub==0) {
					dynamic(fracct,&data,1000,10,0,0,true);
					table[i][2] = fracct->GetEnergy(1) - energy;			
				}
				//---------------------------------------------------------------------------------------------
				//option 2 + refold + consider all suboptimal structures ensemble energy
				else if (Usesub==1) {
					
					//sum of unconstrained energy from energyarray[]
					Q = (PFPRECISION) 0;
					for (k=1;k<=numofstructures;k++) {
	          			fn = exp(-(((long double)energyarray[k]-energyarray[1])/(RT_37C*conversionfactor)));
						Q = Q + fn;
         			}
					//calculate the whole set of suboptimal structures with constrain
					alltrace(false,ctenergy,fracct,&data, percent, delta, NULL,NULL);
					//sum of constrained energy now
					Qc = (PFPRECISION) 0;
					for (k=1;k<=fracct->GetNumberofStructures();k++) {
	          			fn = exp(-(((long double)ctenergy[k]-energyarray[1])/(RT_37C*conversionfactor)));
						Qc = Qc + fn;
         			}
					//final free energy difference					
					table[i][2] = (int)( conversionfactor*RT_37C*log( (long double)Q/(long double)Qc ) );

					numofsubstructures[i][0]=numofstructures; //report in report() function
					numofsubstructures[i][1]=fracct->GetNumberofStructures(); 
									      			
					delete[] energyarray;
				}
				//---------------------------------------------------------------------------------------------
				//option 2 + refold + consider heuristic suboptimal structures average free energy
				else if (Usesub==3) {
      				sum = 0;
					sn = 0;
					for (k=1;k<=numofstructures;k++) {
	          				fn = exp(-(((long double)energyarray[k]-energyarray[1])/(RT_37C*conversionfactor)));
							sn = sn + fn;
         					sum = sum + fn*((long double)(energyarray[k]));
					}
					table[i][2] = (int)(sum/sn);
					delete[] energyarray;

					//dynamic(fracct,&data,1000,10,0,0,true);
      				dynamic(fracct,&data,1000,10,0);
					//save the energies in a save file
					strcpy(savefile, ct->GetCtLabel(1).c_str());
					savefile[strlen(savefile)-1]='\0'; 
					sprintf(pos,"_%d_s3_%d_constrain.sav",i,foldsize+length);
					strcat(savefile,pos);
					if (WRITE) {
						std::ofstream sav(savefile,std::ios::binary);
						int localint = fracct->GetNumberofStructures();
						write(&sav,&(localint));
						for (j=1;j<=fracct->GetNumberofStructures();j++)		{
							int localint=fracct->GetEnergy(j);
							write(&sav, & (localint));
						}
						sav.close();
					}
					
   					sum = 0;
					sn = 0;
					for (k=1;k<=fracct->GetNumberofStructures();k++) {
	          				fn = exp(-(((long double)fracct->GetEnergy(k)-fracct->GetEnergy(1))/(RT_37C*conversionfactor)));
							sn = sn + fn;
         					sum = sum + fn*((long double)fracct->GetEnergy(k));
							//cout<< fracct->energy[k]<<"\n";
					}
					table[i][2] = (int)(sum/sn) - table[i][2];
					numofsubstructures[i][0]=numofstructures; //report in report() function
					numofsubstructures[i][1]=fracct->GetNumberofStructures(); 
				}
				//---------------------------------------------------------------------------------------------
				//option 2 + refold + partionfunction
				else if (Usesub==2) {
					strcpy(pfsfile, ct->GetCtLabel(1).c_str());
					pfsfile[strlen(pfsfile)-1]='\0'; 
					sprintf(pos,"_%d_s2_%d_constrain.pfs",i,foldsize+length);
					strcat(pfsfile,pos);
					//reuse the arrays filled without constrained 
					if (WRITE) {
						target->scanconstrain(fracct,&Qc,foldstart,foldstop,rescaleinrefill,pfsfile);

					}
					else target->scanconstrain(fracct,&Qc,foldstart,foldstop,rescaleinrefill);
					//exchange arrays to be used for next folding site without constrain
					if (prefilter->useit==0) {//not using targetcopy when prefilter is used
					targettemp=target;
					target=targetcopy;
					targetcopy=targettemp;
					}
					table[i][2] = (int)(conversionfactor*RT_37C*
					( (long double)(fracct->GetSequenceLength())*log(rescaleinrefill) 
						+ log((long double)Q/(long double)Qc) ) );
				}
				//---------------------------------------------------------------------------------------------
				//option 2 + refold + stochastic sample
				else if (Usesub==4) {			
					
   					sum = 0;
					for (k=1;k<=numofstructures;k++) {
	          				sum = sum + energyarray[k];
					}
					table[i][2] = (int)(sum/numofstructures);
					delete[] energyarray;

					strcpy(pfsfile, ct->GetCtLabel(1).c_str());
					pfsfile[strlen(pfsfile)-1]='\0'; 
					sprintf(pos,"_%d_s2_%d_constrain.pfs",i,foldsize+length);
					strcat(pfsfile,pos);
					std::ifstream pfs(pfsfile,std::ios::binary);
					if(TESTMODE && pfs) {
						pfs.close();
					}
					else{
						pfs.close();
						//reuse the arrays filled without constrained 
						target->scanconstrain(fracct,&Qc,foldstart,foldstop,rescaleinrefill,pfsfile);
						
					}
					stochastic(fracct,pfsfile,1000,5);
					efn2(&data,fracct);
					sum = 0;
					for (k=1;k<=fracct->GetNumberofStructures();k++) {
	          			sum = sum + fracct->GetEnergy(k);
							
					}
					table[i][2] = (int)(sum/fracct->GetNumberofStructures()) - table[i][2];
					numofsubstructures[i][0]=numofstructures; //report in report() function
					numofsubstructures[i][1]=fracct->GetNumberofStructures(); 

					//exchange arrays to be used for next folding site without constrain
					if (prefilter->useit==0) {//not using targetcopy when prefilter is used
					/*targettemp=target;
					target=targetcopy;
					targetcopy=targettemp;
					*/
						std::cout << "must use prefilter for stochasitic folding\n";
						return;
					}
				}
			
			}
		}//end option 2

		//-----------------------------------------------------------------------------------------------------
		//option 3: no local structure considered

		else {
			table[i][2] = 0;
		}

		
	
		//-----------------------------------------------------------------------------------------------------
		//-----------------------------------------------------------------------------------------------------
		//calculate the stability of the hybrid duplex
		
		//oligo is DNA:
		if (isdna) {
			
			table[i][1] = hybriddata->init;//initiation
			for (j=0;j<(length-1);j++) {
			
				table[i][1] += 
					hybriddata->stack[ct->numseq[i+j]][complement(i+j,ct)][ct->numseq[i+j+1]][complement(i+j+1,ct)];
			}
		}

		//oligo is RNA:
		else {
			table[i][1] = data.init;//initiation
      		for (j=0;j<(length-1);j++) {
          		table[i][1] += 
					data.stack[ct->numseq[i+j]][complement(i+j,ct)][ct->numseq[i+j+1]][complement(i+j+1,ct)];
			 }
			//consider AU end effects for RNA/RNA duplexes
			if (ct->numseq[i]==1||ct->numseq[i]==4) table[i][1]+=data.auend;
			if (ct->numseq[i+length-1]==1||ct->numseq[i+length-1]==4) table[i][1]+=data.auend;

		}

		//calculate the Tm of the duplex
		ds = helixstack->dsi;
		dh = helixstack->dhi;
		for (j=0;j<(length-1);j++) {
	  		dh +=helixstack->dh[ct->numseq[i+j]][complement(i+j,ct)][ct->numseq[i+j+1]][complement(i+j+1,ct)];
			ds +=helixstack->ds[ct->numseq[i+j]][complement(i+j,ct)][ct->numseq[i+j+1]][complement(i+j+1,ct)];
   		}
		if (ct->numseq[i]==1||ct->numseq[i]==4) {
			dh = dh + helixstack->dha;
			ds = ds + helixstack->dsa;
		}
		if (ct->numseq[i+length-1]==1||ct->numseq[i+length-1]==4) {
			dh = dh + helixstack->dha;
			ds = ds + helixstack->dsa;
		}
	      
		table[i][5] = (int) ((conversionfactor*((double)dh*1000)/
								((double)ds+conversionfactor*Rgas*log(c)))-273.15*conversionfactor);



		//-----------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------
		//calculate the free energy of intramolecular folding of oligo

	    
		if (Usesub==2) { //-2 is unavailabe, so not using partition function as interoligo cannot use it
			//reuse some arrays when scanning along the sequence for intramolecule
			//only calculate the partition function for the first one
			//cannot reuse array when prefilter is used
			if (prefilter==0) {
				if (i==start) {	
					intraoligo->reset4oligo(oligo1);
					intraoligo->partition(true,&Qolig);
				
				}
				//reuse some arrays overlapped, only change the index of left binding site
				//be careful that the oligo is complementary to the target, copy direction is reversed
				else {
				
					for (ip=length-1;ip>2;ip--) {
						for (jp=length-1;jp>=ip;jp--) {
								
							intraoligo->wca[ip][jp]=intraoligo->wca[ip-1][jp-1];
							intraoligo->w->f(ip,jp)=intraoligo->w->f(ip-1,jp-1);
							intraoligo->v->f(ip,jp)=intraoligo->v->f(ip-1,jp-1);
							intraoligo->wmb->f(ip,jp)=intraoligo->wmb->f(ip-1,jp-1);
							intraoligo->wl->f(ip,jp)=intraoligo->wl->f(ip-1,jp-1);
							intraoligo->wmbl->f(ip,jp)=intraoligo->wmbl->f(ip-1,jp-1);
							intraoligo->wcoax->f(ip,jp)=intraoligo->wcoax->f(ip-1,jp-1);
						}
					}
					//oligo is complementary to target, set reverse=1 for scanfill()
					intraoligo->scanfill(oligo1,&Qolig,1);

				}
			}
			else {
				intraoligo->reset4oligo(oligo1);
				intraoligo->partition(true,&Qolig);
			}
			oligo1->SetEnergy(1,(int)( conversionfactor*RT_37C*
				((long double)length*log( (intraoligo->data)->scaling) - log((long double)Qolig)) ));

		}
		else {//use the lowest free energy for Usesub 0 and 1 since oligo is small
    		if (isdna)	dynamic(oligo1,&ddata,1000,10,0,0,true);
			else	dynamic(oligo1,&data,1000,10,0,0,true);
		}
	    	
		table[i][3] = oligo1->GetEnergy(1);
		//if the intra structure is unfavorable , not consider it in the total energy
		if (table[i][3]>0) {
			table[i][3]=0;
		}

		//-----------------------------------------------------------------------------------------
		//calculate the free energy of intermolecular folding of oligo
		oligo->intermolecular = true;
		//oligo->GetSequenceLength() = 2*length + 3;
		//if (i==start) oligo->allocate(2*length + 3);
		for (j=1;j<=length;j++) {
			oligo->numseq[j+length+3] = oligo->numseq[j];
		}
		for (j=1;j<=3;j++) {
			oligo->numseq[length+j] = 5;
		}
	    
		if (Usesub==2) {	//-2 is unavailabe, as partition function is not availabe for inter oligo	
				
			interoligo->reset4oligo(oligo);
			interoligo->partition(true,&Qolig);
					
			oligo->SetEnergy(1,(int)( conversionfactor*RT_37C*
				( (long double)(oligo->GetSequenceLength())*log((interoligo->data)->scaling)  - log((long double)Qolig) ) )) ;
		}
		else {
			if (isdna)	dynamic(oligo,&ddata,1000,10,0,0,true);
			else	dynamic(oligo,&data,1000,10,0,0,true);
		}
			      	
		table[i][4] = oligo->GetEnergy(1);
		//if the intra structure is unfavorable , not consider it in the total energy
		if (table[i][4]>0) {
			table[i][4]=0;
		}
		


		//-----------------------------------------------------------------------------------------
		//-----------------------------------------------------------------------------------------
		//calculate the overall free energy for the binding
		k1b = exp(((-1)*(long double)(table[i][4]))/((1.9872)*(3.1)));
		k1u = exp(((-1)*(long double)(table[i][3]))/((1.9872)*(3.1)));

		//I don't know what's going on here.:)  ---John , Nov.9,2005
		//if (i==867) {
		//	table[0][0]=56;
		//}

		if ((k1u/(c*k1b))<100) {

      		dgeff = (int)(-((1.9872)*(3.1)*log(((4.0*k1b*(long double)(c))/(-1-k1u+sqrt(pow((long double)1.0+k1u,(long double) 2.0)
							+8.0*k1b*(long double)(c))))-1.0)));
		}
		else {
			dgeff = table[i][3];

		}
		//this line converts the free energy to the convention explained in
		//the written version of the algorithm:
		table[i][2] = -table[i][2];

		if ((table[i][2]<0)&&(dgeff<0)) {
	   		table[i][0] = table[i][1] + (int)((RT_37C*conversionfactor)*log( (exp(-(long double)(table[i][2])/
						  (RT_37C*conversionfactor))+1.0)*  (exp(-(long double)(dgeff)/(RT_37C*conversionfactor))+1.0) )+0.5);
		}
		else if (table[i][2]<0) {
			table[i][0] = table[i][1] + (int)((RT_37C*conversionfactor)*log( (exp(-(long double)(table[i][2])/
						  (RT_37C*conversionfactor))+1.0))+0.5);
		}
		else {
			table[i][0] = table[i][1] + (int)((RT_37C*conversionfactor)*log(  (exp(-(long double)(dgeff)/
						  (RT_37C*conversionfactor))+1.0) )+0.5);
		}

		//table[i][0] = table[i][1] + table[i][2] - dgeff;

	}//end of main loop over all oligonucleotides
	//The scan is finished now
	//-------------------------------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------  



	//clean up the ct file
	ct->SetEnergy(1,energy) ;

	//clean up memory use
	if (option==1) {
		if (Usesub==3 && foldsize==0) {
			ct->RemoveConstraints();//nnopair=0;
      		for (i=1;i<=ct->GetNumberofStructures();i++) delete[] temp2[i];
			delete[] temp2;
			delete[] energyarray;
		}
		else if (Usesub==0) delete[] temp;
		
	}
    else if (option ==2) {
   		ct->RemoveConstraints();//nnopair=0;
   		for (j=1;j<=ct->GetSequenceLength();j++) {
    		ct->SetPair(j,temp[j],1);//basepr[1][j] = temp[j];
		}
		delete[] temp;
		if (Usesub==1)	delete[] ctenergy;
		if ( (Usesub==1||Usesub==3||Usesub==4) && foldsize==0 )	delete[] energyarray;
		else if (Usesub==2 || Usesub ==4) {
			if(isdna)	delete dpfdata;
			delete pfdata;
			delete target;
			if (foldsize>0 && prefilter->useit == 0)	delete targetcopy;
			if (Usesub==2) {
				delete intraoligo ;
				delete interoligo ;
			}
		}
	}

	if (foldsize > 0) delete fracct;

	delete oligo;
	delete oligo1;
//cout <<"\n"<< numofstructures<<"\n";	
}


//=======================================================================
int readrd (rddata* data,char* dnarna) {
	
	int count,i,k,j,l;
	std::ifstream dr;
	char lineoftext[100];

	//make sure the file exists
	FILE *check;
	if ((check = fopen(dnarna, "r"))== NULL) {
	return 0;
	}

	fclose(check);
	dr.open(dnarna);
	/* Read info from stackdr */
	//add to the stack table the case where X (represented as 0) is looked up:
	for (count=1;count<=2;count++) dr >> lineoftext;//get past text in file
	dr >> lineoftext;
	data->init =(int)floor(conversionfactor*(atof(lineoftext)));

	for (count=1;count<=42;count++) dr >> lineoftext;//get past text in file
	for (i=0;i<=4;i++) {
		if (i!=0) for (count=1;count<=60;count++) dr >> lineoftext;
		for (k=0;k<=4;k++) {
			for (j=0;j<=4;j++) {
				for (l=0;l<=4;l++) {
					if ((i==0)||(j==0)||(k==0)||(l==0)) {
						data->stack[i][j][k][l]=0;
					}
					else {
						dr >> lineoftext;
						if (strcmp(lineoftext,".")){
							data->stack[i][j][k][l] =(int)floor(conversionfactor*(atof(lineoftext))+.5);
						}
						else data->stack[i][j][k][l] = INFINITE_ENERGY;
					}
					//cout <<"stack "<<i<<" "<<j<<" "<<k<<" "<<l<<"  "<<data->stack[i][j][k][l]<<"\n";
				}
				//cin >> m;
			}
		}
	}

	return 1;
}



//=======================================================================
//return the numerical equivalent of the complementary base to nucleotide i
int complement(int i, structure *ct) {
	
	int a;

	if (ct->numseq[i] == 0) return 0;
	else {
	   	a = 5 - ct->numseq[i];
		return a;
    }
}



//=======================================================================
//Function gets the names of data files to open for DNA-DNA parameters
void dgetdat(char *loop, char *stackf, char *tstackh, char *tstacki,
		     char *tloop, char *miscloop, char *danglef, char *int22,
			 char *int21,char *coax, char *tstackcoax,
		     char *coaxstack, char *tstack, char *tstackm,
			 char *triloop, char *int11, char *Path) {
	strcpy (loop,Path);
	strcpy (stackf,Path);
	strcpy (tstackh,Path);
	strcpy (tstacki,Path);
	strcpy (tloop,Path);
	strcpy (miscloop,Path);
	strcpy (danglef,Path);
	strcpy (int22,Path);
	strcpy (int21,Path);
	strcpy (triloop,Path);
	strcpy (coax,Path);
	strcpy (tstackcoax,Path);
	strcpy (coaxstack,Path);
	strcpy (tstack,Path);
	strcpy (tstackm,Path);
	strcpy (int11,Path);

	strcat (loop,"dnaloop.dat");
	strcat (stackf,"dnastack.dat");
	strcat (tstackh,"dnatstackh.dat");
	strcat (tstacki,"dnatstacki.dat");
	strcat (tloop,"dnatloop.dat");
	strcat (miscloop,"dnamiscloop.dat");
	strcat (danglef,"dnadangle.dat");
	strcat (int22,"dnaint22.dat");
	strcat (int21,"dnaint21.dat");
	strcat (triloop,"dnatriloop.dat");
	strcat (coax,"dnacoaxial.dat");
	strcat (tstackcoax,"dnatstackcoax.dat");
	strcat (coaxstack,"dnacoaxstack.dat");
	strcat (tstack,"dnatstack.dat");
	strcat (tstackm,"dnatstackm.dat");
	strcat (int11,"dnaint11.dat");
}



//=======================================================================
//This is a older post-filter of siRNA 
//use siRNA selection criteria to filter the output
//mask will contain true for those that meet the criteria
void filterbysirna ( structure *ct, int **table, int length, datatable *data, 
				     bool *mask, double asuf, double tofe, double fnnfe) {
	
	int iasuf,itofe,ifnnfe;
	int i,j,*k;
	k = new int [length];

	iasuf = (int) (conversionfactor*asuf);
	itofe = (int) (conversionfactor*tofe);
	ifnnfe = (int) (conversionfactor*fnnfe);


	for (i=1;i<=(ct->GetSequenceLength()-length+1); i++) {
		//filter all oligos
		mask[i]=true;
		if (iasuf>table[i][3]) {
			mask[i]=false;
		}
		if (itofe>table[i][2]) {
			mask[i]=false;
		}

		//filter out those that have a repeat of A, G, or U of more than 4. 
		if (length>4) {
			for (j=length-1;j>=0;j--) {
   				k[j] = complement(i+j,ct);
      		
			}
			for (j=0;j<length-3;j++) {
				if (k[j]==1) {
					if (k[j+1]==1&&k[j+2]==1&&k[j+3]==1)	mask[i]=false;
				}
				else if (k[j]==3) {
					if (k[j+1]==3&&k[j+2]==3&&k[j+3]==3)	mask[i]=false;
				}
				else if (k[j]==4) {
					if (k[j+1]==4&&k[j+2]==4&&k[j+3]==4)	mask[i]=false;
				}
			}
		}
		table[i][6]= data->stack[complement(i+length-1,ct)][ct->numseq[i+length-1]]
								[complement(i+length-2,ct)][ct->numseq[i+length-2]];
		//account for change in AU end if necessary:
		if (ct->numseq[i+length-1]==1||ct->numseq[i+length-1]==4) {
			if (ct->numseq[i+length-1]==2||ct->numseq[i+length-1]==3) {
				table[i][6]+=data->auend;
			}
		}
		else {
			if (ct->numseq[i+length-1]==1||ct->numseq[i+length-1]==4) {
				table[i][6]-=data->auend;
			}
		}
		if (ifnnfe>table[i][6])
			mask[i]=false;

	}

	delete[] k;
}


//=======================================================================
//output a report
void report(const char* filename, structure *ct, int **table,int **numofsubstructures, const int length, const bool isdna,
			const double conc, const int Usesub,const int start,const int stop,const siPREFILTER *prefilter,const int foldsize, 
			const bool *mask, const double asuf, const double tofe, const double fnnfe,const bool isHTML) {
	
	int i,j,k;
	bool use;

	//If this is not supposed to be html, use the alternative report below
	if (!isHTML) {
		report(filename,ct,table,length,isdna,conc,Usesub,mask,asuf,tofe,fnnfe);
		return;
	}

	//#if defined (debugmode)
	std::ofstream out(filename);
	//#endif
	
	if (mask!=NULL) {
		out << "Oligonucleotides filtered by siRNA selection criteria: ";
		out << "\tAntisense strand unimolecular folding free energy >= "<<asuf<<"<br>\n";
		out << "\tTarget opening free energy >= "<<tofe<<"<br>\n";
		out << "\tFirst nearest neighbor free energy >= "<<fnnfe<<"<br>\n";
		out << "\tSequences with more than three G's, U's, or A's in a row removed.<br>\n";
		out << "Antisense strand shown 5' to 3'.<br>\n";
	}
	out <<"<br><br>\n";
 	out <<"<h3>Energy table:</h3><br>\n";
	out <<"<table>\n";
	out <<"<tr> <td>Pos.</td>\t<td>Oligo(5\'->3\')</td><td> </td>\t\t\t\t\t<td>Overall</td>\t";
	out <<"<td>Duplex</td>\t<td>Tm-Dup</td>\t<td>Break-targ.</td>";
	out <<"\t<td>Intraoligo</td>\t<td>Interoligo</td>\t<td>End_diff</td>\t<td>prefilter_score</td>";
	if (Usesub==1 || Usesub==3||Usesub==4)	out <<"\t<td>structures #</td>\t<td>constrained structures #</td>";
	if (mask!=NULL) std::cout << "\t<td>First NN Opening</td>";
	out <<"</tr><tr> </tr>\n\n";
	out <<"\t\t\t\t\t\t<tr> <td> </td><td> </td> <td><td>kcal/mol</td>\t<td>kcal/mol</td>\t";
	out <<"<td>degC</td>\t\t<td>kcal/mol</td>\t\t<td>kcal/mol</td>\t<td>kcal/mol</td>\t<td>kcal/mol</td>";
	if (mask!=NULL) "\t <td>kcal/mol</td>"; 
	out << "</tr>\n";

	

	for (i=start;i<=stop; i++) {

		if (prefilter->useit != 0) {
			if (prefilter->score[i] < FILTER_PASS) {
				continue;
			}
		}

 
		//output each possible oligo
		if (mask!=NULL) use = mask[i];
		else use = true;
		
		if (use) {
			out << "<tr><td>"<<i <<"</td>\t";
			out << "<td>";
			for (j=length-1;j>=0;j--) {
   				k = complement(i+j,ct);
      			if (isdna) {
         			if (k==0) out << "X";
					else if (k==1) out << "A";
					else if (k==2) out << "C";
					else if (k==3) out << "G";
					else if (k==4) out << "T";
				}
				else {
		     		if (k==0) out << "X";
					else if (k==1) out << "A";
					else if (k==2) out << "C";
					else if (k==3) out << "G";
					else if (k==4) out << "U";
				}
			}
			out << "</td><td> </td>\t\t";
		
			if (prefilter->useit != 0) {
		
				if (prefilter->score[i] < FILTER_PASS) {
					out<<"<td>-</td>\t\t<td>-</td>\t\t<td>-</td>\t\t<td>-</td>\t\t\t<td>-</td>\t\t<td>-</td>\t\t</tr>\n";
					continue;
				}
			}

			out <<"<td>"<< ((double) (table[i][0]))/conversionfactor << "</td>\t\t";
			out <<"<td>"<< ((double) (table[i][1]))/conversionfactor << "</td>\t\t";
			out <<"<td>"<< ((double) (table[i][5]))/conversionfactor << "</td>\t\t";
			out <<"<td>"<< ((double) (table[i][2]))/conversionfactor << "</td>\t\t\t";
			out <<"<td>"<< ((double) (table[i][3]))/conversionfactor << "</td>\t\t";
			out <<"<td>"<< ((double) (table[i][4]))/conversionfactor << "</td>\t\t";
			
			



			if (prefilter->useit != 0)	{
				out <<"<td>"<<  prefilter->enddiff[i] <<"</td>\t\t";
				out <<"<td>"<<  prefilter->score[i] <<"</td>\t\t";
			}
			if (Usesub==1 || Usesub==3||Usesub==4)	 {
				out <<"<td>"<< numofsubstructures[i][0]  <<"</td>\t\t";
				out <<"<td>"<< numofsubstructures[i][1]  <<"</td>\t\t";
			}
			
			if (mask!=NULL) {
		
				out <<"<td>"<< ((double) (table[i][6])/conversionfactor) << "</td>\t\t";
		
			}
			out << "</tr>\n";

			
		}
	}
	out<< "</table>\n";
//	cout<< "</body></html>\n";

	



}

//output tab delimited report
void report(const char* filename, structure *ct, int **table, const int length, const bool isdna,
	const double conc, const int Usesub,const bool *mask, const double asuf, const double tofe, const double fnnfe) {
std::ofstream out(filename);
int i,j,k;
bool use;


out << "Oligo walk for:  "<<ct->GetSequenceLabel()<<"\n\n";
if (mask!=NULL) {
	out << "Oligonucleotides filtered by siRNA selection criteria: ";
	out << "\tAntisense strand unimolecular folding free energy >= "<<asuf<<"\n";
	out << "\tTarget opening free energy >= "<<tofe<<"\n";
	out << "\tFirst nearest neighbor free energy >= "<<fnnfe<<"\n";
	out << "\tSequences with more than three G's, U's, or A's in a row removed.\n";
	out << "Antisense strand shown 5' to 3'.\n";
}
else {
	if (isdna) out << "DNA oligo";
	else out << "RNA oligo";
}
out << "       length: "<<length<<"   concentration: "<<conc<<"\n\n";
if (Usesub) out << "Suboptimal structures were used.  There are "<<ct->GetNumberofStructures()<<" suboptimal structures.";
else out << "Suboptimal structures were not used.\n\n";
out <<"position\toligo\ttotal binding\tduplex formation\tTm of Duplex\ttarget structure";
out <<"\tintramolecular oligo\tintermolecular oligo";
if (mask!=NULL) out << "\tFirst NN Opening";
out <<"\n";
out <<"\t\tkcal/mol\tkcal/mol\tdeg C\tkcal/mol\tkcal/mol\tkcal/mol";
if (mask!=NULL) "\tkcal/mol"; 
out << "\n\n";


for (i=1;i<=(ct->GetSequenceLength()-length+1); i++) {
	//output each possible oligo
	if (mask!=NULL) use = mask[i];
	else use = true;
	if (use) {
		out << i <<"\t";
		for (j=length-1;j>=0;j--) {
   			k = complement(i+j,ct);
      		if (isdna) {
         		if (k==0) out << "X";
				else if (k==1) out << "A";
				else if (k==2) out << "C";
				else if (k==3) out << "G";
				else if (k==4) out << "T";
			 }
			else {
         		if (k==0) out << "X";
				else if (k==1) out << "A";
				else if (k==2) out << "C";
				else if (k==3) out << "G";
				else if (k==4) out << "U";
			}
		}
		out << "\t";
		out << ((double) (table[i][0]))/conversionfactor << "\t";
		out << ((double) (table[i][1]))/conversionfactor << "\t";
		out << ((double) (table[i][5]))/conversionfactor << "\t";
		out << ((double) (table[i][2]))/conversionfactor << "\t";
		out << ((double) (table[i][3]))/conversionfactor << "\t";
		out << ((double) (table[i][4]))/conversionfactor << "\t";
		if (mask!=NULL) {
		
			out << ((double) (table[i][6])/conversionfactor) << "\t";
		
		}
		out << "\n";
	}
}


}


//=================================================================================================
//copy the arrays to be reused with different index when the folded region move one nucleotide to the right
inline void scancopy(OligoPclass *region, OligoPclass *copyregion) {

	int j,i;
	//only copy the overlapped region which is in the middle of sequence
	for (i=2;i<=(copyregion->number)-2;i++) {
		for (j=i;j<=(copyregion->number)-2;j++) {
			//w5 is not copied as it will be recalc when i==1						
			copyregion->wca[i][j]=region->wca[i+1][j+1];
			copyregion->w->f(i,j)=region->w->f(i+1,j+1);
			copyregion->v->f(i,j)=region->v->f(i+1,j+1);
			copyregion->wmb->f(i,j)=region->wmb->f(i+1,j+1);
			copyregion->wl->f(i,j)=region->wl->f(i+1,j+1);
			copyregion->wmbl->f(i,j)=region->wmbl->f(i+1,j+1);
			copyregion->wcoax->f(i,j)=region->wcoax->f(i+1,j+1);
			
			

		}
	}
}
//copy every arrays to be reused when the folded region did not move right yet 
inline void scancopyend(OligoPclass *region, OligoPclass *copyregion) {

	int j,i;
	for (i=1;i<=copyregion->number;i++) {
		for (j=i;j<=copyregion->number;j++) {
						
			if(i==1)	copyregion->w5[j]=region->w5[j];
			copyregion->wca[i][j]=region->wca[i][j];
			copyregion->w->f(i,j)=region->w->f(i,j);
			copyregion->v->f(i,j)=region->v->f(i,j);
			copyregion->wmb->f(i,j)=region->wmb->f(i,j);
			copyregion->wl->f(i,j)=region->wl->f(i,j);
			copyregion->wmbl->f(i,j)=region->wmbl->f(i,j);
			copyregion->wcoax->f(i,j)=region->wcoax->f(i,j);
		
		}
	}
}

