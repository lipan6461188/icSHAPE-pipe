

#include "stochastic.h"
#include "random.h"
#include "pfunction.h"
#include "stackclass.h"

#include <iostream>
#include <cstdlib>
#include <cstring>
using namespace std;


//Define 1 for the purposes of taking an inverse
static PFPRECISION ONE=1; //include code for extended double if needed


 
//register a base pair between two nucleotides
inline void regbp(structure *ct, int structurenumber, short i, short j) {
	ct->SetPair(i,j,structurenumber);

}



int stochastictraceback(DynProgArray<PFPRECISION> *w,DynProgArray<PFPRECISION> *wmb,DynProgArray<PFPRECISION> *wmbl,DynProgArray<PFPRECISION> *wcoax,DynProgArray<PFPRECISION> *wl,DynProgArray<PFPRECISION> *v,
	forceclass *fce, PFPRECISION *w3,PFPRECISION *w5,PFPRECISION scaling, bool *lfce, bool *mod, pfdatatable *data, int numberofstructures, 
	structure *ct, int randomseed, ProgressHandler *progress) {
	double scalinginv,twoscaling;
	int number;
	int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
	{0,1,0,1,0,0},{0,0,0,0,0,0}};
	scalinginv = ONE/data->scaling;
	twoscaling=data->scaling*data->scaling;
	int tracebackerror=0;

	//allocation of space for the sampled structures
	//must be done OUTSIDE the main loop (which might
	//be run in parallel)
	for(int i = 1; i <= numberofstructures; i++){
		ct->AddStructure();
		ct->SetCtLabel(ct->GetSequenceLabel(), i);
	}

	#ifdef SMP
	#pragma omp parallel for
	#endif
	for (number = 1; number <= numberofstructures; number++) {
		randomnumber rand;
		rand.seed(randomseed+number);
		double roll;
		double cumulative, denominator;
		stackclass stack;
		integersize dummy1;
		short dummy2;
		bool found;
		short switchcase,i,j,k,ip,jp,d;
		if (progress!=NULL) {

			progress->update((int) (((float) 100*number)/((float) numberofstructures)));

		}
		//start by putting the whole fragment on the stack:
		stack.push(1,ct->GetSequenceLength(),0,0,0);

		while (stack.pull(&i,&j,&switchcase,&dummy1,&dummy2)) {
			roll= rand.roll();
			cumulative = 0;
			found = false;
			switch(switchcase) {
				case 0: //switchcase=0, dealing with w5 fragment

					//Try adding a nucleotide to existing w5 fragment:
					if (j==0) found = true;

					if (!lfce[j]&&!found) {
						cumulative = ((w5[j-1]*scaling)/w5[j]);
						if (cumulative > roll) {
							stack.push(1,j-1,0,0,0);
							found=true;
						}
					}
					
					for (k=0;k<=j-4&&!found;k++) {

						#ifdef SIMPLEMBLOOP
						cumulative+=(w5[k]*erg4(j-1,k+1,j,1,ct,data,lfce[j])*v->f(k+1,j-1)*penalty(j-1,k+1,ct,data))/w5[j];
						if (cumulative>roll&&!found) {
							stack.push(1,k,0,0,0);
							stack.push(k+1,j-1,1,0,0);
							found = true;
							

						}
						if (!found&&(mod[k+1]||mod[j-1])) if(inc[ct->numseq[k+2]][ct->numseq[j-2]]&&notgu(k+1,j-1,ct)&&!(fce->f(k+1,j-1)&SINGLE)) {

							cumulative+=(w5[k]*erg4(j-1,k+1,j,1,ct,data,lfce[j])*v->f(k+2,j-2)
								*penalty(j-1,k+1,ct,data)*erg1(k+1,j-1,k+2,j-2,ct,data))/w5[j];
							if (cumulative>roll&&!found) {
								stack.push(1,k,0,0,0);
								stack.push(k+2,j-2,1,0,0);
								regbp(ct,number,k+1,j-1);
								found=true;
								
							}
						}
						#else  //!SIMPLEMBLOOP
						cumulative +=(w5[k]*v->f(k+1,j)*penalty(j,k+1,ct,data))/w5[j];
						if (cumulative > roll) {
							stack.push(1,k,0,0,0);
							stack.push(k+1,j,1,0,0);
							found = true;
						}
						if (!found&&(mod[k+1]||mod[j])) if(inc[ct->numseq[k+2]][ct->numseq[j-1]]&&notgu(k+1,j,ct)&&!(fce->f(k+1,j)&SINGLE)) {
							cumulative +=(w5[k]*v->f(k+2,j-1)*penalty(j,k+1,ct,data)
								*erg1(k+1,j,k+2,j-1,ct,data))/w5[j];
							if (cumulative > roll) {
								stack.push(1,k,0,0,0);
								stack.push(k+2,j-1,1,0,0);
								regbp(ct,number,k+1,j);
								found=true;
							}
						}

						cumulative+=(w5[k]*erg4(j,k+2,k+1,2,ct,data,lfce[k+1])*v->f(k+2,j)*penalty(j,k+2,ct,data))/w5[j];
						if (cumulative > roll&&!found) {
							stack.push(1,k,0,0,0);
							stack.push(k+2,j,1,0,0);
							found=true;
			
						}
						if (!found&&(mod[k+2]||mod[j])) if(inc[ct->numseq[k+3]][ct->numseq[j-1]]&&notgu(k+1,j-1,ct)&&!(fce->f(k+1,j-1)&SINGLE)) {
							cumulative += (w5[k]*erg4(j,k+2,k+1,2,ct,data,lfce[j])*v->f(k+3,j-1)
								*penalty(j,k+2,ct,data)*erg1(k+2,j,k+3,j-1,ct,data))/w5[j];
							if (cumulative>roll) {
								stack.push(1,k,0,0,0);
								stack.push(k+2,j,1,0,0);
								regbp(ct,number,k+3,j-1);
								found=true;

							}
						}

						cumulative+=(w5[k]*erg4(j-1,k+1,j,1,ct,data,lfce[j])*v->f(k+1,j-1)*penalty(j-1,k+1,ct,data))/w5[j];
						if (cumulative>roll&&!found) {
							stack.push(1,k,0,0,0);
							stack.push(k+1,j-1,1,0,0);
							found = true;
							

						}
						if (!found&&(mod[k+1]||mod[j-1])) if(inc[ct->numseq[k+2]][ct->numseq[j-2]]&&notgu(k+1,j-1,ct)&&!(fce->f(k+1,j-1)&SINGLE)) {

							cumulative+=(w5[k]*erg4(j-1,k+1,j,1,ct,data,lfce[j])*v->f(k+2,j-2)
								*penalty(j-1,k+1,ct,data)*erg1(k+1,j-1,k+2,j-2,ct,data))/w5[j];
							if (cumulative>roll&&!found) {
								stack.push(1,k,0,0,0);
								stack.push(k+2,j-2,1,0,0);
								regbp(ct,number,k+1,j-1);
								found=true;
								
							}
						}
						cumulative+=(w5[k]*data->tstack[ct->numseq[j-1]][ct->numseq[k+2]][ct->numseq[j]][ct->numseq[k+1]] 
									*pfchecknp(lfce[j],lfce[k+1]) * v->f(k+2,j-1)*
									penalty(j-1,k+2,ct,data))/w5[j];
						if (cumulative>roll&&!found) {
							stack.push(1,k,0,0,0);
							stack.push(k+2,j-1,1,0,0);
							found=true;
							


						}

						if (!found&&(mod[k+2]||mod[j-1])) if(inc[ct->numseq[k+3]][ct->numseq[j-2]]&&notgu(k+2,j-1,ct)&&!(fce->f(k+2,j-1)&SINGLE)) {

							cumulative+=(w5[k]*data->tstack[ct->numseq[j-1]][ct->numseq[k+2]][ct->numseq[j]][ct->numseq[k+1]] 
									*pfchecknp(lfce[j],lfce[k+1]) * v->f(k+3,j-2)*
									penalty(j-1,k+2,ct,data)*erg1(k+2,j-1,k+3,j-2,ct,data))/w5[j];

							if (cumulative>roll&&!found) {
								stack.push(1,k,0,0,0);
								stack.push(k+3,j-2,1,0,0);
								regbp(ct,number,k+2,j-1);
								found=true;
								
							}

						}

						//recheck all the coaxial stacking possibilities:
						
						i = k+1;
						for (ip=i+minloop+1;ip<j-minloop-1&&!found;ip++) {
						
							//first consider flush stacking
							cumulative+=w5[k]*v->f(i,ip)*v->f(ip+1,j)*penalty(i,ip,ct,data)
								*penalty(ip+1,j,ct,data)*ergcoaxflushbases(i,ip,ip+1,j,ct,data)/w5[j];

							if (cumulative > roll) {
								stack.push(i,ip,1,0,0);
								stack.push(ip+1,j,1,0,0);
								found=true;
								stack.push(1,k,0,0,0);
							}


							if ((mod[i]||mod[ip]||mod[ip+1]||mod[j])) {

								if ((mod[i]||mod[ip])&&(mod[ip+1]||mod[j])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]
									&&inc[ct->numseq[ip+2]][ct->numseq[j-1]]&&notgu(i,ip,ct)&&notgu(ip+1,j,ct)
									&&!(fce->f(ip+1,j)&SINGLE)&&!(fce->f(i,ip)&SINGLE)) {

									cumulative+=w5[k]*v->f(i+1,ip-1)*v->f(ip+2,j-1)*penalty(i,ip,ct,data)
										*penalty(ip+1,j,ct,data)*ergcoaxflushbases(i,ip,ip+1,j,ct,data)
										*erg1(i,ip,i+1,ip-1,ct,data)*erg1(ip+1,j,ip+2,j-1,ct,data)/w5[j];

									if (!found&&cumulative>roll) {
										stack.push(i+1,ip-1,1,0,0);
										stack.push(ip+2,j-1,1,0,0);
										found=true;
										regbp(ct,number,i,ip);
										regbp(ct,number,ip+1,j);
										stack.push(1,k,0,0,0);
									}


								}

								if ((mod[i]||mod[ip])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]&&notgu(i,ip,ct)&&!(fce->f(i,ip)&SINGLE)) {
							
									cumulative+=w5[k]*v->f(i+1,ip-1)*v->f(ip+1,j)*penalty(i,ip,ct,data)
										*penalty(ip+1,j,ct,data)*ergcoaxflushbases(i,ip,ip+1,j,ct,data)
										*erg1(i,ip,i+1,ip-1,ct,data)/w5[j];

									if (!found&&cumulative>roll) {
										stack.push(i+1,ip-1,1,0,0);
										stack.push(ip+1,j,1,0,0);
										found=true;
										regbp(ct,number,i,ip);
										stack.push(1,k,0,0,0);

									}


								}

								if ((mod[ip+1]||mod[j])&&inc[ct->numseq[ip+2]][ct->numseq[j-1]]&&notgu(ip+1,j,ct)&&!(fce->f(ip+1,j)&SINGLE)) {


									cumulative+=w5[k]*v->f(i,ip)*v->f(ip+2,j-1)*penalty(i,ip,ct,data)
										*penalty(ip+1,j,ct,data)*ergcoaxflushbases(i,ip,ip+1,j,ct,data)
										*erg1(ip+1,j,ip+2,j-1,ct,data)/w5[j];

									if (!found&&cumulative>roll) {
										stack.push(i,ip,1,0,0);
										stack.push(ip+2,j-1,1,0,0);
										found=true;
										regbp(ct,number,ip+1,j);
										stack.push(1,k,0,0,0);
			
									}

								}


							}

							if (!lfce[ip+1]&&!lfce[j]) {
								//now consider an intervening mismatch
								cumulative+=w5[k]*v->f(i,ip)*v->f(ip+2,j-1)*penalty(i,ip,ct,data)
									*penalty(ip+2,j-1,ct,data)*ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)/w5[j];

								if (cumulative>roll&&!found) {
									stack.push(i,ip,1,0,0);
									stack.push(ip+2,j-1,1,0,0);
									found=true;
									stack.push(1,k,0,0,0);

								}

								if (mod[i]||mod[ip]||mod[ip+2]||mod[j-1]) {
									if ((mod[i]||mod[ip])&&(mod[ip+2]||mod[j-1])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]
										&&inc[ct->numseq[ip+3]][ct->numseq[j-2]]&&notgu(i,ip,ct)&&notgu(ip+2,j-1,ct)
											&&!(fce->f(i,ip)&SINGLE)&&!(fce->f(ip+2,j-1)&SINGLE)) {

										 cumulative+=w5[k]*v->f(i+1,ip-1)*v->f(ip+3,j-2)*penalty(i,ip,ct,data)
											*penalty(ip+2,j-1,ct,data)*ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)
											*erg1(i,ip,i+1,ip-1,ct,data)*erg1(ip+2,j-1,ip+3,j-2,ct,data)/w5[j];

										 if (!found&&cumulative>roll) {
											stack.push(i+1,ip-1,1,0,0);
											stack.push(ip+3,j-2,1,0,0);
											regbp(ct,number,i,ip);
											regbp(ct,number,ip+2,j-1);
											found=true;
											stack.push(1,k,0,0,0);

										 }


									}

									if ((mod[i]||mod[ip])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]&&notgu(i,ip,ct)&&!(fce->f(i,ip)&SINGLE)) {
							
										cumulative+=w5[k]*v->f(i+1,ip-1)*v->f(ip+2,j-1)*penalty(i,ip,ct,data)
											*penalty(ip+2,j-1,ct,data)*ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)
											*erg1(i,ip,i+1,ip-1,ct,data)/w5[j];

										if(!found&&cumulative>roll) {
											stack.push(i+1,ip-1,1,0,0);
											stack.push(ip+2,j-1,1,0,0);
											regbp(ct,number,i,ip);
											found=true;
											stack.push(1,k,0,0,0);

										}


									}

									if ((mod[ip+2]||mod[j-1])&&inc[ct->numseq[ip+3]][ct->numseq[j-2]]&&notgu(ip+2,j-1,ct)&&!(fce->f(ip+2,j-1)&SINGLE)) {


										cumulative+=w5[k]*v->f(i,ip)*v->f(ip+3,j-2)*penalty(i,ip,ct,data)
											*penalty(ip+2,j-1,ct,data)*ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)
											*erg1(ip+2,j-1,ip+3,j-2,ct,data)/w5[j];

										if (!found&&cumulative>roll) {
											stack.push(i,ip,1,0,0);
											stack.push(ip+3,j-2,1,0,0);
											regbp(ct,number,ip+2,j-1);
											found=true;
											stack.push(1,k,0,0,0);

										}

									}
								}
							}

							if(!lfce[i]&&!lfce[ip+1]) {
								cumulative+=w5[k]*v->f(i+1,ip)*v->f(ip+2,j)*penalty(i+1,ip,ct,data)
									*penalty(ip+2,j,ct,data)*ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)/w5[j];


								if (!found&&cumulative>roll) {
									stack.push(i+1,ip,1,0,0);
									stack.push(ip+2,j,1,0,0);
									stack.push(1,k,0,0,0);
									found=true;

								}

								if (mod[i+1]||mod[ip]||mod[ip+2]||mod[j]) {
									if ((mod[i+1]||mod[ip])&&(mod[ip+2]||mod[j])&&inc[ct->numseq[i+2]][ct->numseq[ip-1]]
										&&inc[ct->numseq[ip+3]][ct->numseq[j-1]]&&notgu(i+1,ip,ct)&&notgu(ip+2,j,ct)
										&&!(fce->f(i+1,ip)&SINGLE)&&!(fce->f(ip+2,j)&SINGLE)	) {

										cumulative+=w5[k]*v->f(i+2,ip-1)*v->f(ip+3,j-1)*penalty(i+1,ip,ct,data)
											*penalty(ip+2,j,ct,data)*ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)
											*erg1(i+1,ip,i+2,ip-1,ct,data)*erg1(ip+2,j,ip+3,j-1,ct,data)/w5[j];	

										if (!found&&cumulative>roll) {
											stack.push(i+2,ip-1,1,0,0);
											stack.push(ip+3,j-1,1,0,0);
											regbp(ct,number,i+1,ip);
											regbp(ct,number,ip+2,j);
											found=true;
											stack.push(1,k,0,0,0);

										}

							
									}
									if ((mod[i+1]||mod[ip])&&inc[ct->numseq[i+2]][ct->numseq[ip-1]]&&notgu(i+1,ip,ct)&&!(fce->f(i+1,ip)&SINGLE)) {
							
										cumulative+=w5[k]*v->f(i+2,ip-1)*v->f(ip+2,j)*penalty(i+1,ip,ct,data)
											*penalty(ip+2,j,ct,data)*ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)
											*erg1(i+1,ip,i+2,ip-1,ct,data)/w5[j];

										if (!found&&cumulative>roll) {
											stack.push(i+2,ip-1,1,0,0);
											stack.push(ip+2,j,1,0,0);
											regbp(ct,number,i+1,ip);
											stack.push(1,k,0,0,0);
										
											found=true;

										}


									}

									if ((mod[ip+2]||mod[j])&&inc[ct->numseq[ip+3]][ct->numseq[j-1]]&&notgu(ip+2,j,ct)&&!(fce->f(ip+2,j)&SINGLE)) {


										cumulative+=w5[k]*v->f(i+1,ip)*v->f(ip+3,j-1)*penalty(i+1,ip,ct,data)
											*penalty(ip+2,j,ct,data)*ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)
											*erg1(ip+2,j,ip+3,j-1,ct,data)/w5[j];

										if (!found&&cumulative>roll) {
											stack.push(i+1,ip,1,0,0);
											stack.push(ip+3,j-1,1,0,0);
											stack.push(1,k,0,0,0);
										
											regbp(ct,number,ip+2,j);
											found=true;

										}

									}
								}
						

							}
						}

						#endif //SIMPLEMBLOOP


					}

					if (!found) {
						cout << "Traceback error at w5\n";
						tracebackerror=14;

					}
         
      		
					break;
				case 1: //switchcase=1, dealing with a v fragment
					regbp(ct,number,i,j);

					//check to see if constant is used.  
					//If it is, then v->f(i,j) was multiplied by ct->constant[j][i] and
					//this influence needs to be removed.
					if (ct->constant!=NULL) denominator = v->f(i,j) / ct->constant[j][i];
					else denominator = v->f(i,j);

					//try closing a hairpin
					cumulative += (erg3(i,j,ct,data,fce->f(i,j)))/denominator;
					if (cumulative>roll) {
						found = true; //nothing to put on the stack

					}

					//try stacking on a previous pair
					if (!mod[i]&&!mod[j]) {
						cumulative+=erg1(i,j,i+1,j-1,ct,data)*v->f(i+1,j-1)/denominator;
						if (!found&&cumulative>roll) {
							stack.push(i+1,j-1,1,0,0);
							found=true;

						}

					}
					else {
						if ((ct->numseq[i]==3&&ct->numseq[j]==4)||(ct->numseq[i]==4&&ct->numseq[j]==3)) {
							cumulative+=erg1(i,j,i+1,j-1,ct,data)*v->f(i+1,j-1)/denominator;
							if (cumulative>roll) {
								stack.push(i+1,j-1,1,0,0);
								found=true;

							}

						}
						else if ((ct->numseq[i+1]==3&&ct->numseq[j-1]==4)||(ct->numseq[i+1]==4&&ct->numseq[j-1]==3)) {

							cumulative+=erg1(i,j,i+1,j-1,ct,data)*v->f(i+1,j-1)/denominator;
							if (cumulative>roll) {
								stack.push(i+1,j-1,1,0,0);
								found=true;

							}

						}
						else if (i-1>0) {
							if ((ct->numseq[i-1]==3&&ct->numseq[j+1]==4)||(ct->numseq[i-1]==4&&ct->numseq[j+1]==3)) {

								cumulative+=erg1(i,j,i+1,j-1,ct,data)*v->f(i+1,j-1)/denominator;
								if (cumulative>roll) {
									stack.push(i+1,j-1,1,0,0);
									found=true;

								}
				
							}

						}


					}
					for (d=(j-i-3);d>=1&&!found;d--) {
						for (ip=(i+1);ip<=(j-1-d)&&!found;ip++) {
            				jp = d+ip;
							if ((j-i-2-d)>(data->maxintloopsize)) break;
								if (abs(ip-i+jp-j)<=(data->maxintloopsize)) {

									cumulative+=erg2(i,j,ip,jp,ct,data,fce->f(i,ip),fce->f(jp,j))*v->f(ip,jp)/denominator;

									if (!found&&cumulative>roll) {

										stack.push(ip,jp,1,0,0);
										found=true;

									}

									if (!found&&(mod[ip]||mod[jp])) if (inc[ct->numseq[ip]][ct->numseq[jp]]&&notgu(ip,jp,ct)&&!(fce->f(ip,jp)&SINGLE)) {
										//i or j is modified
										cumulative+=erg2(i,j,ip,jp,ct,data,fce->f(i,ip),fce->f(jp,j))*
                  							v->f(ip+1,jp-1)*erg1(ip,jp,ip+1,jp-1,ct,data)/denominator;
										if (cumulative>roll) {

											regbp(ct,number,ip,jp);
											stack.push(ip+1,jp-1,1,0,0);
											found=true;
										}

									}
								}
						
						}
					}
					//consider the multiloop closed by i,j
					if (!found&&(j-i)>(2*minloop+4)) {
          				
						#ifdef SIMPLEMBLOOP
						//i+1 dangles
						cumulative+=erg4(i,j,i+1,1,ct,data,lfce[i+1])*penalty(i,j,ct,data)*
            				wmb->f(i+2,j-1)* data->eparam[5] * data->eparam[6] * data->eparam[10]*twoscaling/v->f(i,j);
						if (!found&&cumulative>roll) {
							stack.push(i+2,j-1,3,0,0);
							found=true;
			
						}
						#else  //!SIMPLEMBLOOP

						
						//no dangling ends on i-j pair:
						cumulative+=wmb->f(i+1,j-1)*data->eparam[5]*data->eparam[10]
            				*penalty(i,j,ct,data)*twoscaling/denominator;

						if (cumulative>roll) {
							stack.push(i+1,j-1,3,0,0);
							found=true;
			
						}

						//i+1 dangles
						cumulative+=erg4(i,j,i+1,1,ct,data,lfce[i+1])*penalty(i,j,ct,data)*
            				wmb->f(i+2,j-1)* data->eparam[5] * data->eparam[6] * data->eparam[10]*twoscaling/denominator;
						if (!found&&cumulative>roll) {
							stack.push(i+2,j-1,3,0,0);
							found=true;
			
						}

						//j-1 dangles
						cumulative+=erg4(i,j,j-1,2,ct,data,lfce[j-1]) * penalty(i,j,ct,data) *
            				wmb->f(i+1,j-2) * data->eparam[5] * data->eparam[6] * data->eparam[10]*twoscaling/denominator;
						if (!found&&cumulative>roll) {
							stack.push(i+1,j-2,3,0,0);
							found=true;
			
						}

						//both i+1 and j-1 dangle
				
            			cumulative+=data->tstkm[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]]*
								pfchecknp(lfce[i+1],lfce[j-1])*
								wmb->f(i+2,j-2) * data->eparam[5] * data->eparam[6] * data->eparam[6]* data->eparam[10]
								*penalty(i,j,ct,data)*twoscaling/denominator;

						if (!found&&cumulative>roll) {
							stack.push(i+2,j-2,3,0,0);
							found=true;
			
						}

						#ifndef disablecoax //a flag to turn off coaxial stacking

						//consider the coaxial stacking of a helix from i to j onto helix i+1 or i+2 to ip:
						for (ip=i+1;(ip<j)&&!found;ip++) {
							
							//first consider flush stacking
							cumulative+=penalty(i,j,ct,data)*v->f(i+1,ip)*
								penalty(i+1,ip,ct,data)*data->eparam[5]
								*data->eparam[10]*data->eparam[10]*(w->f(ip+1,j-1))*ergcoaxflushbases(j,i,i+1,ip,ct,data)
								*twoscaling/denominator;

							if (!found&&cumulative>roll) {
								stack.push(ip+1,j-1,4,0,0);
								stack.push(i+1,ip,1,0,0);
								found=true;
			
							}

							cumulative+=penalty(i,j,ct,data)*v->f(i+1,ip)*
								penalty(i+1,ip,ct,data)*data->eparam[5]
								*data->eparam[10]*data->eparam[10]*(wmb->f(ip+1,j-1))*ergcoaxflushbases(j,i,i+1,ip,ct,data)
								*twoscaling/denominator;

							if (!found&&cumulative>roll) {
								stack.push(ip+1,j-1,3,0,0);
								stack.push(i+1,ip,1,0,0);
								found=true;
			
							}

							if((mod[i+1]||mod[ip])) if (inc[ct->numseq[i+2]][ct->numseq[ip-1]]&&notgu(i+1,ip,ct)&&!(fce->f(i+1,ip)&SINGLE)) {

								cumulative+=penalty(i,j,ct,data)*v->f(i+2,ip-1)*
									penalty(i+1,ip,ct,data)*data->eparam[5]
									*data->eparam[10]*data->eparam[10]*(w->f(ip+1,j-1))*ergcoaxflushbases(j,i,i+1,ip,ct,data)
									*erg1(i+1,ip,i+2,ip-1,ct,data)*twoscaling/denominator;

								if (!found&&cumulative>roll) {
									stack.push(ip+1,j-1,4,0,0);
									stack.push(i+2,ip-1,1,0,0);
									regbp(ct,number,i+1,ip);
									found=true;
			
								}



								cumulative+=penalty(i,j,ct,data)*v->f(i+2,ip-1)*
									penalty(i+1,ip,ct,data)*data->eparam[5]
									*data->eparam[10]*data->eparam[10]*(wmb->f(ip+1,j-1))*ergcoaxflushbases(j,i,i+1,ip,ct,data)
									*erg1(i+1,ip,i+2,ip-1,ct,data)*twoscaling/denominator;

								if (!found&&cumulative>roll) {
									stack.push(ip+1,j-1,3,0,0);
									stack.push(i+2,ip-1,1,0,0);
									regbp(ct,number,i+1,ip);
									found=true;
			
								}

							}

							//Now calculate ca stacki8ng with intervening mismatch
							if ((ip+2<j-1)) {
								cumulative+=penalty(i,j,ct,data)*v->f(i+2,ip)*
									penalty(i+2,ip,ct,data)*data->eparam[5]
									*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
									(w->f(ip+2,j-1))
									*ergcoaxinterbases2(j,i,i+2,ip,ct,data)*twoscaling*pfchecknp(lfce[i+1],lfce[ip+1])/denominator;

								if (!found&&cumulative>roll) {
									stack.push(ip+2,j-1,4,0,0);
									stack.push(i+2,ip,1,0,0);
									
									found=true;
			
								}


								cumulative+=penalty(i,j,ct,data)*v->f(i+2,ip)*
									penalty(i+2,ip,ct,data)*data->eparam[5]
									*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
									(wmb->f(ip+2,j-1))
									*ergcoaxinterbases2(j,i,i+2,ip,ct,data)*twoscaling*pfchecknp(lfce[i+1],lfce[ip+1])/denominator;

								if (!found&&cumulative>roll) {
									stack.push(ip+2,j-1,3,0,0);
									stack.push(i+2,ip,1,0,0);
									
									found=true;
			
								}



							}

							if((mod[i+2]||mod[ip])) if (inc[ct->numseq[i+3]][ct->numseq[ip-1]]&&notgu(i+2,ip,ct)
								&&!(fce->f(i+2,ip)&SINGLE)) {

								cumulative+=penalty(i,j,ct,data)*v->f(i+3,ip-1)*
									penalty(i+2,ip,ct,data)*data->eparam[5]
									*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
									(w->f(ip+2,j-1))
									*ergcoaxinterbases2(j,i,i+2,ip,ct,data)
									*erg1(i+2,ip,i+3,ip-1,ct,data)*twoscaling*pfchecknp(lfce[i+1],lfce[ip+1])/denominator;

								if (!found&&cumulative>roll) {
									stack.push(ip+2,j-1,4,0,0);
									stack.push(i+3,ip-1,1,0,0);
									regbp(ct,number,i+2,ip);
									
									found=true;
			
								}


								cumulative+=penalty(i,j,ct,data)*v->f(i+3,ip-1)*
									penalty(i+2,ip,ct,data)*data->eparam[5]
									*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
									(wmb->f(ip+2,j-1))
									*ergcoaxinterbases2(j,i,i+2,ip,ct,data)
									*erg1(i+2,ip,i+3,ip-1,ct,data)*twoscaling*pfchecknp(lfce[i+1],lfce[ip+1])/denominator;

								if (!found&&cumulative>roll) {
									stack.push(ip+2,j-1,3,0,0);
									stack.push(i+3,ip-1,1,0,0);
									regbp(ct,number,i+2,ip);
									
									found=true;
			
								}

							}

							if (ip+1<j-2) {
								cumulative+=penalty(i,j,ct,data)*v->f(i+2,ip)*
									penalty(i+2,ip,ct,data)*data->eparam[5]
									*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
									(w->f(ip+1,j-2))
									*ergcoaxinterbases1(j,i,i+2,ip,ct,data)*twoscaling
									*pfchecknp(lfce[i+1],lfce[j-1])/denominator;

								if (!found&&cumulative>roll) {
									stack.push(ip+1,j-2,4,0,0);
									stack.push(i+2,ip,1,0,0);
									
									found=true;
			
								}


								cumulative+=penalty(i,j,ct,data)*v->f(i+2,ip)*
									penalty(i+2,ip,ct,data)*data->eparam[5]
									*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
									(wmb->f(ip+1,j-2))
									*ergcoaxinterbases1(j,i,i+2,ip,ct,data)*twoscaling
									*pfchecknp(lfce[i+1],lfce[j-1])/denominator;

								if (!found&&cumulative>roll) {
									stack.push(ip+1,j-2,3,0,0);
									stack.push(i+2,ip,1,0,0);
									
									found=true;
			
								}


								if(!found&&(mod[i+2]||mod[ip])) if (inc[ct->numseq[i+3]][ct->numseq[ip-1]]&&notgu(i+2,ip,ct)
									&&!(fce->f(i+2,ip)&SINGLE)) {

									cumulative+=penalty(i,j,ct,data)*v->f(i+3,ip-1)*
										penalty(i+2,ip,ct,data)*data->eparam[5]
										*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
										(w->f(ip+1,j-2))
										*ergcoaxinterbases1(j,i,i+2,ip,ct,data)
										*erg1(i+2,ip,i+3,ip-1,ct,data)*twoscaling*pfchecknp(lfce[i+1],lfce[j-1])/denominator;

									if (!found&&cumulative>roll) {
										stack.push(ip+1,j-2,4,0,0);
										stack.push(i+3,ip-1,1,0,0);
										regbp(ct,number,i+2,ip);
										found=true;
			
									}

									cumulative+=penalty(i,j,ct,data)*v->f(i+3,ip-1)*
										penalty(i+2,ip,ct,data)*data->eparam[5]
										*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
										(wmb->f(ip+1,j-2))
										*ergcoaxinterbases1(j,i,i+2,ip,ct,data)
										*erg1(i+2,ip,i+3,ip-1,ct,data)*twoscaling*pfchecknp(lfce[i+1],lfce[j-1])/denominator;
									if (!found&&cumulative>roll) {
										stack.push(ip+1,j-2,3,0,0);
										stack.push(i+3,ip-1,1,0,0);
										regbp(ct,number,i+2,ip);
										found=true;
			
									}

								}
							}
						}
						//consider the coaxial stacking of a helix from i to j onto helix ip to j-2 or j-1:
						for (ip=j-1;ip>i&&!found;ip--) {
				
				
							//first consider flush stacking
							cumulative+=penalty(i,j,ct,data)*v->f(ip,j-1)*
								penalty(j-1,ip,ct,data)*data->eparam[5]
								*data->eparam[10]*data->eparam[10]*(w->f(i+1,ip-1))*ergcoaxflushbases(ip,j-1,j,i,ct,data)
								*twoscaling/denominator;
									
							if (!found&&cumulative>roll) {
								stack.push(ip,j-1,1,0,0);
								stack.push(i+1,ip-1,4,0,0);
								found=true;
							}


							cumulative+=penalty(i,j,ct,data)*v->f(ip,j-1)*
								penalty(j-1,ip,ct,data)*data->eparam[5]
								*data->eparam[10]*data->eparam[10]*(wmb->f(i+1,ip-1))*ergcoaxflushbases(ip,j-1,j,i,ct,data)
								*twoscaling/denominator;

							if (!found&&cumulative>roll) {
								stack.push(ip,j-1,1,0,0);
								stack.push(i+1,ip-1,3,0,0);
								found=true;
							}

							if((mod[ip]||mod[j-1])) if(inc[ct->numseq[ip+1]][ct->numseq[j-2]]&&notgu(ip,j-1,ct)&&!(fce->f(ip,j-1)&SINGLE)) {
									
								cumulative+=penalty(i,j,ct,data)*v->f(ip+1,j-2)*
									penalty(j-1,ip,ct,data)*data->eparam[5]
									*data->eparam[10]*data->eparam[10]*(w->f(i+1,ip-1))
									*ergcoaxflushbases(ip,j-1,j,i,ct,data)
									*erg1(ip,j-1,ip+1,j-2,ct,data)*twoscaling/denominator;

								if (!found&&cumulative>roll) {
									stack.push(ip+1,j-2,1,0,0);
									stack.push(i+1,ip-1,4,0,0);
									regbp(ct,number,ip,j-1);
									found=true;
								}


								cumulative+=penalty(i,j,ct,data)*v->f(ip+1,j-2)*
									penalty(j-1,ip,ct,data)*data->eparam[5]
									*data->eparam[10]*data->eparam[10]*(wmb->f(i+1,ip-1))
									*ergcoaxflushbases(ip,j-1,j,i,ct,data)
									*erg1(ip,j-1,ip+1,j-2,ct,data)*twoscaling/denominator;
								if (!found&&cumulative>roll) {
									stack.push(ip+1,j-2,1,0,0);
									stack.push(i+1,ip-1,3,0,0);
									regbp(ct,number,ip,j-1);
									found=true;
								}

							}

							//now consider an intervening nuc
								
							if (ip-2>i+1) {
								cumulative+=penalty(i,j,ct,data)*v->f(ip,j-2)*
									penalty(j-2,ip,ct,data)*data->eparam[5]
									*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
									(w->f(i+1,ip-2))
									*ergcoaxinterbases1(ip,j-2,j,i,ct,data)*twoscaling*pfchecknp(lfce[j-1],lfce[ip-1])/denominator;

								if (!found&&cumulative>roll) {
									stack.push(ip,j-2,1,0,0);
									stack.push(i+1,ip-2,4,0,0);
									found=true;
								}

								cumulative+=penalty(i,j,ct,data)*v->f(ip,j-2)*
									penalty(j-2,ip,ct,data)*data->eparam[5]
									*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
									(wmb->f(i+1,ip-2))
									*ergcoaxinterbases1(ip,j-2,j,i,ct,data)*twoscaling*pfchecknp(lfce[j-1],lfce[ip-1])/denominator;

								if (!found&&cumulative>roll) {
									stack.push(ip,j-2,1,0,0);
									stack.push(i+1,ip-2,3,0,0);
									found=true;
								}

								if((mod[ip]||mod[j-2])) if(inc[ct->numseq[ip+1]][ct->numseq[j-3]]&&notgu(ip,j-2,ct)&&!(fce->f(ip,j-2)&SINGLE)) {
										
									cumulative+=penalty(i,j,ct,data)*v->f(ip+1,j-3)*
										penalty(j-2,ip,ct,data)*data->eparam[5]
										*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
										(w->f(i+1,ip-2))
										*ergcoaxinterbases1(ip,j-2,j,i,ct,data)
										*erg1(ip,j-2,ip+1,j-3,ct,data)*twoscaling*pfchecknp(lfce[j-1],lfce[ip-1])/denominator;

									if (!found&&cumulative>roll) {
										stack.push(ip+1,j-3,1,0,0);
										stack.push(i+1,ip-2,4,0,0);
										regbp(ct,number,ip,j-2);
										found=true;
									}


									cumulative+=penalty(i,j,ct,data)*v->f(ip+1,j-3)*
										penalty(j-2,ip,ct,data)*data->eparam[5]
										*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
										(wmb->f(i+1,ip-2))
										*ergcoaxinterbases1(ip,j-2,j,i,ct,data)
										*erg1(ip,j-2,ip+1,j-3,ct,data)*twoscaling*pfchecknp(lfce[j-1],lfce[ip-1])/denominator;

									if (!found&&cumulative>roll) {
										stack.push(ip+1,j-3,1,0,0);
										stack.push(i+1,ip-2,3,0,0);
										regbp(ct,number,ip,j-2);
										found=true;
									}

								}

							}

							if ((ip-1>i+2)) {
								cumulative+=penalty(i,j,ct,data)*v->f(ip,j-2)*
									penalty(j-2,ip,ct,data)*data->eparam[5]
									*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
									(w->f(i+2,ip-1))
									*ergcoaxinterbases2(ip,j-2,j,i,ct,data)*twoscaling*pfchecknp(lfce[j-1],lfce[i+1])/denominator;

								if (!found&&cumulative>roll) {
									stack.push(ip,j-2,1,0,0);
									stack.push(i+2,ip-1,4,0,0);
									found=true;

								}


								cumulative+=penalty(i,j,ct,data)*v->f(ip,j-2)*
									penalty(j-2,ip,ct,data)*data->eparam[5]
									*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
									(wmb->f(i+2,ip-1))
									*ergcoaxinterbases2(ip,j-2,j,i,ct,data)*twoscaling*pfchecknp(lfce[j-1],lfce[i+1])/denominator;

								if (!found&&cumulative>roll) {
									stack.push(ip,j-2,1,0,0);
									stack.push(i+2,ip-1,3,0,0);
									found=true;

								}

								if((mod[ip]||mod[j-2])) if(inc[ct->numseq[ip+1]][ct->numseq[j-3]]&&notgu(ip,j-2,ct)
									&&!(fce->f(ip,j-2)&SINGLE)) {
									cumulative+=penalty(i,j,ct,data)*v->f(ip+1,j-3)*
										penalty(j-2,ip,ct,data)*data->eparam[5]
										*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
										(w->f(i+2,ip-1))
										*ergcoaxinterbases2(ip,j-2,j,i,ct,data)
										*erg1(ip,j-2,ip+1,j-3,ct,data)*twoscaling*pfchecknp(lfce[j-1],lfce[i+1])/denominator;

									if (!found&&cumulative>roll) {
										stack.push(ip+1,j-3,1,0,0);
										stack.push(i+2,ip-1,4,0,0);
										regbp(ct,number,ip,j-2);
										found=true;

									}


									cumulative+=penalty(i,j,ct,data)*v->f(ip+1,j-3)*
										penalty(j-2,ip,ct,data)*data->eparam[5]
										*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
										(wmb->f(i+2,ip-1))
										*ergcoaxinterbases2(ip,j-2,j,i,ct,data)
										*erg1(ip,j-2,ip+1,j-3,ct,data)*twoscaling*pfchecknp(lfce[j-1],lfce[i+1])/denominator;

									if (!found&&cumulative>roll) {
										stack.push(ip+1,j-3,1,0,0);
										stack.push(i+2,ip-1,3,0,0);
										regbp(ct,number,ip,j-2);
										found=true;

									}

								}
							}


						
						}			
						#endif  //disable coax stacking
						#endif  //SIMPLEMBLOOP
						
					}

					
					if (!found) {
						cout << "Traceback error in v!\n";
						tracebackerror=14;
					}
		
					

					break;

				case 2: //switchcase = 2, dealing with a wcoax fragment
					#ifndef SIMPLEMBLOOP
					#ifndef disablecoax

					for (ip=i+minloop+1;ip<j-minloop-1&&!found;ip++) {
						
						//first consider flush stacking
						cumulative+=v->f(i,ip)*v->f(ip+1,j)*penalty(i,ip,ct,data)
							*penalty(ip+1,j,ct,data)*ergcoaxflushbases(i,ip,ip+1,j,ct,data)*data->eparam[10]*data->eparam[10]/wcoax->f(i,j);

						if (cumulative > roll) {
							stack.push(i,ip,1,0,0);
							stack.push(ip+1,j,1,0,0);
							found=true;
						}


						if ((mod[i]||mod[ip]||mod[ip+1]||mod[j])) {

							if ((mod[i]||mod[ip])&&(mod[ip+1]||mod[j])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]
								&&inc[ct->numseq[ip+2]][ct->numseq[j-1]]&&notgu(i,ip,ct)&&notgu(ip+1,j,ct)
								&&!(fce->f(ip+1,j)&SINGLE)&&!(fce->f(i,ip)&SINGLE)) {

								cumulative+=v->f(i+1,ip-1)*v->f(ip+2,j-1)*penalty(i,ip,ct,data)
									*penalty(ip+1,j,ct,data)*ergcoaxflushbases(i,ip,ip+1,j,ct,data)
									*erg1(i,ip,i+1,ip-1,ct,data)*erg1(ip+1,j,ip+2,j-1,ct,data)*data->eparam[10]*data->eparam[10]/wcoax->f(i,j);

								if (!found&&cumulative>roll) {
									stack.push(i+1,ip-1,1,0,0);
									stack.push(ip+2,j-1,1,0,0);
									found=true;
									regbp(ct,number,i,ip);
									regbp(ct,number,ip+1,j);
								}


							}

							if ((mod[i]||mod[ip])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]&&notgu(i,ip,ct)&&!(fce->f(i,ip)&SINGLE)) {
							
								cumulative+=v->f(i+1,ip-1)*v->f(ip+1,j)*penalty(i,ip,ct,data)
									*penalty(ip+1,j,ct,data)*ergcoaxflushbases(i,ip,ip+1,j,ct,data)
									*erg1(i,ip,i+1,ip-1,ct,data)*data->eparam[10]*data->eparam[10]/wcoax->f(i,j);

								if (!found&&cumulative>roll) {
									stack.push(i+1,ip-1,1,0,0);
									stack.push(ip+1,j,1,0,0);
									found=true;
									regbp(ct,number,i,ip);

								}


							}

							if ((mod[ip+1]||mod[j])&&inc[ct->numseq[ip+2]][ct->numseq[j-1]]&&notgu(ip+1,j,ct)&&!(fce->f(ip+1,j)&SINGLE)) {


								cumulative+=v->f(i,ip)*v->f(ip+2,j-1)*penalty(i,ip,ct,data)
									*penalty(ip+1,j,ct,data)*ergcoaxflushbases(i,ip,ip+1,j,ct,data)
									*erg1(ip+1,j,ip+2,j-1,ct,data)*data->eparam[10]*data->eparam[10]/wcoax->f(i,j);

								if (!found&&cumulative>roll) {
									stack.push(i,ip,1,0,0);
									stack.push(ip+2,j-1,1,0,0);
									found=true;
									regbp(ct,number,ip+1,j);
			
								}

							}


						}

						if (!lfce[ip+1]&&!lfce[j]) {
							//now consider an intervening mismatch
							cumulative+=v->f(i,ip)*v->f(ip+2,j-1)*penalty(i,ip,ct,data)
								*penalty(ip+2,j-1,ct,data)*ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)*data->eparam[10]*data->eparam[10]*data->eparam[6]*data->eparam[6]/wcoax->f(i,j);

							if (cumulative>roll&&!found) {
								stack.push(i,ip,1,0,0);
								stack.push(ip+2,j-1,1,0,0);
								found=true;

							}

							if (mod[i]||mod[ip]||mod[ip+2]||mod[j-1]) {
								if ((mod[i]||mod[ip])&&(mod[ip+2]||mod[j-1])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]
									&&inc[ct->numseq[ip+3]][ct->numseq[j-2]]&&notgu(i,ip,ct)&&notgu(ip+2,j-1,ct)
										&&!(fce->f(i,ip)&SINGLE)&&!(fce->f(ip+2,j-1)&SINGLE)) {

									 cumulative+=v->f(i+1,ip-1)*v->f(ip+3,j-2)*penalty(i,ip,ct,data)
										*penalty(ip+2,j-1,ct,data)*ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)
										*erg1(i,ip,i+1,ip-1,ct,data)*erg1(ip+2,j-1,ip+3,j-2,ct,data)*data->eparam[10]*data->eparam[10]*data->eparam[6]*data->eparam[6]/wcoax->f(i,j);

									 if (!found&&cumulative>roll) {
										stack.push(i+1,ip-1,1,0,0);
										stack.push(ip+3,j-2,1,0,0);
										regbp(ct,number,i,ip);
										regbp(ct,number,ip+2,j-1);
										found=true;

									 }


								}

								if ((mod[i]||mod[ip])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]&&notgu(i,ip,ct)&&!(fce->f(i,ip)&SINGLE)) {
							
									cumulative+=v->f(i+1,ip-1)*v->f(ip+2,j-1)*penalty(i,ip,ct,data)
										*penalty(ip+2,j-1,ct,data)*ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)
										*erg1(i,ip,i+1,ip-1,ct,data)*data->eparam[10]*data->eparam[10]*data->eparam[6]*data->eparam[6]/wcoax->f(i,j);

									if(!found&&cumulative>roll) {
										stack.push(i+1,ip-1,1,0,0);
										stack.push(ip+2,j-1,1,0,0);
										regbp(ct,number,i,ip);
										found=true;

									}


								}

								if ((mod[ip+2]||mod[j-1])&&inc[ct->numseq[ip+3]][ct->numseq[j-2]]&&notgu(ip+2,j-1,ct)&&!(fce->f(ip+2,j-1)&SINGLE)) {


									cumulative+=v->f(i,ip)*v->f(ip+3,j-2)*penalty(i,ip,ct,data)
										*penalty(ip+2,j-1,ct,data)*ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)
										*erg1(ip+2,j-1,ip+3,j-2,ct,data)*data->eparam[10]*data->eparam[10]*data->eparam[6]*data->eparam[6]/wcoax->f(i,j);

									if (!found&&cumulative>roll) {
										stack.push(i,ip,1,0,0);
										stack.push(ip+3,j-2,1,0,0);
										regbp(ct,number,ip+2,j-1);
										found=true;

									}

								}
							}
						}

						if(!lfce[i]&&!lfce[ip+1]) {
							cumulative+=v->f(i+1,ip)*v->f(ip+2,j)*penalty(i+1,ip,ct,data)
								*penalty(ip+2,j,ct,data)*ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)*data->eparam[10]*data->eparam[10]*data->eparam[6]*data->eparam[6]/wcoax->f(i,j);


							if (!found&&cumulative>roll) {
								stack.push(i+1,ip,1,0,0);
								stack.push(ip+2,j,1,0,0);
								
								found=true;

							}

							if (mod[i+1]||mod[ip]||mod[ip+2]||mod[j]) {
								if ((mod[i+1]||mod[ip])&&(mod[ip+2]||mod[j])&&inc[ct->numseq[i+2]][ct->numseq[ip-1]]
									&&inc[ct->numseq[ip+3]][ct->numseq[j-1]]&&notgu(i+1,ip,ct)&&notgu(ip+2,j,ct)
									&&!(fce->f(i+1,ip)&SINGLE)&&!(fce->f(ip+2,j)&SINGLE)	) {

									cumulative+=v->f(i+2,ip-1)*v->f(ip+3,j-1)*penalty(i+1,ip,ct,data)
										*penalty(ip+2,j,ct,data)*ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)
										*erg1(i+1,ip,i+2,ip-1,ct,data)*erg1(ip+2,j,ip+3,j-1,ct,data)*data->eparam[10]*data->eparam[10]*data->eparam[6]*data->eparam[6]/wcoax->f(i,j);	

									if (!found&&cumulative>roll) {
										stack.push(i+2,ip-1,1,0,0);
										stack.push(ip+3,j-1,1,0,0);
										regbp(ct,number,i+1,ip);
										regbp(ct,number,ip+2,j);
										found=true;

									}

							
								}
								if ((mod[i+1]||mod[ip])&&inc[ct->numseq[i+2]][ct->numseq[ip-1]]&&notgu(i+1,ip,ct)&&!(fce->f(i+1,ip)&SINGLE)) {
							
									cumulative+=v->f(i+2,ip-1)*v->f(ip+2,j)*penalty(i+1,ip,ct,data)
										*penalty(ip+2,j,ct,data)*ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)
										*erg1(i+1,ip,i+2,ip-1,ct,data)*data->eparam[10]*data->eparam[10]*data->eparam[6]*data->eparam[6]/wcoax->f(i,j);

									if (!found&&cumulative>roll) {
										stack.push(i+2,ip-1,1,0,0);
										stack.push(ip+2,j,1,0,0);
										regbp(ct,number,i+1,ip);
										
										found=true;

									}


								}

								if ((mod[ip+2]||mod[j])&&inc[ct->numseq[ip+3]][ct->numseq[j-1]]&&notgu(ip+2,j,ct)&&!(fce->f(ip+2,j)&SINGLE)) {


									cumulative+=v->f(i+1,ip)*v->f(ip+3,j-1)*penalty(i+1,ip,ct,data)
										*penalty(ip+2,j,ct,data)*ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)
										*erg1(ip+2,j,ip+3,j-1,ct,data)*data->eparam[10]*data->eparam[10]*data->eparam[6]*data->eparam[6]/wcoax->f(i,j);

									if (!found&&cumulative>roll) {
										stack.push(i+1,ip,1,0,0);
										stack.push(ip+3,j-1,1,0,0);
										
										regbp(ct,number,ip+2,j);
										found=true;

									}

								}
							}
						

						}
					}

					if (!found) {

						cout << "Traceback error at wcoax/n";
						tracebackerror=14;

					}

					#endif
					#endif //SIMPLEMBLOOP
					break;

				case 3: //switchcase = 3, dealing with wmb fragment

					cumulative = wmbl->f(i,j)/wmb->f(i,j);

					if (cumulative>roll) {
						found=true;
						stack.push(i,j,5,0,0);

					}
			
					if (!lfce[j]) {
						
						cumulative+=wmb->f(i,j-1)*data->eparam[6]*data->scaling/wmb->f(i,j);

						if (!found&&cumulative>roll) {
							found=true;
							stack.push(i,j-1,3,0,0);

						}
					

					}

					if (!found) {
						cout << "Traceback error at wmb\n";
						tracebackerror=14;

					}


					break;

				case 4: //switchcase = 4, dealing with a w fragment

					cumulative = wl->f(i,j)/w->f(i,j);
					if (cumulative>roll) {
						stack.push(i,j,6,0,0);
						found=true;

					}

					if (!lfce[j]) {
             	
               			cumulative+= w->f(i,j-1) * data->eparam[6]*data->scaling/w->f(i,j);
						if (cumulative>roll&&!found) {
							stack.push(i,j-1,4,0,0);
							found=true;

						}
               

					}

					if (!found) {
						cout << "Traceback error at w\n";
						tracebackerror=14;

					}
	  

					break;

				case 5:  //switchcase = 5, wmbl fragment
					
					cumulative=wcoax->f(i,j)/wmbl->f(i,j);
					if (cumulative>roll&&!found) {
						stack.push(i,j,2,0,0);
								
						found=true;

					}



					for (k=i+1;k<j&&!found;k++) {
					
						
						if (!lfce[i]) {
							cumulative+=(wl->f(i,k)-wl->f(i+1,k)*data->eparam[6]*data->scaling)*(wl->f(k+1,j))/wmbl->f(i,j);
							if (cumulative>roll&&!found) {
								stack.push(i,k,7,0,0);
								stack.push(k+1,j,6,0,0);
								found=true;

							}

							cumulative+=(wcoax->f(i,k))*(wl->f(k+1,j))/wmbl->f(i,j);
							if (cumulative>roll&&!found) {
								stack.push(i,k,2,0,0);
								stack.push(k+1,j,6,0,0);
								found=true;

							}

							cumulative+=(wl->f(i,k)-wl->f(i+1,k)*data->eparam[6]*data->scaling)*(wmbl->f(k+1,j))/wmbl->f(i,j);
							if (cumulative>roll&&!found) {
								stack.push(i,k,7,0,0);
								stack.push(k+1,j,5,0,0);
								found=true;

							}

							cumulative+=(wcoax->f(i,k))*(wmbl->f(k+1,j))/wmbl->f(i,j);
							if (cumulative>roll&&!found) {
								stack.push(i,k,2,0,0);
								stack.push(k+1,j,5,0,0);
								found=true;

							}

						}

						else {
							
							cumulative+=(wl->f(i,k))*(wl->f(k+1,j))/wmbl->f(i,j);
							if (cumulative>roll&&!found) {
								stack.push(i,k,6,0,0);
								stack.push(k+1,j,6,0,0);
								found=true;

							}

							cumulative+=(wcoax->f(i,k))*(wl->f(k+1,j))/wmbl->f(i,j);
							if (cumulative>roll&&!found) {
								stack.push(i,k,2,0,0);
								stack.push(k+1,j,6,0,0);
								found=true;

							}

							cumulative+=(wl->f(i,k)*wmbl->f(k+1,j))/wmbl->f(i,j);
							if (cumulative>roll&&!found) {
								stack.push(i,k,6,0,0);
								stack.push(k+1,j,5,0,0);
								found=true;

							}

							cumulative+=(wcoax->f(i,k))*(wmbl->f(k+1,j))/wmbl->f(i,j);
							if (cumulative>roll&&!found) {
								stack.push(i,k,2,0,0);
								stack.push(k+1,j,5,0,0);
								found=true;

							}


						}

						
					}


					
					if (!lfce[i]) {
						cumulative+=wmbl->f(i+1,j)*data->eparam[6]*data->scaling/wmbl->f(i,j);
						if (cumulative>roll&&!found) {
								
								stack.push(i+1,j,5,0,0);
								found=true;

							}

					}

					if (!found) {
						cout << "Traceback error at wmbl\n";
						tracebackerror=14;
			
					}

				break;

				case 6: //switchcase 6, wl fragment

					#ifdef SIMPLEMBLOOP
					if (!lfce[i]) {
         		
         				cumulative=  wl->f(i+1,j)*data->eparam[6]*data->scaling/wl->f(i,j);

						if (cumulative>roll) {
							stack.push(i+1,j,6,0,0);
							found=true;

						}
            	
					}
					cumulative+= v->f(i,j-1)* data->eparam[10] * data->eparam[6] *
         				erg4(j-1,i,j,1,ct,data,lfce[j])*penalty(i,j-1,ct,data)/wl->f(i,j);

					if (cumulative>roll&&!found) {
						stack.push(i,j-1,1,0,0);
						found=true;

					}

					if ((mod[i]||mod[j-1])) if(inc[ct->numseq[i+1]][ct->numseq[j-2]]&&notgu(i,j-1,ct)&&!(fce->f(i,j-1)&SINGLE)) {

						cumulative+= v->f(i+1,j-2) * data->eparam[10] * data->eparam[6] *
         					erg4(j-1,i,j,1,ct,data,lfce[j])*penalty(i,j-1,ct,data)
							*erg1(i,j-1,i+1,j-2,ct,data)/wl->f(i,j);

						if (cumulative>roll&&!found) {
							stack.push(i+1,j-2,1,0,0);
							regbp(ct,number,i,j-1);
							found=true;

						}
				


					}
					#else //!SIMPLEMBLOOP

					if (!lfce[i]) {
         		
         				cumulative=  wl->f(i+1,j)*data->eparam[6]*data->scaling/wl->f(i,j);

						if (cumulative>roll) {
							stack.push(i+1,j,6,0,0);
							found=true;

						}
            	
					}
					cumulative+= data->eparam[10]*v->f(i,j)*penalty(j,i,ct,data)/wl->f(i,j);
					if (cumulative>roll&&!found) {
						found=true;
						stack.push(i,j,1,0,0);

					}

					if ((mod[i]||mod[j])) if (inc[ct->numseq[i+1]][ct->numseq[j-1]]&&notgu(i,j,ct)) {

						cumulative+= data->eparam[10]*v->f(i+1,j-1)*penalty(j,i,ct,data)*erg1(i,j,i+1,j-1,ct,data)/wl->f(i,j);	

						if (!found&&cumulative>roll) {
							stack.push(i+1,j-1,1,0,0);
							regbp(ct,number,i,j);
							found = true;

						}

					} 

					cumulative+= v->f(i+1,j)*data->eparam[10]*data->eparam[6]*
         				erg4(j,i+1,i,2,ct,data,lfce[i])*penalty(i+1,j,ct,data)/wl->f(i,j);

					if (!found&&cumulative>roll) {
						stack.push(i+1,j,1,0,0);
						found=true;

					}

					if ((mod[i+1]||mod[j])) if(inc[ct->numseq[i+2]][ct->numseq[j-1]]&&notgu(i+1,j,ct)&&!(fce->f(i+1,j)&SINGLE)) {

						cumulative+= v->f(i+2,j-1) * data->eparam[10] *data->eparam[6] *
         					erg4(j,i+1,i,2,ct,data,lfce[i])*penalty(i+1,j,ct,data)
							*erg1(i+1,j,i+2,j-1,ct,data)/wl->f(i,j);

						if (!found&&cumulative>roll) {
							stack.push(i+2,j-1,1,0,0);
							found=true;
							regbp(ct,number,i+1,j);

						}


					}

					cumulative+= v->f(i,j-1)* data->eparam[10] * data->eparam[6] *
         				erg4(j-1,i,j,1,ct,data,lfce[j])*penalty(i,j-1,ct,data)/wl->f(i,j);

					if (cumulative>roll&&!found) {
						stack.push(i,j-1,1,0,0);
						found=true;

					}

					if ((mod[i]||mod[j-1])) if(inc[ct->numseq[i+1]][ct->numseq[j-2]]&&notgu(i,j-1,ct)&&!(fce->f(i,j-1)&SINGLE)) {

						cumulative+= v->f(i+1,j-2) * data->eparam[10] * data->eparam[6] *
         					erg4(j-1,i,j,1,ct,data,lfce[j])*penalty(i,j-1,ct,data)
							*erg1(i,j-1,i+1,j-2,ct,data)/wl->f(i,j);

						if (cumulative>roll&&!found) {
							stack.push(i+1,j-2,1,0,0);
							regbp(ct,number,i,j-1);
							found=true;

						}
				


					}

					if (j!=1&&!lfce[i]&&!lfce[j]) {
         				cumulative+= v->f(i+1,j-1) *data->eparam[10] * (data->eparam[6]*data->eparam[6]) *
         				data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
						*penalty(j-1,i+1,ct,data)/wl->f(i,j);

						if (cumulative>roll&&!found) {
							stack.push(i+1,j-1,1,0,0);
							found=true;

						}



						if ((mod[i+1]||mod[j-1])) if((j-2>0)&&!(fce->f(i+1,j-1)&SINGLE)) {
							if(inc[ct->numseq[i+2]][ct->numseq[j-2]]&&notgu(i+1,j-1,ct)) {

								cumulative+= v->f(i+2,j-2) * data->eparam[10] * (data->eparam[6]*data->eparam[6]) *
         							data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
									*penalty(j-1,i+1,ct,data)*erg1(i+1,j-1,i+2,j-2,ct,data)/wl->f(i,j);

								if (cumulative>roll&&!found) {
									stack.push(i+2,j-2,1,0,0);
									regbp(ct,number,i+1,j-1);
									found=true;

								}

							}
						}
					}
					#endif  //SIMPLEMBLOOP
					if (!found) {
						cout << "Traceback errors at wl\n";
						tracebackerror=14;

					}

				break;

				case 7: //switchcase 7, a wl seed (helix only)
						
					#ifdef SIMPLEMBLOOP

					cumulative+= v->f(i,j-1)* data->eparam[10] * data->eparam[6] *
         				erg4(j-1,i,j,1,ct,data,lfce[j])*penalty(i,j-1,ct,data)/(wl->f(i,j)-wl->f(i+1,j)*data->eparam[6]*data->scaling);

					if (cumulative>roll&&!found) {
						stack.push(i,j-1,1,0,0);
						found=true;

					}

					if ((mod[i]||mod[j-1])) if(inc[ct->numseq[i+1]][ct->numseq[j-2]]&&notgu(i,j-1,ct)&&!(fce->f(i,j-1)&SINGLE)) {

						cumulative+= v->f(i+1,j-2) * data->eparam[10] * data->eparam[6] *
         					erg4(j-1,i,j,1,ct,data,lfce[j])*penalty(i,j-1,ct,data)
							*erg1(i,j-1,i+1,j-2,ct,data)/(wl->f(i,j)-wl->f(i+1,j)*data->eparam[6]*data->scaling);

						if (cumulative>roll&&!found) {
							stack.push(i+1,j-2,1,0,0);
							regbp(ct,number,i,j-1);
							found=true;

						}
				


					}
					#else   //!SIMPLEMBLOOP

					cumulative+= data->eparam[10]*v->f(i,j)*penalty(j,i,ct,data)/(wl->f(i,j)-wl->f(i+1,j)*data->eparam[6]*data->scaling);
					if (cumulative>roll&&!found) {
						found=true;
						stack.push(i,j,1,0,0);

					}

					if ((mod[i]||mod[j])) if (inc[ct->numseq[i+1]][ct->numseq[j-1]]&&notgu(i,j,ct)) {

						cumulative+= data->eparam[10]*v->f(i+1,j-1)*penalty(j,i,ct,data)*erg1(i,j,i+1,j-1,ct,data)/(wl->f(i,j)-wl->f(i+1,j)*data->eparam[6]*data->scaling);	

						if (!found&&cumulative>roll) {
							stack.push(i+1,j-1,1,0,0);
							regbp(ct,number,i,j);
							found = true;

						}

					} 

					cumulative+= v->f(i+1,j)*data->eparam[10]*data->eparam[6]*
         				erg4(j,i+1,i,2,ct,data,lfce[i])*penalty(i+1,j,ct,data)/(wl->f(i,j)-wl->f(i+1,j)*data->eparam[6]*data->scaling);

					if (!found&&cumulative>roll) {
						stack.push(i+1,j,1,0,0);
						found=true;

					}

					if ((mod[i+1]||mod[j])) if(inc[ct->numseq[i+2]][ct->numseq[j-1]]&&notgu(i+1,j,ct)&&!(fce->f(i+1,j)&SINGLE)) {

						cumulative+= v->f(i+2,j-1) * data->eparam[10] *data->eparam[6] *
         					erg4(j,i+1,i,2,ct,data,lfce[i])*penalty(i+1,j,ct,data)
							*erg1(i+1,j,i+2,j-1,ct,data)/(wl->f(i,j)-wl->f(i+1,j)*data->eparam[6]*data->scaling);

						if (!found&&cumulative>roll) {
							stack.push(i+2,j-1,1,0,0);
							found=true;
							regbp(ct,number,i+1,j);

						}


					}

					cumulative+= v->f(i,j-1)* data->eparam[10] * data->eparam[6] *
         				erg4(j-1,i,j,1,ct,data,lfce[j])*penalty(i,j-1,ct,data)/(wl->f(i,j)-wl->f(i+1,j)*data->eparam[6]*data->scaling);

					if (cumulative>roll&&!found) {
						stack.push(i,j-1,1,0,0);
						found=true;

					}

					if ((mod[i]||mod[j-1])) if(inc[ct->numseq[i+1]][ct->numseq[j-2]]&&notgu(i,j-1,ct)&&!(fce->f(i,j-1)&SINGLE)) {

						cumulative+= v->f(i+1,j-2) * data->eparam[10] * data->eparam[6] *
         					erg4(j-1,i,j,1,ct,data,lfce[j])*penalty(i,j-1,ct,data)
							*erg1(i,j-1,i+1,j-2,ct,data)/(wl->f(i,j)-wl->f(i+1,j)*data->eparam[6]*data->scaling);

						if (cumulative>roll&&!found) {
							stack.push(i+1,j-2,1,0,0);
							regbp(ct,number,i,j-1);
							found=true;

						}
				


					}

					if (j!=1&&!lfce[i]&&!lfce[j]) {
         				cumulative+= v->f(i+1,j-1) *data->eparam[10] * (data->eparam[6]*data->eparam[6]) *
         				data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
						*penalty(j-1,i+1,ct,data)/(wl->f(i,j)-wl->f(i+1,j)*data->eparam[6]*data->scaling);

						if (cumulative>roll&&!found) {
							stack.push(i+1,j-1,1,0,0);
							found=true;

						}



						if ((mod[i+1]||mod[j-1])) if((j-2>0)&&!(fce->f(i+1,j-1)&SINGLE)) {
							if(inc[ct->numseq[i+2]][ct->numseq[j-2]]&&notgu(i+1,j-1,ct)) {

								cumulative+= v->f(i+2,j-2) * data->eparam[10] * (data->eparam[6]*data->eparam[6]) *
         							data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
									*penalty(j-1,i+1,ct,data)*erg1(i+1,j-1,i+2,j-2,ct,data)/(wl->f(i,j)-wl->f(i+1,j)*data->eparam[6]*data->scaling);

								if (cumulative>roll&&!found) {
									stack.push(i+2,j-2,1,0,0);
									regbp(ct,number,i+1,j-1);
									found=true;

								}

							}
						}
					}
					#endif //SIMPLEMBLOOP
					if (!found) {
						cout << "Traceback error at wl seed\n";
						tracebackerror=14;

					}


				break;
			}
			if (cumulative>(1.0+1e-5)) {
				cout << "Over 1 probability error\n";
				tracebackerror=21;

			}
			if (!found) {
				cout << "Overall traceback error\n";
				tracebackerror=14;

			}
		}

	}

	return tracebackerror;

}


//Stochastic() starts stochastic traceback from the point of reading a save file from disk.
//Stochastic() calls stochastictraceback() to do the actual work of sampling.
void stochastic(structure *ct, char *savefilename, int numberofstructures, int randomseed, ProgressHandler *progress) {
	
	
	short vers;

	DynProgArray<PFPRECISION> *w,*wmb,*wmbl,*wcoax,*wl,*v;
	forceclass *fce;
	PFPRECISION *w3,*w5,scaling;
	bool *lfce,*mod;
	pfdatatable *data;
	datatable *data2;
	

	ifstream sav(savefilename,ios::binary);

	read(&sav,&vers);//read the save file version
		//right now there is no infrastructure to indicate an error in file version
		//that should be added at some point
	
	int sequencelength;
	read(&sav,&(sequencelength));

	sav.close();
	//allocate everything
		
	data = new pfdatatable;
	data2 = new datatable;
	
	ct->allocate(sequencelength);

	w = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
	v = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
	wmb = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
	fce = new forceclass(ct->GetSequenceLength());
	wl = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
	wcoax = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
	wmbl = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());

	w5 = new PFPRECISION [ct->GetSequenceLength()+1];
	w3 = new PFPRECISION [ct->GetSequenceLength()+2];

	lfce = new bool [2*ct->GetSequenceLength()+1];
	mod = new bool [2*ct->GetSequenceLength()+1];

	
	

	//load all the data from the pfsavefile:
	readpfsave(savefilename, ct, w5, w3,v, w, wmb,wl, wmbl, wcoax, fce,&scaling,mod,lfce,data,data2);

	data->scaling = scaling;
	
	//now that partition function data was read from disk, do the sampling:
	stochastictraceback(w,wmb,wmbl,wcoax,wl,v,
		fce, w3,w5,scaling, lfce, mod, data, numberofstructures, 
		ct, randomseed, progress);


	//delete everything
	delete data;
	delete data2;
	delete w;
	delete v;
	delete wmb;
	delete fce;
	delete wl;
	delete wcoax;
	delete wmbl;
	delete[] w5;
	delete[] w3;
	delete[] lfce;
	delete[] mod;


}
