#if !defined(BOLTZMANN_H)
#define BOLTZMANN_H

/*======================================================================
// define boltzman() *once* in the code -- this is the place right now. 
// is partition.h is a better place?
//
// also removed a version of boltzman that accepts ints -- this
//  causes weird results, e.g., in calculating partition function
//  and base pair probabilities with SHAPE or 'EX' data. (rhiju das, 2011)
//
========================================================================*/

inline PFPRECISION boltzman(double i, PFPRECISION temp) {

	if (i==INFINITE_ENERGY) return 0;
	else return exp((-((PFPRECISION) i)/((PFPRECISION)conversionfactor))/(RKC*temp));

}
#endif

