
#include "forceclass.h"

// forceclass encapsulates a large 2-d arrays of char used by the
// dynamic algorithm to enforce folding constraints

forceclass::forceclass(int size) {
	

  Size = size;
  register int i,j;
  dg = new char *[size+1];

	for (i=0;i<=(size);i++)  {
    dg[i] = new char [size+1];
  }
  for (i=0;i<=size;i++) {
    for (j=0;j<size+1;j++) {
      dg[i][j] = 0;
             
    }
  }
}

forceclass::~forceclass() {
  for (int i = 0; i <= Size; i++) {
    delete[] dg[i];
  }

  delete[] dg;
}
