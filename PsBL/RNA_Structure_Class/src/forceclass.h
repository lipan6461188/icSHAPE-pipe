
#ifndef FORCECLASS_H
#define FORCECLASS_H

// forceclass encapsulates a large 2-d arrays used by the dynamic
// algorithm to enforce folding constraints
class forceclass {
private:
  int Size;

public:
  int k;
  char **dg;
  
  // the constructor allocates the space needed by the arrays
  forceclass(int size);

  // the destructor deallocates the space used
  ~forceclass();
  
  // f is an integer function that references the correct element of
  // the array
  char &f(int i, int j);
};

inline char &forceclass::f(int i, int j) {
      	
  if (i > j) {
    // swap them
    int temp = i;
    i = j;
    j = temp;
  }

  if (i > Size) {
    i -= Size;
    j -= Size;
  }

  return dg[i][j-i];
}

#endif
