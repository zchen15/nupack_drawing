/*
  UTILS.C 

  Utility files for array and string processing.  They do various
  matrix operations, dot products, etc.  A useful toolchest to have
  around.

  For use with Concentations.c.

  Justin Bois, Caltech, 1 September 2006
  bois@caltech.edu
*/

#include "CirclePolyGraphHeaderFile.h" // Concentrations header file.



/* ******************************************************************************** */
double str2double (char *str) {
  /* 
     Converts a string to a double.  The string may either be
     an integer or float.  E.g., 56 or 0.45 or 654.234.  It may 
     also be in scientific notation, e.g., 1.5e-12, 43.54e5,
     3.32E-8, 4E4, 2.45e+15, or 8.55E+12.  No spaces are allowed.
     This is easier than having a bunch of sscanf's with if state-
     ments in the main code.
   */

  int i,k; // counters
  int noE; // Haven't encountered an e or E yet.
  char *MantissaStr; // string storing the mantissa
  char *ExpStr; // string storing the exponent
  int Len; // length of string
  double mantissa; // number is mantissa * 10^exponent 
  int exponent; 

  Len = strlen(str);

  noE = 1;
  k = 0;
  while (k < Len && noE) {
    if (str[k] == 'e' || str[k] == 'E') {
      noE = 0;
    }
    k++;
  }

  if (k == Len) { // Not in scientific notation
    return atof(str);
  }

  // k is now the index of the start of the exponent
  ExpStr = (char *) malloc((Len-k+1) * sizeof(char));
  MantissaStr = (char *) malloc(k * sizeof(char));
  strncpy(MantissaStr,str,k-1);
  MantissaStr[k-1] = '\0';

  for (i = 0; i < Len-k; i++) {
    ExpStr[i] = str[k+i];
  }
  ExpStr[Len-k] = '\0';

  mantissa = atof(MantissaStr);
  exponent = atoi(ExpStr);

  free(MantissaStr);
  free(ExpStr);
  
  return (mantissa * pow(10,exponent));

}
/* ******************************************************************************** */

/* ******************************************************************************** */
void MatrixVectorMult(double *c, double **A, double *b, int n) {
  /*
    Performs the multiplication of the n x n matrix A with n-vector b
    and stores the result in vector c, which must be pre-allocated.

    All entries are doubles.
  */

  int i; // Counter

  for (i = 0; i < n; i++) {
      c[i] = dot(A[i],b,n);
  }

}
/* ******************************************************************************** */

/* ******************************************************************************** */
double dot(double *v1, double *v2, int len) {
  /*
    Computes dot product of v1 and v2 (v1^T v2) and
    returns the result.

    v1 and v2 must be doubles.  They must have len entries in them.
  */
  
  int i; // Counter
  double dotprod; // The dot product
  
  dotprod = 0.0;
  for (i = 0; i < len; i++) {
    dotprod += v1[i]*v2[i];
  }

  return dotprod;
}
/* ******************************************************************************** */

