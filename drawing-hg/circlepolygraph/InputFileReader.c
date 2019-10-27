/*
  INPUTFILEREADER.C

*/


#include "CirclePolyGraphHeaderFile.h" // File with important definitions



/* ******************************************************************************** */
void ReadInputFile(int *K, int *nStrands, int *nBases, int **SeqLen, int **Perm,
		   int **pairs, int *nPairs, char *InputFile) {

  /*
  */
  
  int i,j; // Counters
  char line[MAXLINE]; // A line in the input file
  char linecpy[MAXLINE]; // A copy of a line in the input file
  char *tok; // Token
  char tokseps[] = " ,;\t\n"; // Token separators 
  char  **Sequences; // The sequences of the strands
  FILE *fp;

  // Open the input file
  if ((fp = fopen(InputFile,"r")) == NULL) {
    printf("Error opening %s!\n\nExiting....\n",InputFile);
    exit(ERR_INPUTFILE);
  }

  // Blow through comments
  if (fgets(line,MAXLINE,fp) != NULL) {
    while (line[0] == '%' || line[0] == '\n' || line[0] == '\0') {
      fgets(line,MAXLINE,fp);
    }
  }
  else {
    printf("Error in input file [comments] %s!\n",InputFile);
    printf("\nExiting....\n\n");
    exit(ERR_INPUT);
  }

  // The first line is the number of strands, K
  // Make sure it's number of strands
  if (!isdigit(line[0])) {
    printf("Error in input file [num strands] %s!\n",InputFile);
    printf("\nExiting...\n\n");
    exit(ERR_INPUT);
  }
  *K = atoi(line);

  // Allocate sequences array
  Sequences = (char **) malloc((*K) * sizeof(char *));
  (*SeqLen) = (int *) malloc((*K) * sizeof(int));

  // The next K lines are the sequences
  for (i = 0; i < *K; i++) {
    if (fgets(line,MAXLINE,fp) == NULL) {
      printf("Error in input file [sequences] %s!\n",InputFile);
      printf("\nExiting...\n\n");
      exit(ERR_INPUT);
    }
	
    // Make sure line contains a sequence (check first character)
    if (isalpha(line[0]) == 0) {
      printf("Error in input file [sequence]%s!\n",InputFile);
      printf("\nExiting...\n\n");
      exit(ERR_INPUT);
    } 
    
    (*SeqLen)[i] = strlen(line);  // sequence length
    Sequences[i] = (char *) malloc(((*SeqLen)[i]+1) * sizeof(char));
    
    strcpy(Sequences[i],line);
    
    // Chop off null character if it's there
    if (Sequences[i][(*SeqLen)[i]-1] == '\n') {
      Sequences[i][(*SeqLen)[i]-1] = '\0';
      ((*SeqLen)[i])--;
    }
	}
  
  // Get the permuation
  if (fgets(line,MAXLINE,fp) == NULL) {
    printf("Error in input file [permutation] %s!\n",InputFile);
    printf("\nExiting...\n\n");
    exit(ERR_INPUT);
  }
  // Make a copy of the line
  strcpy(linecpy,line);

  // Count how many strands in the complex
  tok = strtok(line,tokseps);
  *nStrands = 1;
  while ((tok = strtok(NULL,tokseps)) != NULL) {
    (*nStrands)++;
  }

  // Allocate memory for the ordered complex
  (*Perm) = (int *) malloc(*nStrands * sizeof(int));

  // Now put the entries into the perm
  tok = strtok(linecpy,tokseps);
  i = 0;
  (*Perm)[i++] = atoi(tok) - 1;
  while ((tok = strtok(NULL,tokseps)) != NULL) {
    (*Perm)[i++] = atoi(tok) - 1;
  }

  // Get total number of bases
  *nBases = 0;
  for (i = 0; i < *nStrands; i++) {
    *nBases += (*SeqLen)[(*Perm)[i]];
  }

  // Allocate memory for pairs
  *pairs = (int *) malloc((*nBases) * sizeof(int));

  // Now read in the base pairs
  if (fgets(line,MAXLINE,fp) == NULL) { // The means there are no base pairs
    for (i = 0; i < *nBases; i++) {
      (*pairs)[i] = -1;
    }
  }
  else if (isdigit(line[0]) == 0) { // Dot-paren format
    getStructureFromParens(line,(*pairs),*nBases);
  }
  else { // Pair list format
    // Initialize pairs
    for (i = 0; i < *nBases; i++) {
      (*pairs)[i] = -1;
    }
    tok = strtok(line,tokseps);
    i = atoi(tok)-1;
    tok = strtok(NULL,tokseps);
    j = atoi(tok)-1;
    (*pairs)[i] = j;
    (*pairs)[j] = i;
    while (fgets(line,MAXLINE,fp) != NULL) {
      if (line[0] != '%' && line[0] != '\n' && line[0] != '\0') {
	tok = strtok(line,tokseps);
	i = atoi(tok)-1;
	tok = strtok(NULL,tokseps);
	j = atoi(tok)-1;
	(*pairs)[i] = j;
	(*pairs)[j] = i;
      }
    }
  }

  fclose(fp);


  // Get number of base pairs
  *nPairs = 0;
  for (i = 0; i < *nBases; i++) {
    if ((*pairs)[i] > -1) {
      (*nPairs)++;
    }
  }
  (*nPairs) /= 2;

  for (i = 0; i < *K; i++) {
    free(Sequences[i]);
  }
  free(Sequences);
  
}
/* ******************************************************************************** */




/* ******************************************************************************** */
void getStructureFromParens( char *line, int *pairs, int seqlength) {
  
  int i, j;
  int braces[seqlength];
  int leftParenIndex;
  
  char pairSymbols[] = { '(', ')', '{','}', '[', ']', '<', '>' };
  int type = 0;
  int nTypes = 4;
  
  for( i = 0; i <= seqlength-1; i++) {
    pairs[i] = -1;
  }
  
  
  for( type = 0; type < nTypes; type++) {
    leftParenIndex = 0;		
    for( i = 0; i <= seqlength-1; i++) {
      braces[i] = -5;
    }
    
    i = 0;
    j = 0;
    
    while( i <= seqlength - 1) {
      if( leftParenIndex < 0 || leftParenIndex >= seqlength) {
	printf("Too many %c, not enough %c!\n", pairSymbols[2*type+1],
	       pairSymbols[2*type]);
	exit(ERR_INVALIDSTRUCTURE);
      }
      if( line[j] == pairSymbols[ 2*type]) {
	braces[ leftParenIndex++] = i;
      } 
      else if( line[j] == pairSymbols[ 2*type+1]) {
	pairs[ braces[ --leftParenIndex]] = i;
	pairs[ i] = braces[ leftParenIndex];
      }
      j++;
      if( line[j] != '+') {
	i++;
      }
    }
  }
}
/* ******************************************************************************** */


