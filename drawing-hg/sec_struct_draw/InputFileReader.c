/*
  INPUTFILEREADER.C

*/


#include "SecStructDrawHeader.h" // File with important definitions



/* ******************************************************************************** */
void ReadInputFile( struct base **bases, int *nbases, char *InputFile, int *prob, 
                   char *probFile) {
 /*
   Takes data from input file and converts it to the bases struct.
     */
 
 int i,j,k; // Counters
 int nSeqs; // Number of sequences
 int nStrands; // Number of strands in ordered complex
 int *pairs; // pairs[i] is the index of the base i is paired to
 int *seqlen; // seqlen[i] is the length of sequence i
 int *perm; // An array representing the permutation for the ordered complex
 char **seqs; // An array of the input sequences
 char line[MAXLINE]; // A line in the input file
 char linecpy[MAXLINE]; // A copy of a line in the input file
 char *tok; // token
 char tokseps[] = " ,;\t\n"; // Token separators 
 FILE *fp;
 
 // Open the input file
 if ((fp = fopen(InputFile,"r")) == NULL) {
   printf("ReadInputFile1: Error opening %s!\n\nExiting....\n",InputFile);
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
 
 // The first line is the number of sequences, nSeqs
 // Make sure it's number of strands
 if (!isdigit(line[0])) {
   printf("Error in input file [num strands] %s!\n",InputFile);
   printf("\nExiting...\n\n");
   exit(ERR_INPUT);
 }
 nSeqs = atoi(line);
 
 // Allocate sequences array
 seqs = (char **) malloc(nSeqs * sizeof(char *));
 seqlen = (int *) malloc(nSeqs * sizeof(int));
 
 // The next nSeqs lines are the sequences
 for (i = 0; i < nSeqs; i++) {
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
   
   seqlen[i] = strlen(line);  // sequence length
   seqs[i] = (char *) malloc((seqlen[i]+1) * sizeof(char));
   
   strcpy(seqs[i],line);
   
   // Chop off null character if it's there
   if (seqs[i][seqlen[i]-1] == '\n') {
     seqs[i][seqlen[i]-1] = '\0';
     (seqlen[i])--;
   }
 }
 
 // Go through sequences and check validity
 for (i = 0; i < nSeqs; i++) {
   for (j = 0; j < seqlen[i]; j++) {
     seqs[i][j] = toupper(seqs[i][j]);
     if (seqs[i][j] != 'A' && seqs[i][j] != 'C' && seqs[i][j] != 'G'
         && seqs[i][j] != 'T' && seqs[i][j] != 'U') {
           printf("Error in input sequences!  %c is not a valid base.\n",seqs[i][j]);
           printf("\nExiting....\n\n");
           exit(ERR_INPUT);
         }
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
 nStrands = 1;
 while ((tok = strtok(NULL,tokseps)) != NULL) {
   nStrands++;
 }
 
 // Allocate memory for the ordered complex
 perm = (int *) malloc(nStrands * sizeof(int));
 
 // Now put the entries into the perm
 tok = strtok(linecpy,tokseps);
 i = 0;
 perm[i++] = atoi(tok) - 1;
 while ((tok = strtok(NULL,tokseps)) != NULL) {
   perm[i++] = atoi(tok) - 1;
 }
 
 // Get total number of bases
 (*nbases) = 0;
 for (i = 0; i < nStrands; i++) {
   (*nbases) += seqlen[perm[i]];
 }
 
 // Allocate memory for pairs
 pairs = (int *) malloc((*nbases) * sizeof(int));
 
 // Now read in the base pairs
 if (fgets(line,MAXLINE,fp) == NULL) { // The means there are no base pairs
   for (i = 0; i < (*nbases); i++) {
     pairs[i] = -1;
   }
 }
 else if (isdigit(line[0]) == 0) { // Dot-paren format
   getStructureFromParens(line,pairs,(*nbases));
 }
 else { // Pair list format
   // Initialize pairs
   for (i = 0; i < *nbases; i++) {
     pairs[i] = -1;
   }
   tok = strtok(line,tokseps);
   i = atoi(tok)-1;
   tok = strtok(NULL,tokseps);
   j = atoi(tok)-1;
   pairs[i] = j;
   pairs[j] = i;
   while (fgets(line,MAXLINE,fp) != NULL) {
     if (line[0] != '%' && line[0] != '\n' && line[0] != '\0') {
       tok = strtok(line,tokseps);
       i = atoi(tok)-1;
       tok = strtok(NULL,tokseps);
       j = atoi(tok)-1;
       pairs[i] = j;
       pairs[j] = i;
     }
   }
 }
 
 fclose(fp);
 
 
 // Build the bases struct
 (*bases) = (struct base *) malloc((*nbases) * sizeof(struct base));
 i = 0; // This is the index of the base
 for (j = 0; j < nStrands; j++) {
   (*bases)[i].pair = pairs[i];
   (*bases)[i].prevConnection = 0;
   (*bases)[i].nextConnection = 1;
   (*bases)[i].strandID = perm[j];
   (*bases)[i].x[0] = 0.0; // Just to initialize
   (*bases)[i].x[1] = 0.0; // Just to initialize
   (*bases)[i].baseType = seqs[perm[j]][0];
   i++;
   for (k = 1; k < seqlen[perm[j]]-1; k++) {
     (*bases)[i].pair = pairs[i];
     (*bases)[i].prevConnection = 1;
     (*bases)[i].nextConnection = 1;
     (*bases)[i].strandID = perm[j];
     (*bases)[i].x[0] = 0.0; // Just to initialize
     (*bases)[i].x[1] = 0.0; // Just to initialize
     (*bases)[i].baseType = seqs[perm[j]][k];
     i++;
   }
   if (seqlen[perm[j]] == 1) { // A strand of 1 base is silly, but just in case
     (*bases)[i-1].nextConnection = 0;
   }
   else {
     (*bases)[i].pair = pairs[i];
     (*bases)[i].prevConnection = 1;
     (*bases)[i].nextConnection = 0;
     (*bases)[i].strandID = perm[j];
     (*bases)[i].x[0] = 0.0; // Just to initialize
     (*bases)[i].x[1] = 0.0; // Just to initialize
     (*bases)[i].baseType = seqs[perm[j]][seqlen[perm[j]]-1];
     i++;
   }
 }
 
 // Get data from base prob file
 if (*prob) {
   // Open the input file
   if ((fp = fopen(probFile,"r")) == NULL) {
     printf("ReadInputFile2: Error opening %s!\n\nExiting....\n",probFile);
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
   
   // The first line is prob for first base.  Read through all
   for (i = 0; i < *nbases && (*prob) == 1; i++) {
     (*bases)[i].p = str2double(line);
     if ((*bases)[i].p < 0.0 || (*bases)[i].p > 1.0) {
       printf("Error in prob file.  Probability not between 0 and 1.\n");
       printf("Ignoring probs.\n");
       *prob = 0;
     }
     if (i < *nbases - 1 && fgets(line,MAXLINE,fp) == NULL) {
       printf("Error in prob file.  Not enough entries.\n");
       printf("Ignoring probs.\n");
       *prob = 0;
     }
   }
   fclose(fp);
 }
 
 // Free memory
 free(perm);
 free(pairs);
 free(seqlen);
 for (i = 0; i < nSeqs; i++) {
   free(seqs[i]);
 }
 free(seqs);
 
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


