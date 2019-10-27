/*
  READCOMMANDLINE.C

  Justin Bois, Caltech, 2 September 2006
  bois@caltech.edu
*/

#include "SecStructDrawHeader.h" // Secondary structure drawing header file

/* ******************************************************************************** */
void ReadCommandLine(int nargs, char **args, int *drawBases, int *ColorStrands, int *AllowOverlap, 
                     int *doCircle, int *white, int *prob, int *opacity, int *coords, 
                     int *stackNicks, int *material, int *showFreeEnergy,
                     double *freeEnergy, char *units, char *freeEnergyPreamble,
                     int *key, int *thymine, int *keyBuffer, int *colorbarBuffer,
                     int *freeEnergyBuffer, int *squareBuffer, int *frameBuffer,
                     char *InputFile, char *OutputFile, char *probFile, 
                     char *coordFile) {
 
 int options;  // Counters used in getting flags
 int UserOutput = 0; // = 1 if user gives output file name
 char UserInput[MAXLINE]; // The input the user gives for the input file
 FILE *fp; // The cx file, used to check if we can open it.
 
 if (nargs == 1) {
   //    printf("For instructions on running this program, run it with the ");
   //    printf("-help flag.\n\nExiting....\n\n");
   printf("Error!  No input!\n");
   exit(ERR_NOINPUT);
 }

 // Initialize input flags to defaults
 *ColorStrands = 0;
 *AllowOverlap = 0;
 *doCircle = 0;
 *white = 0;
 *prob = 0;
 *opacity = 0;
 strcpy(units,"kcal/mol");
 strcpy(freeEnergyPreamble,"");
 *key = 0; 
 *thymine = 0;
 *coords = 0;
 *showFreeEnergy = 0;
 *freeEnergy = 0.0;
 *material = 2; // Default is DNA
 *keyBuffer = 0;
 *colorbarBuffer = 0;
 *freeEnergyBuffer = 0;
 *squareBuffer = 0;
 *frameBuffer = 0;
 
 
 // Get the option flags
 while (1)
 {
   static struct option long_options [] =
   {
     {"colorstrands", no_argument,           0, 'a'},
     {"allowoverlap", no_argument,           0, 'b'},
     {"docircle", no_argument,               0, 'c'},
     {"svgfile", required_argument,          0, 'd'},
     {"white", no_argument,                  0, 'e'},
     {"prob", required_argument,             0, 'f'},
     {"coords", required_argument,           0, 'g'},
     {"material", required_argument,         0, 'h'},
     {"opacity", no_argument,                0, 'i'},
     {"stacknicks", no_argument,             0, 'j'},
     {"energy", required_argument,           0, 'k'},
     {"units", required_argument,            0, 'l'},
     {"key",   no_argument,                  0, 'm'},
     {"thymine", no_argument,                0, 'n'},
     {"energypreamble", required_argument,   0, 'o'},
     {"keybuffer", no_argument,              0, 'p'},
     {"colorbarbuffer", no_argument,         0, 'q'},
     {"freeenergybuffer", no_argument,       0, 'r'},
     {"squarebuffer", no_argument,           0, 's'},
     {"framebuffer",  no_argument,           0, 't'},
     {"drawbases",  no_argument,             0, 'u'},
     {0, 0, 0, 0}
   };
   /* getopt_long stores the option index here. */
   int option_index = 0;
   
   options = getopt_long_only (nargs, args, 
                               "abcd:ef:g:h:ijk:l:mno:pqrstu", long_options, 
                               &option_index);

   // Detect the end of the options.
   if (options == -1)
     break;
   printf("options=%c, optarg=%s\n",options,optarg);
   switch (options)
   {
     case 'a':
         *ColorStrands = 1;
         break;
         
       case 'b':
         *AllowOverlap = 1;
         break;
         
       case 'c':
         *doCircle = 1;
         break;
         
       case 'd':
         strcpy(OutputFile,optarg);
         UserOutput = 1;
         break;
         
       case 'e':
         *white = 1;
         break;
         
       case 'f':
         strcpy(probFile,optarg);
         *prob = 1;
         break;
         
       case 'g':
         strcpy(coordFile,optarg);
         *coords = 1;
         break;
         
       case 'h':
         strcpy(UserInput,optarg);
         if (strlen(UserInput) >= 3) {
           if (strncmp(UserInput,"rna",3) == 0) {
             *material = 1;
           }
           else if (strncmp(UserInput,"dna",3) != 0) {
             printf("Invalid material parameter.  Using B-DNA.\n");
           }
         }
         else {
           printf("Invalid material parameter.  Using B-DNA.\n");
         }
         break;
         
       case 'i':
         *opacity = 1;
         break;
         
       case 'j':
         *stackNicks = 1;
         break;
         
       case 'k':
         *showFreeEnergy = 1;
         strcpy(UserInput,optarg);
         *freeEnergy = str2double(UserInput);
         break;
         
       case 'l':
         strcpy(units,optarg);
         break;
         
       case 'm':
         *key = 1;
         break;
         
       case 'n':
         *thymine = 1;
         break;
         
       case 'o':
         strcpy(freeEnergyPreamble,optarg);
         break;
         
       case 'p':
         *keyBuffer = 1;
         break;
         
       case 'q':
         *colorbarBuffer = 1;
         break;
         
       case 'r':
         *freeEnergyBuffer = 1;
         break;
         
       case 's':
         *squareBuffer = 1;
         break;
         
       case 't':
         *frameBuffer = 1;
         break;

       case 'u':
         *drawBases = 1;
         break;
         
       case '?':
         // getopt_long already printed an error message.
         break;
         
       default:
         abort ();
   }
 }
 
 
 // Get the the input file
 if (optind == nargs) { // There's no input from the user
   printf("You must have a prefix or an input file on the command line.\n");
   exit(ERR_NOINPUT);
 }
 else {
   strcpy(UserInput,args[optind]);
 }
 
 printf("UserInput=args[%d]=%s\n",optind,args[optind]);
 // Name the files
 strcpy(InputFile,UserInput);
 strcat(InputFile,".in");
 if (!UserOutput) {
   strcpy(OutputFile,UserInput);
   strcat(OutputFile,"-mfe_drawing.svg");
 }
 
 // Check for consitency of flags
 if (*coords) {
   if (*prob || *white || *doCircle || *ColorStrands) {
     printf("-coord selected: no SVG will be written.\n");
     *prob = 0;
     *white = 0;
     *doCircle = 0;
     *ColorStrands = 0;
   }
   if (*AllowOverlap == 0) {
     printf("-coord selected: overlaps allowed.\n");
     *AllowOverlap = 1;
   }
 }
 
 if (*key) {
   if (*prob) { // ignore key
     *key = 0;
   }
   else { // Key is shown and no colorbar
     *keyBuffer = 1;
   }
 }
 
 if (*prob) {
   *colorbarBuffer = 1;
 }
 
 if (*showFreeEnergy) {
   *freeEnergyBuffer = 1;
 }
 
 
 if (*squareBuffer) { // Square buffer overrides all others
   if (*keyBuffer || *colorbarBuffer || *freeEnergyBuffer || *frameBuffer) {
     *keyBuffer = 0;
     *colorbarBuffer = 0;
     *freeEnergyBuffer = 0;
     *frameBuffer = 0;
   }
 }
 else if (*frameBuffer) {
   // Alert the user that right now square buffer is not possible and exit.
   printf("FRAME BUFFERING NOT YET AVAILABLE!  EXITING....\n\n");
   exit(1);
   if (*keyBuffer || *colorbarBuffer || *freeEnergyBuffer) {
     *keyBuffer = 0;
     *colorbarBuffer = 0;
     *freeEnergyBuffer = 0;
   }
 }
 
 printf("ReadCommandLine: Inputfile=%s\n",InputFile);
 // Do a quick check to make sure the input file exists before we proceed
 if ((fp = fopen(InputFile,"r")) == NULL) {
   printf("ReadCommandLine: Error opening %s!\n\nExiting....\n",InputFile);
   exit(ERR_INPUTFILE);
 }
 fclose(fp);
 
} 
/* ******************************************************************************** */

