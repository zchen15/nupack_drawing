/*
  READCOMMANDLINE.C

  Justin Bois, Caltech, 2 September 2006
  bois@caltech.edu
*/

#include "CirclePolyGraphHeaderFile.h" // Concentrations header file

/* ******************************************************************************** */
void ReadCommandLine(int nargs, char **args, int *ShowSequences, int *ColorStrands,
		     int *LabelStrands, char *InputFile, char *OutputFile) {

  int options;  // Counters used in getting flags
  int ShowHelp; // ShowHelp = 1 if help option flag is selected
  char UserInput[MAXLINE]; // The input the user gives for the input file
  FILE *fp; // The cx file, used to check if we can open it.

  if (nargs == 1) {
    printf("For instructions on running this program, run it with the ");
    printf("-help flag.\n\nExiting....\n\n");
    exit(ERR_NOINPUT);
  }


  // Initialize input flags to defaults
  *ShowSequences = 0;
  *ColorStrands = 1;
  *LabelStrands = 1;
 
  // Get the option flags
  while (1)
    {
      static struct option long_options [] =
	{
          {"sequences", no_argument,         0, 'a'},
	  {"nocolor", no_argument,           0, 'b'},
	  {"nolabel", no_argument,           0, 'c'},
	  {"help", no_argument,              0, 'd'},
          {0, 0, 0, 0}
        };
      /* getopt_long stores the option index here. */
      int option_index = 0;

      options = getopt_long_only (nargs, args, 
				  "abcd", long_options, 
				  &option_index);

      // Detect the end of the options.
      if (options == -1)
        break;

      switch (options)
        {
        case 'a':
	  *ShowSequences = 1;
          break;

	case 'b':
	  *ColorStrands = 0;
	  break;

	case 'c':
	  *LabelStrands = 0;
	  break;
	  
	case 'd':
	  ShowHelp = 1;
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
    printf("For instructions on running this program, run it with the ");
    printf("-help flag.\n\nExiting....\n\n");
    exit(ERR_NOINPUT);
  }
  else {
    strcpy(UserInput,args[optind]);
  }


  // Name the files
  strcpy(InputFile,UserInput);
  strcat(InputFile,".in");
  strcpy(OutputFile,UserInput);
  strcat(OutputFile,"-mfe_polygraph.svg");


  // Do a quick check to make sure the cx file exists before we proceed
  if ((fp = fopen(InputFile,"r")) == NULL) {
    printf("Error opening %s!\n\nExiting....\n",InputFile);
    exit(ERR_INPUTFILE);
  }
  fclose(fp);

} 
/* ******************************************************************************** */

