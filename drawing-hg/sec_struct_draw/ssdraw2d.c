/*
  SSDRAW2D.C

  This is the main routine for secondary structure drawing in 2D.  Generates an SVG
  illustration of a secondary structure disctated by the input file.  The drawing
  is colored in various ways, illustrating sequence and pair probability.

  ssdraw2d [-allowoverlap] [-docircle] [-svgfile filename] [-white] 
           [-prob filename] [-coords filename] [-material rna|dna] [-opacity]
           [-energy free_energy] [-energypreambe preable_text] [-units units_text]
           [-key] [-thymine] [-keybuffer] [-colorbarbuffer] [-freeenergybuffer]
           [-squarebuffer] prefix

  The input file is read from prefix.in.  It is formatted as follows:
  Line 1: The number of sequences, Nseq
  Line 2-Nseq+1: Each line contains a sequence
  Line Nseq+2: The ordering of the strands for the ordered complex, whitespace
               separated.  This is just a list of integers.  E.g. "1 2 3".
  Line Nseq+3: This may be a dot-paren structure.  Or this and subsequen lines
               are a list of pairs, (i j), whitespace separated.

  Option flags:
  -allowoverlap: Generate the canonical 2D structure, even if overlaps occur.
  -docircle: Do circle polymer graph.  This overrides allowoverlap.
  -svgfile: Allows specification of the name of the SVG file to which the graphic
            is written.  The required argument following the flag is the file name.
  -white: Swap black lines for white in the rendering.
  -prob: Enables shading of the bases based on how likely they are to be in the 
         depicted state at equilibrium.  The argument for the flag is the name of
         a file with the probabilities.  It is simply a list (one per line) of 
         the probability that a given base is in the state in the structure being
         predicted.  The number of lines must equal the number of bases and appear
         in the order the bases appear in the specifies ordering of the complex.
  -coords: Output a file containing the coordinates of the bases and do not write
           and SVG file.
  -material: Required argument.  Either rna or dna.  Uses actual geometry of these
             materials.  Default is to use standard non-physical drawing.
  -opacity: Probability is indicated with opacity as opposed to colormap.
  -stacknicks: Nicked helices are stacked.  CURRENTLY UNAVAILABLE.
  -energy: Argument of flag is printed on plot as free energy in %.2f format.
  -energypreamble: Argument of flag is printed before the free energy, followed 
                   immediately by a single whitespace and then the free energy.
                   Example: -energypreamble "Free energy:"
  -units: A string giving units of free energy, default is kcal/mol.
  -key: Selected if key for base coloring is shown on plot if -prob flag is inactive.
        If -prob flag is active, it shows a colorbar showing the base pair probability
        key.  If opacity is also selected, this is ignored.
  -thymine: Use T is used instead of U, only active with -key flag
  -keybuffer: Makes whitespace on the drawing of the same size the key would be on
              the right side of the drawing.
  -colorbarbuffer: Makes whitespace on the drawing of the same size the colorbar would
                   be on the right side of the drawing.
  -freeenergybuffer: Makes whitespace on the drawing of the same size the free 
                     energy text would be on the bottom of the drawing
  -squarebuffer: Makes white space around the entire figure such that the overall
                 dimension is still square.  The size of the boundary is given by
                 whichever would be biggest: the space taken by the colorbar, base 
                 key, or free energy text.  This flag supercedes all other buffering
                 flags.
  -drawbases Puts base letters on circles
  Justin Bois, September 2007
*/


#include "SecStructDrawHeader.h" // Header file for Concentrations

/* ******************************************************************************** */
/*                                BEGIN MAIN                                        */
/* ******************************************************************************** */

int main(int argc, char *argv[]) {

  int i; // Counter
  int nbases; // Total number of bases in structure
  int drawCircle = 0; // = 1 is ladder-hose drawing failed and we render circle
  int drawBases = 0; // Draws base letters on circle
  int ColorStrands; // = 1 for coloring of the backbone
  int AllowOverlap; // = 1 if overlap is allowed
  int doCircle; // = 1 if only circle drawing is to be done
  int white; // == 1 if lines are white
  int prob; // == 1 if shade bases with probability
  int opacity; // == 1 if probability is given by opacity
  int key; // == 1 if the color key for bases is to be included
  int thymine; //  == 1 if use T instead of U in base color key
  int stackNicks; // == 1 if nicked helices are stacked
  int showFreeEnergy; // == 1 if free energy is displayed on SVG graphic
  int coords; // == 1 if write out coordinates of bases  (and not SVG)
  int material; // == 1 for RNA, == 2 for DNA
  int keyBuffer; // == 1 if white space is added if there's no key
  int colorbarBuffer; // == 1 if white space is added if there's no colorbar
  int freeEnergyBuffer; // == 1 if white space is added if there's no colorbar
  int squareBuffer; // == 1 if white space is added to make square drawings with
                    // The drawing in the center and enough space for colorbar or
                    // base key and free energy.  The trim is such that the image 
                    // is square
  int frameBuffer; // == 1 if white space is added to make square drawings with
                    // The drawing in the center and enough space for colorbar or
                    // base key and free energy.
  int canonicalSuccessful; // = 1 if canonical drawing was successful, 0 otherwise
  int isPK; // = 1 if it is a pseudoknot
  double **StandardArrow; // The coordinates for a standard arrow aligned on x-axis
  double ArrowHeight; // The total height of an arrowhead
  double ArrowRadius; // Half the width of the arrowheads
  double freeEnergy; // Value of free energy to be shown on SVG graphic
  char units[MAXLINE]=""; // Units for free energy
  char freeEnergyPreamble[MAXLINE]=""; // Preable for free energy
  char InputFile[MAXLINE]=""; // The name of the file from which the input is read
  char OutputFile[MAXLINE]=""; // The file to which the SVG is written
  char probFile[MAXLINE]=""; // The file to which the SVG is written
  char coordFile[MAXLINE]=""; // The file to which the coordinates are written
  struct base *bases; // Struct of information about the bases.
  

  // Get input data
  ReadCommandLine(argc,argv, &drawBases, &ColorStrands,&AllowOverlap,&doCircle,&white,&prob,
    &opacity,&coords,&stackNicks,&material,&showFreeEnergy, &freeEnergy,
    units,freeEnergyPreamble,&key,&thymine,&keyBuffer,&colorbarBuffer,
    &freeEnergyBuffer,&squareBuffer,&frameBuffer,InputFile,OutputFile,
    probFile,coordFile);
  ReadInputFile(&bases,&nbases,InputFile,&prob,probFile);

  isPK = isPseudoknot(bases,nbases);
 
  if (isPK && coords) {
    printf("Error: Cannot generate coord file for pseudoknotted structure.\n");
    printf("Exiting...\n\n");
    exit(1);
  }

  // Get coordinates for standard arrow
  getStandardArrow(&StandardArrow,&ArrowHeight,&ArrowRadius);
  
  // Determine if loop is pseudoknotted
  if (doCircle || isPK) {
    drawCircle = 1;
  }
  else {
    canonicalSuccessful = DrawCanonicalSS(bases,nbases,StandardArrow,ArrowHeight,
                                          ArrowRadius,ColorStrands,AllowOverlap,
                                          white,prob,opacity,drawBases, showFreeEnergy,
                                          freeEnergy,units,freeEnergyPreamble,key,
                                          thymine,coords,material,keyBuffer,
                                          colorbarBuffer,freeEnergyBuffer,
                                          squareBuffer,frameBuffer,OutputFile,
                                          coordFile);
    
    if (!canonicalSuccessful) {
      drawCircle = 1;
    }
  }

  if (drawCircle) {
    RenderCirclePolyGraph(bases,nbases,StandardArrow,ArrowHeight,ArrowRadius,
                  ColorStrands,white,prob,opacity,showFreeEnergy,freeEnergy,
                  units,freeEnergyPreamble,key,thymine,keyBuffer,
                  colorbarBuffer,freeEnergyBuffer,squareBuffer,frameBuffer, 
                  OutputFile);
  }

  for (i = 0; i < 6; i++) {
    free(StandardArrow[i]);
  }
  free(StandardArrow);
  free(bases);

  return 0;

}
/* ******************************************************************************** */
/*                                  END MAIN                                        */
/* ******************************************************************************** */


