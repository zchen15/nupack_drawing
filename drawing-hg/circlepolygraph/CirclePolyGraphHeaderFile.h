#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <getopt.h>

// Constants
#define MAXLINE 1000 // Maximum characters in a line
#define PI 3.1415926535897932384626
#define SVG_SIZE 500.0 // SVG height in pixels
#define INTERBASE_DISTANCE 50.0 // Distance between bases
#define BACKBONE_WIDTH_FRAC 0.1 // Backbone line width as fraction of interbase dist
#define TICK_WIDTH_FRAC 0.033 // Tick line width as fraction of interbase dist
#define LINE_WIDTH_FRAC 0.05 // bp line width as fraction of interbase dist
#define TICK_LENGTH_FRAC 0.2 // tick length as a fraction of interbase dist
#define GAP_FACT 1.5 // gap = INTERBASE_DISTANCE * GAP_FACT
#define OVERHANG_FRAC 1.0 // Fraction of interbase distance left at ends of strands
#define ARROW_HEIGHT_FRAC 0.33 // arrow height as fraction of interbase dist.
#define ARROW_ANGLE PI/4.0 // Angle arrow wings make with vertical
#define ARROW_FRAC 0.5 // Fraction of total arrow height under wings
#define TOLERANCE 1e-6 // Tolerance for Newton's method to converge
#define MAX_ITERATIONS 1000 // Maximal number of Newton iterations

// Colors for strands
static int StrandColors[10][3] = {{0,0,128}, // Navy blue
				  {164,0,0}, // Red
				  {34,139,34}, // Forrest green
				  {79,0,147}, // Purple
				  {255,73,0}, // Orange
				  {135,206,235}, // Light blue
				  {255,169,204}, // Light red
				  {128,255,89}, // Light green
				  {181,145,209}, // Lavender
				  {0,0,0}}; // black

#define MAXSTRANDS 10 // Maximal number of strand sequences allowed


// Where the help file to print is
#define CIRCLEPOLYGRAPH_HELP_FILE "CirclePolyGraph.help"

// Error codes (start with 30 so as not to overlap with Complexes error code)
#define ERR_INPUT 1 // Error in input file
#define ERR_INPUTFILE 2 // Error in opening the input file
#define ERR_NOINPUT 3 // User supplied no input
#define ERR_INVALIDSTRUCTURE 4 // User supplied invalid dot-paren structure
#define ERR_NEWTON_FAIL 5 // Newton's method failed to converge
#define ERR_OUTPUTFILE 6 // Error opening output file
#define ERR_INVALID_ID 7 // Invalid comp or perm id in command line input

// Struct that carries base information
typedef struct base {
  int pair;  // index of this base's pair, -1 if none
  int next;  // index of next base
  int prev;  // index of previous base
  int nextConnection; // = 1 if base is connected to next one
  int prevConnection; // = 1 if base is connected to previous one
  double x; // Position of base
  double y; // Position of base
  char baseType; // 'A', 'T', 'U', 'G', or 'C'
} base;

// Struct that has the arc informtion
typedef struct arc {
  double x0[2];
  double x1[2];
  double radius;
  int LargeArc; // = 1 if arc has angle > PI
  int Sweep; // == 1 if arc is goes clockwise start to end
  int RGB[3]; // Color
} arc;

// Struct that has the arrow informtion
typedef struct arrow {
  // Positions
  double x0[2];
  double x1[2];
  double x2[2];
  double x3[2];
  double x4[2];
  double x5[2];
  int RGB[3]; // Color
} arrow;

// Struct for lines between base pairs
typedef struct line {
  // Start and end point
  double x0[2];
  double x1[2];
  int RGB[3]; // Color
} line;

// Struct that has the tick information
typedef struct tick {
  // Positions of start and end
  double x0[2];
  double x1[2];
  int RGB[3]; // Color
} tick;


// Function prototypes
/* ************************ IN READCOMMANDLINE.C ********************************** */
void ReadCommandLine(int nargs, char **args, int *ShowSequences, int *ColorStrands,
		     int *LabelStrands, char *InputFile, char *OutputFile);
/* ******************************************************************************** */

/* ************************* IN UTILS.C ******************************************* */
double str2double (char *str);
void MatrixVectorMult(double *c, double **A, double *b, int n);
double dot(double *v1, double *v2, int len);
/* ******************************************************************************** */

/* ************************** IN INPUTFILEREADER.C ******************************** */
void ReadInputFile(int *K, int *nStrands, int *nBases, int **SeqLen, int **Perm,
		   int **pairs, int *nPairs, char *InputFile);
void getStructureFromParens( char *line, int *pairs, int seqlength);
/* ******************************************************************************** */

/* ************************* IN CIRCLEPOLYGRAPH.C ********************************* */
double CircleRadius(int nBases, int nStrands, double b0, double gap, double fb0,
		    double a, double tol, int MaxIters);
void SVGHeader(char *OutputFile, double SVGWidth, double SVGHeight, double vbw,
	       double vbh);
void SVGTicks(char *OutputFile, struct tick *ticks, int nBases, double TickWidth);
void SVGArcs(char *OutputFile, struct arc *arcs, struct arrow *arrows, int nStrands, 
	     double ArcWidth);
void SVGLines(char *OutputFile, struct line *lines, int nPairs, double LineWidth);
void SVGEnd(char *OutputFile);
void getArrow(struct arrow *arrows, int i, double **StandardArrow, double theta, 
	      double r);
/* ******************************************************************************** */


