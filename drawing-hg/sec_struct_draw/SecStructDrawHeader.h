#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <getopt.h>

// Definition constants for drawing
#define SVG_SIZE 500.0 // SVG height in pixels
#define INTERBASE_DISTANCE 50.0 // Distance between bases in ss region
#define STACK_DIST_FRAC 3.4/4.3 // Distance between stacked bases as a fraction of interbase dist.
#define BASE_PAIR_DIST_FRAC STACK_DIST_FRAC*1.618
              // interbase distance of paired bases as frac of above
//#define STACK_DIST_FRAC 1.0
//#define BASE_PAIR_DIST_FRAC 1.0
#define N_PHANTOM_BASES_PER_NICK 1 // Number of phantom bases per nick
#define N_PHANTOM_1 3 // Number of phantom bases for sticky end length 1
#define N_PHANTOM_2 2 // Number of phantom bases for sticky end length 2
#define NICK_PHANTOM_DIST_FRAC 1.0 // Interbase distance for phantom base at nick
#define EDGE_BUFFER_FRAC 0.001 // Buffer on image border as frac of total width
#define BASE_RADIUS_FRAC 0.175 // Radius of dot as frac of interbase dist.
#define BASE_STROKE_FRAC 0.5 // Width of solid ring in circle as fraction of total radius
#define BACKBONE_WIDTH_FRAC 0.1 // Backbone line width as fraction of interbase dist
#define BP_WIDTH_FRAC 0.05 // bp line width as fraction of interbase dist
#define CLASH_BUFFER_FRAC BASE_RADIUS_FRAC*0.1 // buffer for clash as a fraction of interbase dist
#define GAP_FACT 1.5 // gap = INTERBASE_DISTANCE * GAP_FACT
#define OVERHANG_FRAC 1.0 // Fraction of interbase distance left at ends of strands for cirle plot
#define ARROW_LINE_FRAC 0.1 // Arrow line with as fraction of interbase distance
#define ARROW_OVERHANG_FRAC 0.55 // Fraction of interbase distance at end of strand in ladder-hose
// Value of ARROW_OVERHANG_FRAC of 0.5 works if no stacking of helices
#define ARROW_HEIGHT_FRAC 0.33 // arrow height as fraction of interbase dist.
#define ARROW_ANGLE PI/4.0 // Angle arrow wings make with vertical
#define ARROW_FRAC 0.5 // Fraction of total arrow height under wings
#define COLORBAR_HEIGHT_TO_WIDTH 67.0 // Height to width ratio of colorbar (roughly same as dot plot)
#define COLORBAR_HEIGHT_FRAC 0.8 // Fraction of total SVG height taken by color bar
#define N_COLORBAR_COLORS 256 // Number of color increments in the colorbar
#define N_COLORBAR_LABELS 6 // Number of labels on the colorbar
#define COLORBAR_LABEL_FONT_SIZE_FRAC 0.035 // label font size as fraction of total colorbar height
#define COLORBAR_LABEL_OFFSET_FRAC 0.5 // Offset from colorbar edge of labels as fraction of colorbar width
#define COLOR_BAR_VERT_FUDGE_FACT_FRAC 1.0/3.0 // Fraction of colorbar label font size to shift color bar down
#define COLORBAR_OFFSET_FRAC 0.05 // Fraction of total image width the colorbar is offset
#define COLORBAR_OVERSHOOT_FRAC 2.0 // Amount to buffer right margin as fraction of colorbar label font size
#define BASE_KEY_FONT_SIZE_FRAC_TOTAL 0.035 // Base key font size as fraction of total image size
#define BASE_KEY_DOT_RADIUS_FRAC 0.0125 // Dot radius of base key as fraction of height of view box
#define BASE_KEY_SEP_FRAC 5.0 // Separation between base key entries as fraction of base key dot radius
#define BASE_KEY_FONT_SIZE_FRAC 3.0 // Base key font size as fraction of base key dot radius
#define BASE_KEY_LABEL_OFFSET_FRAC 2.0 // Base key label offset from dot as fraction of base key dot radius
#define BASE_KEY_VERT_FUDGE_FACT_FRAC 1.0 // Fraction of base key font size to shift base key label down
#define BASE_KEY_OFFSET_FRAC 0.05 // Fraction of total image width the base key is offset
#define BASE_KEY_OVERSHOOT_FRAC 0.0 // Amount to buffer right margin as fraction of base key label font size
#define FREE_ENERGY_FONT_SIZE_FRAC 0.035 // Free energy font size as fraction of total image height
#define FREE_ENERGY_HEIGHT_FRAC 1.1 // Size of space for free energy as fraction of text size. Should be > 1. 
#define FREE_ENERGY_SIDE_BUFFER_FRAC 0.025 // Buffer on left of free energy
#define FREE_ENERGY_TOP_BUFFER_FRAC 0.025 // Space from bottom of drawing to top of free energy space as fraction of image height
#define FREE_ENERGY_BOTTOM_BUFFER_FRAC 0.025 // Buffer on bottom of free energy
#define FONT "sans-serif" // Font family

// Constants for outputting coordinates (based on values in Niles's NUDRAW)
// DNA:
#define INTERBASE_DISTANCE_DNA 6.11812 // = Niles's dsb
#define STACK_DIST_FRAC_DNA 0.55573 // = Niles's dzb/dsb != 3.4/4.3. which is commonly published value
#define BASE_PAIR_DIST_FRAC_DNA 2.77863 // = Niles's 2*rdh/dsb
// RNA:
#define INTERBASE_DISTANCE_RNA 7.06463 // = Niles's dsb
#define STACK_DIST_FRAC_RNA 0.36803 // = Niles's dzb/dsb
#define BASE_PAIR_DIST_FRAC_RNA 3.25566 // = Niles's 2*rdh/dsb


// Shading parameters
// The opacity of the bases is piece-wise linear in the probability that the base
// is in the given configuration at equilibrium.  The equation for the piece-wise
// linear function is: opacity = m1*p for 0 <= p <= P_SLOPE_CHANGE
//                             = m2*p + b2 for P_SLOPE_CHANGE <= p <= 1
//  where m1 = OPACITY_SLOPE_CHANGE / P_SLOPE_CHANGE,
//        m2 = OPACITY_SLOPE_CHANGE / (1 - P_SLOPE_CHANGE),
//  and   b2 = 1 - P_SLOPE_CHANGE / (1 - P_SLOPE_CHANGE)
#define P_SLOPE_CHANGE 0.9
#define OPACITY_SLOPE_CHANGE 0.5

// Constants
#define MAXLINE 100000 // Maximum characters in a line
#define PI 3.1415926535897932384626
#define TOLERANCE 1e-6 // Tolerance for Newton's method to converge
#define MAX_ITERATIONS 1000 // Maximal number of Newton iterations
#define MAXSTRANDS 10 // Maximal number of strand sequences allowed

// Error codes (start with 30 so as not to overlap with Complexes error code)
#define ERR_INPUT 31 // Error in input file
#define ERR_INPUTFILE 32 // Error in opening the input file
#define ERR_NOINPUT 33 // User supplied no input
#define ERR_INVALIDSTRUCTURE 34 // User supplied invalid dot-paren structure
#define ERR_NEWTON_FAIL 35 // Newton's method failed to converge
#define ERR_OUTPUTFILE 36 // Error opening output file
#define ERR_DISCONNECTED 37 // input structure is disconnected

// Colors for strands
static int DefaultStrandColor[3] = {0,0,0}; // Black
static int StrandColorWhite[3] = {255,255,255}; // White

static int BasePairColor[3] = {0,0,0}; // Black
static int BasePairColorWhite[3] = {255,255,255}; // White

// Not used at the moment
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

// Colors for bases (order: A (green), C (blue), G  (black), T/U (red))
static int BaseColors[4][3] = {{34,139,34}, // Forrest green
                               {24,73,255}, // Blue
                               {0,0,0}, // Black
                               {164,0,0}}; // Red

// Base colors when -white flag is active
static int BaseColorsWhite[4][3] = {{34,139,34}, // Forrest green
                                     {24,73,255}, // Blue
                                     {255,255,255}, // White
                                     {164,0,0}}; // Red



// Struct that carries base information
typedef struct base {
  int pair;  // index of this base's pair, -1 if none
  int nextConnection; // = 1 if base is connected to next one
  int prevConnection; // = 1 if base is connected to previous one
  int strandID; // The indentifier number of the strand base belongs to
  int loop; // Index of the loop to which it belongs
  int helix; // Index of helix to which it belongs
  double x[2]; // The x-y position of the base
  double p; // The probability that the base is in state depicted in structure (0 to 1)
  char baseType; // 'A', 'T', 'U', 'G', or 'C'
} base;


// Struct for loop information
typedef struct loop {
  int nbases;  // Number of bases that are members of the loop
  int nnicks;  // Number of nicks in the loop (either 1 or 0 for connected)
  int nhelices;  // Number of helices in loop
  int nphantom; // Number of phantom bases in the loop
  int *bases;  // List of base indicies for loop
  double radius; // Loop radius
  int bigAngle; // = 1 if loop lies entirely above circle enscribing it.
  double center[2];  // Coordinates for loop center
  double **x; // Position of the bases
} loop;


// Struct for helix information
typedef struct helix {
  int nbases; // Number of bases in helix
  int *bases; // List of bases in helix
  int corners[4]; // Helix defined by [i,j,k,l] with i,j paired and k,l paired
                  // with i < k < l < j
  int loops[2]; // loop[0] is 5' loop and loop[1] is 3'loop
  double **x; // Position of the bases
} helix;


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

// Struct for line between base pairs
typedef struct arrowline {
  // Start and end point
  double x0[2];
  double x1[2];
  int RGB[3]; // Color
} arrowline;

// Struct for line between base pairs
typedef struct basepair {
  // Start and end point
  double x0[2];
  double x1[2];
  int RGB[3]; // Color
} basepair;

// Struct for lines connecing bases in backbone
typedef struct backbone {
  // Start and end point
  double x0[2];
  double x1[2];
  double radius; // Radius of loop to which it belongs
  int curved; // = 1 if line is curved
  int RGB[3]; // Color
} backbone;

// Struct for dots
typedef struct basedot {
  double x[2];
  int RGB[3]; // Color
  double p; // Level of opacity
} basedot;


// Function prototypes
/* ************************ IN READCOMMANDLINE.C ********************************** */
void ReadCommandLine(int nargs, char **args, int *drawBases, int *ColorStrands, int *AllowOverlap, 
                     int *doCircle, int *white, int *prob, int *opacity, int *coords, 
                     int *stackNicks, int *material, int *showFreeEnergy,
                     double *freeEnergy, char *units, char *freeEnergyPreamble,
                     int *key, int *thymine, int *keyBuffer, int *colorbarBuffer,
                     int *freeEnergyBuffer, int *squareBuffer, int *frameBuffer,
                     char *InputFile, char *OutputFile, char *probFile, 
                     char *coordFile);
/* ******************************************************************************** */

/* ************************* IN UTILS.C ******************************************* */
int isPseudoknot(struct base *bases, int nbases);
void BaseColor(int *RGB, char baseType,int white);
void ColorMap(int *RGB, double prob);
void getStandardArrow(double ***StandardArrow, double *ArrowHeight, 
		      double *ArrowRadius);
void getArrow(struct arrow *arrows, int i, double **StandardArrow, double pos[2], 
	      double dir[2]);
double CircleRadius(int *bigAngle,double *lengths, int nsegments, double tol, 
		    int MaxIters);
void SVGheader(char *OutputFile, double SVGWidth, double SVGHeight, double leftCorner,
	       double topCorner, double vbw, double vbh);
void SVGarcs(char *OutputFile, struct arc *arcs, struct arrow *arrows, int nStrands, 
	     double ArcWidth);
void SVGarrows(char *OutputFile, struct arrow *arrows, struct arrowline *arrowlines, 
	       int nStrands, double arrowLineWidth);
void SVGbasepairs(char *OutputFile, struct basepair *basepairs, int nPairs,
		  double BPWidth);
void SVGbackbone(char *OutputFile, struct backbone *backbones, int nbackbone, 
		 double BackboneWidth);
void SVGdots(char *OutputFile, struct basedot *dots, struct base *bases, int nbases, double DotRadius, 
	     int prob, int opacity, int drawbases);
void SVGColorBar(char *OutputFile, double colorBarx, double colorBarShifty, 
		 double dim, int white);
void SVGBaseKey(char *OutputFile, double dotRadius,  double baseKeyx, double dim, 
		int thymine, int white);
void SVGEnergy(char *OutputFile, double freeEnergy, char *units, 
	       char *freeEnergyPreamble, double xpos, double ypos, double dim,
	       int white);
void SVGend(char *OutputFile);
void getRfromTheta(double **R, double theta);
void MatrixVectorMult(double *c, double **A, double *b, int n);
double dot(double *v1, double *v2, int len);
double max(double *ar, int len);
double norm(double *ar, int len);
int doubleGreaterCmp(const void *p1, const void *p2);
double min2(double a, double b);
double max2(double a, double b);
double str2double (char *str);
double pointpointdist(double *x1, double *x2, int n);
double pointlinedist(double *x, double *line, int n);
int linecross2d(double line1[4], double line2[4]);
void printBases(struct base *bases,int nbases);
void printLoops(struct loop *loops,int nloops);
void printHelices(struct helix *helices,int nhelices);
/* ******************************************************************************** */

/* ************************** IN INPUTFILEREADER.C ******************************** */
void ReadInputFile(struct base **bases, int *nbases, char *InputFile, int *prob, 
                   char *probFile);
void getStructureFromParens( char *line, int *pairs, int seqlength);
/* ******************************************************************************** */

/* ********************** IN RENDERCIRCLEPOLYGRAPH.C ****************************** */
void RenderCirclePolyGraph(struct base *bases, int nbases, double **StandardArrow,
     double ArrowHeight, double ArrowRadius, int ColorStrands,
     int white, int prob, int opacity, int showFreeEnergy, 
     double freeEnergy, char *units, char *freeEnergyPreamble, 
     int key, int thymine, int keyBuffer, int colorbarBuffer, 
     int freeEnergyBuffer, int squareBuffer, int frameBuffer, 
     char *OutputFile);
/* ******************************************************************************** */

/* ************************* IN DRAWCANONICALSS.C ********************************* */
int DrawCanonicalSS(struct base *bases, int nbases, double **StandardArrow,
    double ArrowHeight, double ArrowRadius, int ColorStrands,
    int AllowOverlap, int white, int prob, int opacity, int drawBases,
    int showFreeEnergy, double freeEnergy, char *units, 
    char *freeEnergyPreamble, int key, int thymine, int coords, 
    int material, int keyBuffer, int colorbarBuffer, 
    int freeEnergyBuffer, int squareBuffer, int frameBuffer,
    char *OutputFile, char *coordFile);

void getLoopsHelices(struct loop **loops, struct helix **helices, int *nloops, 
     int *nhelices, struct base *bases, int nbases);

int detectClashes(struct base *bases, int nbases);
/* ******************************************************************************** */


