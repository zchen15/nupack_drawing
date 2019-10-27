/*
*/


#include "CirclePolyGraphHeaderFile.h" // Header file for Concentrations

/* ******************************************************************************** */
/*                                BEGIN MAIN                                        */
/* ******************************************************************************** */

int main(int argc, char *argv[]) {
  
  int i,j,k; // Counters
  int K; // number of strands
  int ShowSequences; // = 1 if sequence is on circle diagram
  int ColorStrands; // = 1 if strands are colored
  int LabelStrands; // = 1 if strands are to be labeled
  int nStrands; // Number of strands in circle diagram
  int nBases; // Total number of bases in circle diagram
  int nStrandTypes; // Number of different sequences on circle diagram
  int ColorKey[MAXSTRANDS]; // ColorKey[i] is the entry in StrandColors for strand i
  int *SeqLen; // Sequence lengths
  int *Perm; // The permutation
  int  *pairs; // List of base pairs (pairs[i] is the base i is paired to)
  int nPairs; // Number of pairs
  double SVGWidth, SVGHeight; // Dimensions of SVG graphic
  double vbw, vbh; // View box width and height
  double r; // The radius of the circle
  double b0 = INTERBASE_DISTANCE; // inter-base length
  double gap = GAP_FACT * INTERBASE_DISTANCE; // Gap size
  double a = INTERBASE_DISTANCE*ARROW_HEIGHT_FRAC; // The arrow height
  double **StandardArrow; // Coordinates for arrowhead at origin points in x-dir.
  double fb0 = INTERBASE_DISTANCE*OVERHANG_FRAC; // overhang distance for strand ends
  double TickLength = TICK_LENGTH_FRAC*INTERBASE_DISTANCE; // tick lengths
  double ArcWidth = BACKBONE_WIDTH_FRAC*INTERBASE_DISTANCE; // arc line thickness
  double TickWidth = TICK_WIDTH_FRAC*INTERBASE_DISTANCE; // tick line thickness
  double LineWidth = LINE_WIDTH_FRAC*INTERBASE_DISTANCE; // bp line thinkness
  double d; // d = half width of arrow head
  double gapAngle, b0Angle, fb0Angle; // Angles associated with lengths
  double arrowAngle, arrowEndAngle; // Angles associated with lengths
  double theta,theta0; // Angles along circle plot
  double tol = TOLERANCE; // Newton's method tolerance
  double TranslateDist; // Distance to translate structure (a hack right now)
  int MaxIters = MAX_ITERATIONS; //  Maximum allowed Newton iterations
  char InputFile[MAXLINE]; // The name of the input file
  char OutputFile[MAXLINE]; // The name of the output file
  struct arc *arcs; // array of arc structs
  struct arrow *arrows; // array of arrow structs
  struct tick *ticks; // array of tick structs
  struct line *lines; // array of line structs

  // These are the coordinates for a standard arrow.  To get the arrow correct,
  // these must be rotated and translated.
  d = ARROW_FRAC*a*tan(ARROW_ANGLE);
  StandardArrow = (double **) malloc(6 * sizeof(double *));
  for (i = 0; i < 6; i++) {
    StandardArrow[i] = (double *) malloc(2 * sizeof(double));
  }
  StandardArrow[0][0] = 0.0;
  StandardArrow[0][1] = 0.0;
  StandardArrow[1][0] = -ARROW_FRAC*a;
  StandardArrow[1][1] = -d;
  StandardArrow[2][0] = (1 - 2.0*ARROW_FRAC)*a;
  StandardArrow[2][1] = -d;
  StandardArrow[3][0] = (1 - ARROW_FRAC)*a;
  StandardArrow[3][1] = 0.0;
  StandardArrow[4][0] = (1 - 2.0*ARROW_FRAC)*a;
  StandardArrow[4][1] = d;
  StandardArrow[5][0] = -ARROW_FRAC*a;
  StandardArrow[5][1] = d;
  
  ReadCommandLine(argc,argv,&ShowSequences,&ColorStrands,&LabelStrands,InputFile,
		  OutputFile);

  ReadInputFile(&K,&nStrands,&nBases,&SeqLen,&Perm,&pairs,&nPairs,InputFile);

  // Assign colors
  for (i = 0; i < MAXSTRANDS; i++) {
    ColorKey[i] = -1;
  }
  k = 0;
  for (i = 0; i < nStrands; i++) {
    if (ColorKey[Perm[i]] == -1) {
      ColorKey[Perm[i]] = k;
      k++;
    }
  }

  // Number of types of strands present
  nStrandTypes = k;

  // Allocate memory
  arcs = (struct arc *) malloc(nStrands * sizeof(struct arc));
  arrows = (struct arrow *) malloc(nStrands * sizeof(struct arrow));
  ticks = (struct tick *) malloc(nBases * sizeof(struct tick));
  lines = (struct line *) malloc(nPairs * sizeof(struct line));

  // Get the circle radius
  r = CircleRadius(nBases,nStrands,b0,gap,fb0,a,tol,MaxIters);

  // Hack the translation distance
  TranslateDist = r + 1.1*TickLength;

  // SVG sizes (a hack)
  vbw = 2*r + 2.2*TickLength;
  vbh = vbw;
  SVGHeight = SVG_SIZE;
  SVGWidth = SVGHeight*vbw/vbh;

  // Get the associated angles
  gapAngle = 2.0*asin(gap/2.0/r);
  b0Angle = 2.0*asin(b0/2.0/r);
  fb0Angle = 2.0*asin(fb0/2.0/r);
  arrowAngle = 2.0*asin(a/2.0/r);
  arrowEndAngle = 2.0*asin((1-ARROW_FRAC*a)/2.0/r);
  
  // Enter data into structs for SVG rendering
  k = 0; // Base index
  theta = -PI/2.0 - gapAngle/2.0;
  for (i = 0 ; i < nStrands; i++) {
    theta0 = theta; // Arc start theta

    // Arc start
    arcs[i].x0[0] = r*cos(theta);
    arcs[i].x0[1] = r*sin(theta);
    
    // Advance theta to first tick
    theta -= fb0Angle;

    // Do the ticks
    for (j = 0; j < SeqLen[Perm[i]]-1; j++) {
      ticks[k].x0[0] = r*cos(theta);
      ticks[k].x0[1] = r*sin(theta);
      ticks[k].x1[0] = (r + TickLength)*cos(theta);
      ticks[k].x1[1] = (r + TickLength)*sin(theta);
      k++;

      theta -= b0Angle;
    }

    // Last tick
    ticks[k].x0[0] = r*cos(theta);
    ticks[k].x0[1] = r*sin(theta);
    ticks[k].x1[0] = (r + TickLength)*cos(theta);
    ticks[k].x1[1] = (r + TickLength)*sin(theta);
    k++;

    // End of arc
    theta -= (fb0Angle + arrowEndAngle);
    arcs[i].x1[0] = r*cos(theta);
    arcs[i].x1[1] = r*sin(theta);
    if (theta0-theta > PI) {
      arcs[i].LargeArc = 1;
    }
    else {
      arcs[i].LargeArc = 0;
    }
    arcs[i].Sweep = 1;

    // Arrowhead
    getArrow(arrows,i,StandardArrow,theta,r);
  
    // Advance theta past the arrowhead and the gap
    theta -= (arrowAngle + gapAngle);
  }


  // Put in lines for base pairs
  k = 0;
  for (i = 0; i < nBases; i++) {
    if (pairs[i] > -1 && (j = pairs[i]) > i) {
      lines[k].x0[0] = ticks[i].x0[0];
      lines[k].x0[1] = ticks[i].x0[1];
      lines[k].x1[0] = ticks[j].x0[0];
      lines[k].x1[1] = ticks[j].x0[1];
      k++;
    }
  }

  // Go through and put in the colors
  if (nStrandTypes == 1) { // Make them all black
    for (i = 0; i < nStrands; i++) {
      arcs[i].RGB[0] = 0;
      arcs[i].RGB[1] = 0;
      arcs[i].RGB[2] = 0;
      arrows[i].RGB[0] = 0;
      arrows[i].RGB[1] = 0;
      arrows[i].RGB[2] = 0;
    }
  }
  else {
    for (i = 0; i < nStrands; i++) {
      arcs[i].RGB[0] = StrandColors[ColorKey[Perm[i]]][0];
      arcs[i].RGB[1] = StrandColors[ColorKey[Perm[i]]][1];
      arcs[i].RGB[2] = StrandColors[ColorKey[Perm[i]]][2];
      arrows[i].RGB[0] = StrandColors[ColorKey[Perm[i]]][0];
      arrows[i].RGB[1] = StrandColors[ColorKey[Perm[i]]][1];
      arrows[i].RGB[2] = StrandColors[ColorKey[Perm[i]]][2];
    }
  }
  for (i = 0; i < nPairs; i++) { // pair lines are black for now
    lines[i].RGB[0] = 0;
    lines[i].RGB[1] = 0;
    lines[i].RGB[2] = 0;
  }
  for (i = 0; i < nBases; i++) { // tick lines are black for now
    ticks[i].RGB[0] = 0;
    ticks[i].RGB[1] = 0;
    ticks[i].RGB[2] = 0;
  }

  // Put in arc radii
  for (i = 0; i < nStrands; i++) {
    arcs[i].radius = r;
  }

  // Reflect all coordinates through x-axis
  for (i = 0; i < nStrands; i++) {
    arcs[i].x0[1] *= -1;
    arcs[i].x1[1] *= -1;
    arrows[i].x0[1] *= -1;
    arrows[i].x1[1] *= -1;
    arrows[i].x2[1] *= -1;
    arrows[i].x3[1] *= -1;
    arrows[i].x4[1] *= -1;
    arrows[i].x5[1] *= -1;
  }
  for (i = 0; i < nBases; i++) {
    ticks[i].x0[1] *= -1;
    ticks[i].x1[1] *= -1;
  }
  for (i = 0; i < nPairs; i++) {
    lines[i].x0[1] *= -1;
    lines[i].x1[1] *= -1;
  }

  // Translate
  for (i = 0; i < nStrands; i++) {
    arcs[i].x0[0] += TranslateDist;
    arcs[i].x0[1] += TranslateDist;
    arcs[i].x1[0] += TranslateDist;
    arcs[i].x1[1] += TranslateDist;
    arrows[i].x0[0] += TranslateDist;
    arrows[i].x0[1] += TranslateDist;
    arrows[i].x1[0] += TranslateDist;
    arrows[i].x1[1] += TranslateDist;
    arrows[i].x2[0] += TranslateDist;
    arrows[i].x2[1] += TranslateDist;
    arrows[i].x3[0] += TranslateDist;
    arrows[i].x3[1] += TranslateDist;
    arrows[i].x4[0] += TranslateDist;
    arrows[i].x4[1] += TranslateDist;
    arrows[i].x5[0] += TranslateDist;
    arrows[i].x5[1] += TranslateDist;
  }
  for (i = 0; i < nBases; i++) {
    ticks[i].x0[0] += TranslateDist;
    ticks[i].x0[1] += TranslateDist;
    ticks[i].x1[0] += TranslateDist;
    ticks[i].x1[1] += TranslateDist;
  }
  for (i = 0; i < nPairs; i++) {
    lines[i].x0[0] += TranslateDist;
    lines[i].x0[1] += TranslateDist;
    lines[i].x1[0] += TranslateDist;
    lines[i].x1[1] += TranslateDist;
  }

  // Render the SVG
  SVGHeader(OutputFile,SVGWidth,SVGHeight,vbw,vbh);
  SVGTicks(OutputFile,ticks,nBases,TickWidth);
  SVGLines(OutputFile,lines,nPairs,LineWidth);
  SVGArcs(OutputFile,arcs,arrows,nStrands,ArcWidth);
  SVGEnd(OutputFile);

  for (i = 0; i < 6; i++) {
    free(StandardArrow[i]);
  }
  free(StandardArrow);
  free(Perm);
  free(SeqLen);
  free(pairs);
  free(arcs);
  free(lines);
  free(arrows);
  free(ticks);

  return 0;

}
/* ******************************************************************************** */
/*                                  END MAIN                                        */
/* ******************************************************************************** */



/* ******************************************************************************** */
double CircleRadius(int nBases, int nStrands, double b0, double gap, double fb0,
		    double a, double tol, int MaxIters) {
  /*
    Uses Newton's method to find radius of circle circumscribing bases.
  */

  int iters; // Number of Newton iterations
  double r; // The radius, which is returned
  double f,fprime; // function and it's derivative we're finding a root for
  
  // Find rmin
  if (gap > b0) {
    r = gap/2.0 + tol;
  }
  else {
    r = b0/2.0 + tol;
  }

  // Do Newton's method to get radius
  iters = 0;
  f = tol + 1; // Just to get started
  while (f > tol && iters < MaxIters) {
    f = -PI + nStrands*asin(gap/2.0/r) + (nBases-nStrands)*asin(b0/2.0/r)
      + nStrands*asin(a/2.0/r) + 2*nStrands*asin(fb0/2.0/r);

    fprime = -1.0/2.0/pow(r,2)*(nStrands*gap/sqrt(1-pow(gap/2.0/r,2)) 
				+ (nBases-nStrands)*b0/sqrt(1-pow(b0/2.0/r,2))
				+ nStrands*a/sqrt(1-pow(a/2.0/r,2)));

    r -= f/fprime;
    iters++;
  }

  if (iters == MaxIters) {
    printf("Error: Newton's method failed to converge.\n");
    printf("\nExiting....\n\n");
    exit(ERR_NEWTON_FAIL);
  }  

  return r;

}
/* ******************************************************************************** */

/* ******************************************************************************** */
void SVGHeader(char *OutputFile, double SVGWidth, double SVGHeight, double vbw,
	       double vbh) {

  /*
    Opens output file to write and writes the header.
  */

  FILE *fp;

  if ((fp = fopen(OutputFile,"w")) == NULL) {
    printf("Error opening output file %s\n",OutputFile);
    printf("\nExiting....\n\n");
    exit(ERR_OUTPUTFILE);
  }

  fprintf(fp,"<?xml version=\"1.0\" standalone=\"no\"?>\n");
  fprintf(fp,"<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n");
  fprintf(fp,"  \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
  fprintf(fp,"<svg width=\"%8.7fpx\" height=\"%8.7fpx\" ",SVGWidth,SVGHeight);
  fprintf(fp,"viewBox=\"0 0 %8.7f %8.7f\"\n",vbw,vbh);
  fprintf(fp,"    xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">\n");
  fprintf(fp,"  <desc>NUPACK DOT PLOT</desc>\n\n");

  fclose(fp);
  
}
/* ******************************************************************************** */


/* ******************************************************************************** */
void SVGTicks(char *OutputFile, struct tick *ticks, int nBases, double TickWidth) {

  /*
    Opens output file to append and writes the tick marks.
  */

  int i; // Counter
  FILE *fp; // File we're outputting to

  if ((fp = fopen(OutputFile,"a")) == NULL) {
    printf("Error opening output file %s\n",OutputFile);
    printf("\nExiting....\n\n");
    exit(ERR_OUTPUTFILE);
  }

  fprintf(fp,"  <!-- Begin Tick Marks --> \n");
  fprintf(fp,"  <g>\n");
  for (i = 0; i < nBases; i++) {
    fprintf(fp,"    <path d=\"M %8.7f %8.7f L %8.7f %8.7f\"\n",ticks[i].x0[0],
	    ticks[i].x0[1],ticks[i].x1[0],ticks[i].x1[1]);
    fprintf(fp,"    stroke=\"rgb(%d,%d,%d)\" stroke-width=\"%8.7f\" />\n",
	    ticks[i].RGB[0],ticks[i].RGB[1],ticks[i].RGB[2],TickWidth);
   }
  fprintf(fp,"  </g>\n");
  fprintf(fp,"  <!-- End Tick Marks --> \n\n");

  fclose(fp);
  
}
/* ******************************************************************************** */


/* ******************************************************************************** */
void SVGArcs(char *OutputFile, struct arc *arcs, struct arrow *arrows, int nStrands,
	     double ArcWidth) {

  /*
    Opens output file to append and writes backbone of the circle diagram.
  */

  int i; // Counter
  FILE *fp; // File we're outputting to

  if ((fp = fopen(OutputFile,"a")) == NULL) {
    printf("Error opening output file %s\n",OutputFile);
    printf("\nExiting....\n\n");
    exit(ERR_OUTPUTFILE);
  }

  fprintf(fp,"  <!-- Begin Backbone --> \n");
  fprintf(fp,"  <g>\n");
  for (i = 0; i < nStrands; i++) {
    fprintf(fp,"    <path d=\"M%8.7f,%8.7f\n",arcs[i].x0[0],arcs[i].x0[1]);
    fprintf(fp,"              A%8.7f,%8.7f 0 %d,%d %8.7f,%8.7f\"\n",arcs[i].radius,
	    arcs[i].radius,arcs[i].LargeArc,arcs[i].Sweep,arcs[i].x1[0],
	    arcs[i].x1[1]);
    fprintf(fp,"              stroke=\"rgb(%d,%d,%d)\" stroke-width=\"%8.7f\"\n",
	    arcs[i].RGB[0],arcs[i].RGB[1],arcs[i].RGB[2],ArcWidth);
    fprintf(fp,"              fill=\"none\" />\n");

    fprintf(fp,"    <path d=\"M %8.7f %8.7f\n",arrows[i].x0[0],arrows[i].x0[1]);
    fprintf(fp,"              L %8.7f %8.7f\n",arrows[i].x1[0],arrows[i].x1[1]);
    fprintf(fp,"              L %8.7f %8.7f\n",arrows[i].x2[0],arrows[i].x2[1]);
    fprintf(fp,"              L %8.7f %8.7f\n",arrows[i].x3[0],arrows[i].x3[1]);
    fprintf(fp,"              L %8.7f %8.7f\n",arrows[i].x4[0],arrows[i].x4[1]);
    fprintf(fp,"              L %8.7f %8.7f\n",arrows[i].x5[0],arrows[i].x5[1]);
    fprintf(fp,"              L %8.7f %8.7f\"\n",arrows[i].x0[0],arrows[i].x0[1]);
    fprintf(fp,"          stroke=\"none\" fill=\"rgb(%d,%d,%d)\" />\n",
	    arrows[i].RGB[0],arrows[i].RGB[1],arrows[i].RGB[2]);
  }
  fprintf(fp,"  </g>\n");
  fprintf(fp,"  <!-- End Backbone --> \n\n");

  fclose(fp);
  
}
/* ******************************************************************************** */


/* ******************************************************************************** */
void SVGLines(char *OutputFile, struct line *lines, int nPairs, double LineWidth) {

  /*
    Opens output file to append and writes the lines representing base pairs.
  */

  int i; // Counter
  FILE *fp; // File we're outputting to

  if ((fp = fopen(OutputFile,"a")) == NULL) {
    printf("Error opening output file %s\n",OutputFile);
    printf("\nExiting....\n\n");
    exit(ERR_OUTPUTFILE);
  }

  fprintf(fp,"  <!-- Begin Base Pair Lines --> \n");
  fprintf(fp,"  <g>\n");
  for (i = 0; i < nPairs; i++) {
    fprintf(fp,"    <path d=\"M %8.7f %8.7f L %8.7f %8.7f\"\n",lines[i].x0[0],
	    lines[i].x0[1],lines[i].x1[0],lines[i].x1[1]);
    fprintf(fp,"    stroke=\"rgb(%d,%d,%d)\" stroke-width=\"%8.7f\" />\n",
	    lines[i].RGB[0],lines[i].RGB[1],lines[i].RGB[2],LineWidth);
   }
  fprintf(fp,"  </g>\n");
  fprintf(fp,"  <!-- End Base Pair Lines --> \n\n");

  fclose(fp);
  
}
/* ******************************************************************************** */

/* ******************************************************************************** */
void SVGEnd(char *OutputFile) {

  /*
    Opens output file to append and writes the end of the SVG file.
  */

  FILE *fp;

  if ((fp = fopen(OutputFile,"a")) == NULL) {
    printf("Error opening output file %s\n",OutputFile);
    printf("\nExiting....\n\n");
    exit(ERR_OUTPUTFILE);
  }

  fprintf(fp,"</svg>\n");

  fclose(fp);

}
/* ******************************************************************************** */

/* ******************************************************************************** */
void getArrow(struct arrow *arrows, int i, double **StandardArrow, double theta, 
	      double r) {
  /*
    Takes the standard arrow and rotates it such that it points along the tangent 
    vector to a circle at the point on the circle marked by theta.  The arrow is
    then translated to be at pos = (r cos(theta), r sin(theta)).
  */

  double beta; // The angle for the rotation
  double *pos; 
  double **R; // rotation matrix

  pos = (double *) malloc(2 * sizeof(double));
  R = (double **) malloc (2 * sizeof(double *));
  R[0] = (double *) malloc(2 * sizeof(double));
  R[1] = (double *) malloc(2 * sizeof(double));

  pos[0] = r*cos(theta);
  pos[1] = r*sin(theta);

  beta = theta - PI/2.0;

  R[0][0] = cos(beta);
  R[0][1] = -sin(beta);
  R[1][0] = -R[0][1];
  R[1][1] = R[0][0];

  // Perform the rotations
  MatrixVectorMult(arrows[i].x0,R,StandardArrow[0],2);
  MatrixVectorMult(arrows[i].x1,R,StandardArrow[1],2);
  MatrixVectorMult(arrows[i].x2,R,StandardArrow[2],2);
  MatrixVectorMult(arrows[i].x3,R,StandardArrow[3],2);
  MatrixVectorMult(arrows[i].x4,R,StandardArrow[4],2);
  MatrixVectorMult(arrows[i].x5,R,StandardArrow[5],2);


  // Perform the translations
  arrows[i].x0[0] += pos[0];
  arrows[i].x0[1] += pos[1];
  arrows[i].x1[0] += pos[0];
  arrows[i].x1[1] += pos[1];
  arrows[i].x2[0] += pos[0];
  arrows[i].x2[1] += pos[1];
  arrows[i].x3[0] += pos[0];
  arrows[i].x3[1] += pos[1];
  arrows[i].x4[0] += pos[0];
  arrows[i].x4[1] += pos[1];
  arrows[i].x5[0] += pos[0];
  arrows[i].x5[1] += pos[1];

  free(pos);
  free(R[0]);
  free(R[1]);
  free(R);
}
/* ******************************************************************************** */
