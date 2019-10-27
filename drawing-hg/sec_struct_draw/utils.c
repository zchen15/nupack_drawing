/*
  UTILS.C 

  Utilities for use with 2D secondary structure drawing.

  Justin Bois, Caltech, 14 January 2007
  bois@caltech.edu
*/

#include "SecStructDrawHeader.h" // Secondary structure drawing header file.

/* ******************************************************************************** */
int isPseudoknot(struct base *bases, int nbases) {
  /*
    Returns 1 if structure is pseduoknotted and zero otherwise.
  */

  int i,j,k; // Counters
  int npairs; // Total number of base pairs
  int **pairlist; // List of pairs [i j] style with i < j.

  pairlist = (int **) malloc(nbases *sizeof(int*));
  for (i = 0; i < nbases; i++) {
    pairlist[i] = (int *) malloc(nbases * sizeof(int));
  }
  
  npairs = 0;
  for (i = 0; i < nbases; i++) {
    if (bases[i].pair > i) {
      pairlist[npairs][0] = i;
      pairlist[npairs][1] = bases[i].pair;
      npairs++;
    }
  }

  for (i = 0; i < npairs-1; i++) {
    for (k = i+1; (k < npairs && pairlist[k][0] < pairlist[i][1]); k++) {
      if (pairlist[i][1] < pairlist[k][1]) { // pseudoknot
        for (j = 0; j < nbases; j++) {
          free(pairlist[j]);
        }
        free(pairlist);
        return 1;
      }
    }
  }
  
  for (j = 0; j < nbases; j++) {
    free(pairlist[j]);
  }
  free(pairlist);

  return 0; // No pseudknot

}
/* ******************************************************************************** */


/* ******************************************************************************** */
void BaseColor(int *RGB, char baseType, int white) {
  /*
    Converts a baseType (A, C, G, or T/U) to an RGB color code.
  */

  int index = 0; // Index in BaseColors for corresponding to base

  if (baseType == 'A') {
    index = 0;
  }
  else if (baseType == 'C') {
    index = 1;
  }
  else if (baseType == 'G') {
    index = 2;
  }
  else if (baseType == 'U' || baseType == 'T') {
    index = 3;
  }

  if (white) {
    RGB[0] = BaseColorsWhite[index][0];
    RGB[1] = BaseColorsWhite[index][1];
    RGB[2] = BaseColorsWhite[index][2];
  }
  else {
    RGB[0] = BaseColors[index][0];
    RGB[1] = BaseColors[index][1];
    RGB[2] = BaseColors[index][2];
  }

}
/* ******************************************************************************** */


/* ******************************************************************************** */
void ColorMap(int *RGB, double prob) {
  /*
    Takes a value between 0 and 1 (prob) and converts it to an RGB
    value based on the Matlab jet colormap.  Low values are darker
    blue and higher values are red.
  */

  double red, green, blue; // The RGB channels in range 0 to 1

  // Check to make sure probability values are legit
  if (prob < 0.0) {
    printf("Warning: probability value must satisfy 0 <= prob <= 1, using prob = 0\n");
    prob = 0.0;
  }
  else if (prob > 1.0) {
    printf("Warning: probability value must satisfy 0 <= prob <= 1, using prob = 1\n");
    prob = 1.0;
  }

  // Red channel
  if (prob <= 0.375) {
    red = 0.0;
  }
  else if (prob <= 0.675) {
    red = 4.0*prob - 1.5;
  }
  else if (prob <= 0.875) {
    red = 1.0;
  }
  else { // prob > 0.875
    red = -4.0*prob + 4.5;
  }

  // Green channel
  if (prob <= 0.125 || prob >= 0.875) {
    green = 0.0;
  }
  else if (prob <= 0.375) {
    green = 4.0*prob - 0.5;
  }
  else if (prob <= 0.625) {
    green = 1.0;
  }
  else { // 0.625 < prob < 0.875
    green = -4.0*prob + 3.5;
  }

  // Blue channel
  if (prob <= 0.125) {
    blue = 4.0*prob + 0.5;
  }
  else if (prob <= 0.375) {
    blue = 1.0;
  }
  else if (prob <= 0.625) {
    blue = -4.0*prob + 2.5;
  }
  else { // 0.625 < prob <= 1
    blue = 0.0;
  }

  RGB[0] = (int) 255*red;
  RGB[1] = (int) 255*green;
  RGB[2] = (int) 255*blue;

}
/* ******************************************************************************** */


/* ******************************************************************************** */
void getStandardArrow(double ***StandardArrow, double *ArrowHeight, 
                      double *ArrowRadius) {
  
  int i; // Counter
  
  *ArrowHeight = INTERBASE_DISTANCE*ARROW_HEIGHT_FRAC;
  *ArrowRadius = ARROW_FRAC*(*ArrowHeight)*tan(ARROW_ANGLE);
  
  (*StandardArrow) = (double **) malloc(6 * sizeof(double *));
  for (i = 0; i < 6; i++) {
    (*StandardArrow)[i] = (double *) malloc(2 * sizeof(double));
  }
  
  (*StandardArrow)[0][0] = 0.0;
  (*StandardArrow)[0][1] = 0.0;
  (*StandardArrow)[1][0] = -ARROW_FRAC*(*ArrowHeight);
  (*StandardArrow)[1][1] = -(*ArrowRadius);
  (*StandardArrow)[2][0] = (1 - 2.0*ARROW_FRAC)*(*ArrowHeight);
  (*StandardArrow)[2][1] = -(*ArrowRadius);
  (*StandardArrow)[3][0] = (1 - ARROW_FRAC)*(*ArrowHeight);
  (*StandardArrow)[3][1] = 0.0;
  (*StandardArrow)[4][0] = (1 - 2.0*ARROW_FRAC)*(*ArrowHeight);
  (*StandardArrow)[4][1] = (*ArrowRadius);
  (*StandardArrow)[5][0] = -ARROW_FRAC*(*ArrowHeight);
  (*StandardArrow)[5][1] = (*ArrowRadius);
  
}
/* ******************************************************************************** */


/* ******************************************************************************** */
void getArrow(struct arrow *arrows, int i, double **StandardArrow, double pos[2], 
              double dir[2]) {
  /*
  Takes the standard arrow and rotates it such that it points along
  the vector dir and translates it to pos (i.e, the point of
  attachment of the arrow is (x,y) = (pos[0],pos[1]).

  The coordinate for the arrow is stored entry i of the array of
  arrow structs called arrows.
  */
  
  double theta; // The angle for the rotation
  double **R; // rotation matrix
  
  R = (double **) malloc (2 * sizeof(double *));
  R[0] = (double *) malloc(2 * sizeof(double));
  R[1] = (double *) malloc(2 * sizeof(double));
  
  theta = atan2(dir[1],dir[0]);
  
  R[0][0] = cos(theta);
  R[0][1] = -sin(theta);
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
  
  free(R[0]);
  free(R[1]);
  free(R);
}
/* ******************************************************************************** */


/* ******************************************************************************** */
double CircleRadius(int *bigAngle, double *lengths, int nsegments, double tol,
                    int MaxIters) {
  /*
    Uses Newton's method to find radius of circle circumscribing a
    convex polygon of lengths defined by input vector lengths with
    nsegments elements.
  */

  int i; // Counter
  int iters; // Number of Newton iterations
  double r; // The radius, which is returned
  double f,fprime; // function and it's derivative we're finding a root for
  double f0; // f(r) = f0 + \sum_{i=1}^{nsegmenets} asin(lengths[i]/2/r)

  // Determine whether or not one of the angles is greater than pi
  r = max(lengths,nsegments)/2.0;
  f = -PI;
  for (i = 0; i < nsegments; i++) {
    f += asin(lengths[i]/2.0/r);
  }
  if (f < 0.0) {
    *bigAngle = 1;
    f0 = -PI/2.0;
  }
  else {
    *bigAngle = 0;
    f0 = -PI;
  }


  // Do Newton's method to get radius
  iters = 0;
  r += tol;
  f = tol + 1.0; // Just to get started
  while (f > tol && iters < MaxIters) {
    // Compute f
    f = f0;
    fprime = 0;
    for (i = 0; i < nsegments; i++) {
      f += asin(lengths[i]/2.0/r);
      fprime += lengths[i]/sqrt(1.0 - pow(lengths[i]/2.0/r,2));
    }

    fprime *= (-1.0)/2.0/pow(r,2);

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
void SVGheader(char *OutputFile, double SVGWidth, double SVGHeight, double leftCorner,
               double topCorner, double vbw, double vbh) {

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
  fprintf(fp,"viewBox=\"%8.7f %8.7f %8.7f %8.7f\"\n",leftCorner,topCorner,vbw,vbh);
  fprintf(fp,"    xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">\n");
  fprintf(fp,"  <desc>NUPACK SECONDARY STRUCTURE REPRESENTATION</desc>\n\n");

  fclose(fp);
  
}
/* ******************************************************************************** */


/* ******************************************************************************** */
void SVGarcs(char *OutputFile, struct arc *arcs, struct arrow *arrows, int nStrands,
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
void SVGarrows(char *OutputFile, struct arrow *arrows, struct arrowline *arrowlines, 
               int nStrands, double arrowLineWidth) {

  /*
    Opens output file to append and writes arrows to ladder-hose drawing.
  */

  int i; // Counter
  FILE *fp; // File we're outputting to

  if ((fp = fopen(OutputFile,"a")) == NULL) {
    printf("Error opening output file %s\n",OutputFile);
    printf("\nExiting....\n\n");
    exit(ERR_OUTPUTFILE);
  }

  fprintf(fp,"  <!-- Begin 3' Arrowheads --> \n");
  fprintf(fp,"  <g>\n");
  for (i = 0; i < nStrands; i++) {
    fprintf(fp,"    <path d=\"M %8.7f %8.7f L %8.7f %8.7f\"\n",arrowlines[i].x0[0],
            arrowlines[i].x0[1],arrowlines[i].x1[0],arrowlines[i].x1[1]);
    fprintf(fp,"     stroke=\"rgb(%d,%d,%d)\" stroke-width=\"%8.7f\" />\n",
            arrowlines[i].RGB[0],arrowlines[i].RGB[1],arrowlines[i].RGB[2],
            arrowLineWidth);
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
  fprintf(fp,"  <!-- End 3' Arrowheads --> \n\n");

  fclose(fp);
  
}
/* ******************************************************************************** */


/* ******************************************************************************** */
void SVGbackbone(char *OutputFile, struct backbone *backbones, int nbackbone, 
                 double BackboneWidth) {

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

  fprintf(fp,"  <!-- Begin Lines for Backbone --> \n");
  fprintf(fp,"  <g>\n");
  for (i = 0; i < nbackbone; i++) {
    if (backbones[i].curved) {
      fprintf(fp,"    <path d=\"M%8.7f,%8.7f\n",backbones[i].x0[0],
              backbones[i].x0[1]);
      fprintf(fp,"     A%8.7f,%8.7f 0 %d,%d %8.7f,%8.7f\"\n",
              backbones[i].radius,backbones[i].radius,0,1,backbones[i].x1[0],
              backbones[i].x1[1]);
      fprintf(fp,"     stroke=\"rgb(%d,%d,%d)\" stroke-width=\"%8.7f\"\n",
              backbones[i].RGB[0],backbones[i].RGB[1],backbones[i].RGB[2],
              BackboneWidth);
      fprintf(fp,"     fill=\"none\" />\n");
      
    }
    else {
      fprintf(fp,"    <path d=\"M %8.7f %8.7f L %8.7f %8.7f\"\n",backbones[i].x0[0],
              backbones[i].x0[1],backbones[i].x1[0],backbones[i].x1[1]);
      fprintf(fp,"     stroke=\"rgb(%d,%d,%d)\" stroke-width=\"%8.7f\" />\n",
              backbones[i].RGB[0],backbones[i].RGB[1],backbones[i].RGB[2],
              BackboneWidth);
    }
  }
  fprintf(fp,"  </g>\n");
  fprintf(fp,"  <!-- End Lines for Backbone --> \n\n");

  fclose(fp);
  
}
/* ******************************************************************************** */


/* ******************************************************************************** */
void SVGbasepairs(char *OutputFile, struct basepair *basepairs, int nPairs, 
                  double BPWidth) {

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
    fprintf(fp,"    <path d=\"M %8.7f %8.7f L %8.7f %8.7f\"\n",basepairs[i].x0[0],
            basepairs[i].x0[1],basepairs[i].x1[0],basepairs[i].x1[1]);
    fprintf(fp,"     stroke=\"rgb(%d,%d,%d)\" stroke-width=\"%8.7f\" />\n",
            basepairs[i].RGB[0],basepairs[i].RGB[1],basepairs[i].RGB[2],BPWidth);
   }
  fprintf(fp,"  </g>\n");
  fprintf(fp,"  <!-- End Base Pair Lines --> \n\n");

  fclose(fp);
  
}
/* ******************************************************************************** */


/* ******************************************************************************** */
void SVGdots(char *OutputFile, struct basedot *dots, struct base *bases, 
             int nbases, double DotRadius, int prob, int opacity, int drawBases) {
 
  /*
   Opens output file and writes dots for bases.
  */

  int i; // Counter
  double StrokeWidth;
  double plotProb;
  double m1, m2, b2; // Parameters for function mapping opacity and probability
  FILE *fp; // File we're outputting to
  double xShift, yShift;

  StrokeWidth = DotRadius*BASE_STROKE_FRAC;
  if (drawBases)
   DotRadius*=2;

  if ((fp = fopen(OutputFile,"a")) == NULL) {
   printf("Error opening output file %s\n",OutputFile);
   printf("\nExiting....\n\n");
   exit(ERR_OUTPUTFILE);
  }

  fprintf(fp,"  <!-- Begin Dots for Bases --> \n");
  fprintf(fp,"  <g>\n"); 
  if (prob && opacity) {
    m1 = OPACITY_SLOPE_CHANGE / P_SLOPE_CHANGE;
    m2 = (1.0 - OPACITY_SLOPE_CHANGE) / (1.0 - P_SLOPE_CHANGE);
    b2 = 1.0 - m2;
    for (i = 0; i < nbases; i++) {
     if (dots[i].p < 0.9) {
       plotProb = m1*dots[i].p;
     }
     else {
       plotProb = m2*dots[i].p + b2;
     }
     
     // White dot under it
     fprintf(fp,"    <circle cx=\"%8.7f\" cy=\"%8.7f\" r=\"%8.7f\"\n",
             dots[i].x[0],dots[i].x[1],DotRadius);
     fprintf(fp,"     fill=\"white\" stroke=\"none\" />\n");
     
     // Colored dot with transparency
     fprintf(fp,"    <circle cx=\"%8.7f\" cy=\"%8.7f\" r=\"%8.7f\"\n",
             dots[i].x[0],dots[i].x[1],DotRadius);
     fprintf(fp,"     fill=\"rgb(%d,%d,%d)\" stroke=\"rgb(%d,%d,%d)\"\n",
             dots[i].RGB[0],dots[i].RGB[1],dots[i].RGB[2],
             dots[i].RGB[0],dots[i].RGB[1],dots[i].RGB[2]);
     fprintf(fp,"     stroke-width=\"%8.7f\" fill-opacity=\"%8.7f\" />\n",
             StrokeWidth,plotProb);
    }
   
  }
  else { // No need for opacity, show sequence information or prob colormap
    if (drawBases){
      for (i = 0; i < nbases; i++) {
       fprintf(fp,"    <circle cx=\"%8.7f\" cy=\"%8.7f\" r=\"%8.7f\"\n",
               dots[i].x[0],dots[i].x[1],DotRadius*0.85);
       if (prob)
         fprintf(fp,"       stroke=\"rgb(%d,%d,%d)\" stroke-width=\"%8.7f\" "
                 "stroke-opacity=\"1\"  "
                 "fill=\"white\" fill-opacity=\"1\" />\n",
                 dots[i].RGB[0],dots[i].RGB[1],dots[i].RGB[2], DotRadius*0.15);
       else
         fprintf(fp,"       stroke=\"none\" fill=\"white\" fill-opacity=\"1\" />\n");
      }//stroke-miterlimit=\"4\"

      fprintf(fp,"  </g>\n\n");

      fprintf(fp,"  <!-- Begin Letters for Bases --> \n");
      fprintf(fp,"  <g>\n");
      fprintf(fp,"  <text style=\"font-size:%8.7fpx;font-family:'sans-serif';fill:'black';"
             "fill-opacity:1;font-weight:'bold'\">\n",DotRadius*1.4);
      for (i = 0; i < nbases; i++) {
       //~ fprintf(fp,"    <text x=\"%8.7f\" y=\"%8.7f\" font-family=\"sans-serif\"\n",
       //~ dots[i].x[0]-DotRadius,dots[i].x[1]-DotRadius);
       //~ fprintf(fp,"      style=\"font-size:%8.7fpx;fill:black;fill-opacity:1;"
       //~ "font-weight:bold \">\n",DotRadius*1.4);
       
       switch (bases[i].baseType){
         case 'A':
           xShift=DotRadius*0.55;
           yShift=DotRadius*0.45;
           break;
         case 'T':
           xShift=DotRadius*0.475;
           yShift=DotRadius*0.55;
           break;
         default:
           xShift=DotRadius*0.585;
           yShift=DotRadius*0.5;
       }
       fprintf(fp,"      <tspan  x=\"%8.7f\" y=\"%8.7f\">%c</tspan>\n",
               dots[i].x[0]-xShift,dots[i].x[1]+yShift, bases[i].baseType);
      }
      fprintf(fp,"  </text>\n");
    }
    else{
      for (i = 0; i < nbases; i++) {
       fprintf(fp,"    <circle cx=\"%8.7f\" cy=\"%8.7f\" r=\"%8.7f\"\n",
               dots[i].x[0],dots[i].x[1],DotRadius);
       fprintf(fp,"     fill=\"rgb(%d,%d,%d)\" stroke=\"none\" />\n",
               dots[i].RGB[0],dots[i].RGB[1],dots[i].RGB[2]);
      }
    }
  }

  fprintf(fp,"  </g>\n");
  fprintf(fp,"  <!-- End Dots for Bases --> \n\n");

  fclose(fp);
 
}
/* ******************************************************************************** */


/* ******************************************************************************** */
void SVGColorBar (char *OutputFile, double colorBarx, double colorBarShifty, 
                  double dim, int white) {
  /*
    Draws a color bar with labels on the SVG graphic
  */

  int i; // Counter
  char fontcolor[10]; // Color of the font
  double p; // Fractional value we're going to convert
  double y; // y position of rectangle
  double width; // Colorbar width
  double segheight; // Height of each colorbar segment
  double labelheight; // Height of colorbar labels
  double labelx; // x-coordinate of colorbar labels
  double fontSize; // Label font size
  int *RGB; // Color for segment of colorbar
  FILE *fp; // File we're outputting to
  
  RGB = (int *) malloc(3 * sizeof(int));

  if ((fp = fopen(OutputFile,"a")) == NULL) {
    printf("Error opening output file %s\n",OutputFile);
    printf("\nExiting....\n\n");
    exit(ERR_OUTPUTFILE);
  }

  fprintf(fp,"  <!-- Begin Colorbar --> \n");
  width = dim*COLORBAR_HEIGHT_FRAC/COLORBAR_HEIGHT_TO_WIDTH;
  segheight = dim*COLORBAR_HEIGHT_FRAC/N_COLORBAR_COLORS;
  y = dim*(1.0 - COLORBAR_HEIGHT_FRAC)/2.0 + colorBarShifty;
  fprintf(fp,"  <g>\n");
  for (i = 0; i < N_COLORBAR_COLORS; i++) {
    p = 1.0 - ((double) i) / ((double) N_COLORBAR_COLORS);
    ColorMap(RGB,p);    
    fprintf(fp,"    <rect x=\"%8.7f\" y = \"%8.7f\" ",colorBarx,y);
    fprintf(fp,"width=\"%8.7f\" height=\"%8.7f\"\n",width,segheight);
    fprintf(fp,"    fill=\"rgb(%d,%d,%d)\" stroke=\"none\"/>\n",RGB[0],RGB[1],RGB[2]);
    y += segheight;
  }
  fprintf(fp,"  <!-- End Colorbar --> \n\n");


  fprintf(fp,"  <!-- Begin Colorbar Labels --> \n");
  fontSize = dim*COLORBAR_HEIGHT_FRAC*COLORBAR_LABEL_FONT_SIZE_FRAC;
  labelx = colorBarx + width + COLORBAR_LABEL_OFFSET_FRAC*width;
  labelheight = dim*COLORBAR_HEIGHT_FRAC/(N_COLORBAR_LABELS-1);
  y = dim*(1.0 - COLORBAR_HEIGHT_FRAC)/2.0 + colorBarShifty
    + COLOR_BAR_VERT_FUDGE_FACT_FRAC*fontSize;
  if (white) {
    strcpy(fontcolor,"white");
  }
  else {
    strcpy(fontcolor,"black");
  }
  for (i = 0; i < N_COLORBAR_LABELS; i++) {
    p = 1.0 - ((double) i) / ((double) N_COLORBAR_LABELS - 1.0);
    fprintf(fp,"    <text x=\"%8.7f\" y=\"%8.7f\" ",labelx,y);
    fprintf(fp,"font-size=\"%g\" font-family=\"%s\" fill=\"%s\" ",
            fontSize,FONT,fontcolor);
    fprintf(fp,"text-anchor=\"start\">%2.1f</text>\n",p);
    y += labelheight;
  }
  fprintf(fp,"  </g>\n");
  fprintf(fp,"  <!-- End Colorbar Labels --> \n\n");

  free(RGB);
  fclose(fp);
}
/* ******************************************************************************** */


/* ******************************************************************************** */
void SVGBaseKey(char *OutputFile, double dotRadius, double baseKeyx, double dim, 
                int thymine, int white) {
  /*
    Draws a color bar with labels on the SVG graphic
  */

  double fontSize; // Label font size
  double labelx; // Label x position
  double y; // y-position
  double keySep; // Separation between entries in the base key
  FILE *fp; // File we're outputting to
  
  if ((fp = fopen(OutputFile,"a")) == NULL) {
    printf("Error opening output file %s\n",OutputFile);
    printf("\nExiting....\n\n");
    exit(ERR_OUTPUTFILE);
  }

  /*
    This is the old way.  It requires BASE_KEY_VERT_FUDGE_FACT_FRAC to be 1.2/3.0.

  dotRadius = max2(dim*BASE_KEY_DOT_RADIUS_FRAC,dotRadius);
  keySep = dotRadius*BASE_KEY_SEP_FRAC;
  fontSize = dotRadius*BASE_KEY_FONT_SIZE_FRAC;
  labelx = baseKeyx + dotRadius*BASE_KEY_LABEL_OFFSET_FRAC;
  */

  // New spacing with constant font size, improved Feng Shui.
  dotRadius = BASE_KEY_DOT_RADIUS_FRAC*dim;
  fontSize = BASE_KEY_FONT_SIZE_FRAC*dotRadius;
  keySep = dotRadius*BASE_KEY_SEP_FRAC;
  labelx = baseKeyx + dotRadius*BASE_KEY_LABEL_OFFSET_FRAC;

  fprintf(fp,"  <!-- Begin Base Key Dots --> \n");
  fprintf(fp,"  <g>\n");
  // Dot for A
  y = (dim - 3.0*keySep)/2.0;
  fprintf(fp,"   <circle cx=\"%8.7f\" cy=\"%8.7f\" r=\"%8.7f\"\n",
          baseKeyx,y,dotRadius);
  if (white) {
    fprintf(fp,"    fill=\"rgb(%d,%d,%d)\" stroke=\"none\" />\n",
            BaseColorsWhite[0][0],BaseColorsWhite[0][1],BaseColorsWhite[0][2]);
  }
  else {
    fprintf(fp,"    fill=\"rgb(%d,%d,%d)\" stroke=\"none\" />\n",
            BaseColors[0][0],BaseColors[0][1],BaseColors[0][2]);
  }

  // Dot for C
  y += keySep;
  fprintf(fp,"    <circle cx=\"%8.7f\" cy=\"%8.7f\" r=\"%8.7f\"\n",
          baseKeyx,y,dotRadius);
  if (white) {
    fprintf(fp,"    fill=\"rgb(%d,%d,%d)\" stroke=\"none\" />\n",
            BaseColorsWhite[1][0],BaseColorsWhite[1][1],BaseColorsWhite[1][2]);
  }
  else {
    fprintf(fp,"    fill=\"rgb(%d,%d,%d)\" stroke=\"none\" />\n",
            BaseColors[1][0],BaseColors[1][1],BaseColors[1][2]);
  }

  // Dot for G
  y += keySep;
  fprintf(fp,"    <circle cx=\"%8.7f\" cy=\"%8.7f\" r=\"%8.7f\"\n",
          baseKeyx,y,dotRadius);
  if (white) {
    fprintf(fp,"    fill=\"rgb(%d,%d,%d)\" stroke=\"none\" />\n",
            BaseColorsWhite[2][0],BaseColorsWhite[2][1],BaseColorsWhite[2][2]);
  }
  else {
    fprintf(fp,"    fill=\"rgb(%d,%d,%d)\" stroke=\"none\" />\n",
            BaseColors[2][0],BaseColors[2][1],BaseColors[2][2]);
  }

  // Dot for T/U
  y += keySep;
  fprintf(fp,"    <circle cx=\"%8.7f\" cy=\"%8.7f\" r=\"%8.7f\"\n",
          baseKeyx,y,dotRadius);
  if (white) {
    fprintf(fp,"    fill=\"rgb(%d,%d,%d)\" stroke=\"none\" />\n",
            BaseColorsWhite[3][0],BaseColorsWhite[3][1],BaseColorsWhite[3][2]);
  }
  else {
    fprintf(fp,"    fill=\"rgb(%d,%d,%d)\" stroke=\"none\" />\n",
            BaseColors[3][0],BaseColors[3][1],BaseColors[3][2]);
  }

  fprintf(fp,"  <!-- End Base Key Dots --> \n\n");


  fprintf(fp,"  <!-- Begin Base Key Labels --> \n");
  // A label
  y = (dim - 3.0*keySep)/2.0 + BASE_KEY_VERT_FUDGE_FACT_FRAC*dotRadius;
  fprintf(fp,"    <text x=\"%8.7f\" y=\"%8.7f\" ",labelx,y);
  if (white) {
    fprintf(fp,"font-size=\"%g\" font-family=\"%s\" fill=\"rgb(%d,%d,%d)\" ",
            fontSize,FONT,BaseColorsWhite[0][0],BaseColorsWhite[0][1],
            BaseColorsWhite[0][2]);
  }
  else {
    fprintf(fp,"font-size=\"%g\" font-family=\"%s\" fill=\"rgb(%d,%d,%d)\" ",
            fontSize,FONT,BaseColors[0][0],BaseColors[0][1],BaseColors[0][2]);
  }
  fprintf(fp,"text-anchor=\"start\">A</text>\n");

  // C label
  y += keySep;
  fprintf(fp,"    <text x=\"%8.7f\" y=\"%8.7f\" ",labelx,y);
  if (white) {
    fprintf(fp,"font-size=\"%g\" font-family=\"%s\" fill=\"rgb(%d,%d,%d)\" ",
            fontSize,FONT,BaseColorsWhite[1][0],BaseColorsWhite[1][1],
            BaseColorsWhite[1][2]);
  }
  else {
    fprintf(fp,"font-size=\"%g\" font-family=\"%s\" fill=\"rgb(%d,%d,%d)\" ",
            fontSize,FONT,BaseColors[1][0],BaseColors[1][1],BaseColors[1][2]);
  }
  fprintf(fp,"text-anchor=\"start\">C</text>\n");

  // G label
  y += keySep;
  fprintf(fp,"    <text x=\"%8.7f\" y=\"%8.7f\" ",labelx,y);
  if (white) {
    fprintf(fp,"font-size=\"%g\" font-family=\"%s\" fill=\"rgb(%d,%d,%d)\" ",
            fontSize,FONT,BaseColorsWhite[2][0],BaseColorsWhite[2][1],
            BaseColorsWhite[2][2]);
  }
  else {
    fprintf(fp,"font-size=\"%g\" font-family=\"%s\" fill=\"rgb(%d,%d,%d)\" ",
            fontSize,FONT,BaseColors[2][0],BaseColors[2][1],BaseColors[2][2]);
  }
  fprintf(fp,"text-anchor=\"start\">G</text>\n");

  // T/U label
  y += keySep;
  fprintf(fp,"    <text x=\"%8.7f\" y=\"%8.7f\" ",labelx,y);
  if (white) {
    fprintf(fp,"font-size=\"%g\" font-family=\"%s\" fill=\"rgb(%d,%d,%d)\" ",
            fontSize,FONT,BaseColorsWhite[3][0],BaseColorsWhite[3][1],
            BaseColorsWhite[3][2]);
  }
  else {
    fprintf(fp,"font-size=\"%g\" font-family=\"%s\" fill=\"rgb(%d,%d,%d)\" ",
            fontSize,FONT,BaseColors[3][0],BaseColors[3][1],BaseColors[3][2]);
  }
  if (thymine) {
    fprintf(fp,"text-anchor=\"start\">T</text>\n");
  }
  else {
    fprintf(fp,"text-anchor=\"start\">U</text>\n");
  }

  fprintf(fp,"  </g>\n");
  fprintf(fp,"  <!-- End Base Key Labels --> \n\n");

  fclose(fp);
}
/* ******************************************************************************** */


/* ******************************************************************************** */
void SVGEnergy (char *OutputFile, double freeEnergy, char *units, 
                char *freeEnergyPreamble, double xpos, double ypos, double dim,
                int white) {
  /*
    Opens output file to append and write the free energy in the lower
    left corner of the SVG image.
  */

  double fontSize; // Label font size
  FILE *fp; // File we're outputting to
  
  if ((fp = fopen(OutputFile,"a")) == NULL) {
    printf("Error opening output file %s\n",OutputFile);
    printf("\nExiting....\n\n");
    exit(ERR_OUTPUTFILE);
  }

  fontSize = dim*FREE_ENERGY_FONT_SIZE_FRAC;

  fprintf(fp,"  <!-- Begin Free Energy --> \n");
  fprintf(fp,"  <g>\n");

  fprintf(fp,"    <text x=\"%8.7f\" y=\"%8.7f\" ",xpos,ypos);
  if (white) {
    fprintf(fp,"font-size=\"%g\" font-family=\"%s\" fill=\"white\" ",fontSize,FONT);
  }
  else {
    fprintf(fp,"font-size=\"%g\" font-family=\"%s\" fill=\"black\" ",fontSize,FONT);
  }
  if (strlen(freeEnergyPreamble) > 0) {
    fprintf(fp,"text-anchor=\"start\">%s %.2f %s</text>\n",freeEnergyPreamble,
            freeEnergy,units);
  }
  else {
    fprintf(fp,"text-anchor=\"start\">%.2f %s</text>\n",freeEnergy,units);    
  }
  fprintf(fp,"  </g>\n");
  fprintf(fp,"  <!-- End Free Energy --> \n\n"); 

  fclose(fp);

}
/* ******************************************************************************** */


/* ******************************************************************************** */
void SVGend(char *OutputFile) {

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
void getRfromTheta(double **R, double theta) {
  /*
    Generates rotation matrix R from angle theta;
  */

  R[0][0] = cos(theta);
  R[0][1] = -sin(theta);
  R[1][0] = -R[0][1];
  R[1][1] = R[0][0];

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


/* ******************************************************************************** */
double norm(double *ar, int len) {
  /*
    Computes the norm of an array of double of length len.
  */

  return sqrt(dot(ar,ar,len));

}
/* ******************************************************************************** */


/* ******************************************************************************** */
double max(double *ar, int len) {
  /*
    Returns the maximum entry in an array of doubles of length len.
  */

  double MaxEntry; // The max entry
  int i; // Counter

  MaxEntry = ar[0];
  for (i = 0; i < len; i++) {
    if (ar[i] > MaxEntry) {
      MaxEntry = ar[i];
    }
  }

  return MaxEntry;

}
/* ******************************************************************************** */


/* ******************************************************************************** */
int doubleGreaterCmp(const void *p1, const void *p2) {
  /*
    Compares two doubles that are in an array of doubles being sorted using qsort.
    This sorts from biggest to smallest.
  */

  // Need to use pointers to double to access values
  const double *a1 = p1;
  const double *a2 = p2;

  if (*a1 < *a2) {
    return 1;
  }
  else if (*a1 == *a2) {
    return 0;
  }
  else {
    return -1;
  }
}
/* ******************************************************************************** */


/* ******************************************************************************** */
double min2(double a, double b) {
  /*
    Returns the min of two doubles, a and b
  */

  if (a < b) {
    return a;
  }
  else {
    return b;
  }
}
/* ******************************************************************************** */


/* ******************************************************************************** */
double max2(double a, double b) {
  /*
    Returns the min of two doubles, a and b
  */

  if (a > b) {
    return a;
  }
  else {
    return b;
  }
}
/* ******************************************************************************** */


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
double pointpointdist(double *x1, double *x2, int n) {
  /*
    Returns the distance between points x1 and x2 in n-space.
  */
  
  int i; // Counter
  double *diffvec; // Vector containing x2-x1
  double dist; // The value that is returned

  diffvec = (double *) malloc(n * sizeof(double));
  
  for (i = 0; i < n; i++) {
    diffvec[i] = x2[i] - x1[i];
  }

  dist = norm(diffvec,n);

  free(diffvec);
  
  return dist;

}
/* ******************************************************************************** */


/* ******************************************************************************** */
double pointlinedist(double *x, double *line, int n) {
  /*
    Returns the distance between the line segment line and point x in n-space.
    
    line must have 2*n entries.  The first n entries are the
    coordinates for the first point of the line segment and the second
    n entries are for the second point.  The end points of the lines
    must not be equal (this is not checked).
  */
  
  int i; // Counter
  double t; // Distance along line where extremum distance occurs
  double *x0; // A point at the end of the lines to check extremum distances
  double *diffvec1; // Vector containing x-x1
  double *diffvec2; // Vector containing x2-x1
  double dist1; // Distance from point to one end of line segment
  double dist2; // Distance of point to other end of line segment
  double dist; // The returned point-line distance

  diffvec1 = (double *) malloc(n * sizeof(double));
  diffvec2 = (double *) malloc(n * sizeof(double));
  x0 = (double *) malloc(n * sizeof(double));
 
  for (i = 0; i < n; i++) {
    diffvec1[i] = x[i] - line[i];
    diffvec2[i] = line[i+n] - line[i];
  }

  t = dot(diffvec1,diffvec2,n)/dot(diffvec2,diffvec2,n);

  if (t > 0.0 && t < 1.0) {
    for (i = 0; i < n; i++) {
      x0[i] = x[i] - (line[i] + t*diffvec2[i]);
    }
    dist = norm(x0,2);
  }
  else {
    for (i = 0; i < n; i++) {
      x0[i] = line[i];
    }
    dist1 = pointpointdist(x0,x,2);

    for (i = 0; i < n; i++) {
      x0[i] = line[i+n];
    }
    dist2 = pointpointdist(x0,x,2);
    dist = min2(dist1,dist2);
  }

  free(diffvec1);
  free(diffvec2);
  free(x0);
  
  return dist;

}
/* ******************************************************************************** */


/* ******************************************************************************** */
int linecross2d(double line1[4], double line2[4]) {
  /*
    Returns 1 if line1 and line 2 cross and 0 otherwise.
    
    lines must have 4 entries.  The first 2 entries are the
    coordinates for the first point of the line segment and the second
    n entries are for the second point.  The end points of the lines
    must not be equal (this is not checked).
  */
  
  double t1, t2; // Distance along line where they might cross
  double u1[2], u2[2]; // Unit vectors pointing along lines
  double l1,l2; // line lengths
  double diffvec1[2]; // Vector connecting line end points
  double diffvec2[2]; // Vector connecting line end points

  diffvec1[0] = line1[2] - line1[0];
  diffvec1[1] = line1[3] - line1[1];
  diffvec2[0] = line2[2] - line2[0];
  diffvec2[1] = line2[3] - line2[1];
  l1 = norm(diffvec1,2);
  l2 = norm(diffvec2,2);
  u1[0] = diffvec1[0]/l1;
  u1[1] = diffvec1[1]/l1;
  u2[0] = diffvec2[0]/l2;
  u2[1] = diffvec2[1]/l2;
  
  t1 = (u2[1]*(line2[0]-line1[0]) - u2[0]*(line2[1]-line1[1])) /
    (l1*(u1[0]*u2[1] - u1[1]*u2[0]));

  t2 = (u1[1]*(line2[0]-line1[0]) - u1[0]*(line2[1]-line1[1])) /
    (l2*(u1[0]*u2[1] - u1[1]*u2[0]));

  if (t1 >= 0.0 && t1 <= 1.0 && t2 >= 0.0 && t2 <= 1.0) { // They cross
    return 1;
  }

  return 0;

}
/* ******************************************************************************** */

/* ******************************************************************************** */
void printBases(struct base *bases,int nbases) {
  /*
    Utility function to print the content of the array of base structs.

    All entries except position, probability, and baseType must be filled in.

    Useful in debugging.
  */

  int i; //Counter

  printf("Bases:\nbase\tpair\tnext\tprev\tstrand\tloop\thelix\n");
  for (i = 0; i < nbases; i++) {
    printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
            i,bases[i].pair,bases[i].nextConnection,bases[i].prevConnection,
            bases[i].strandID,bases[i].loop,bases[i].helix);
  }

}
/* ******************************************************************************** */


/* ******************************************************************************** */
void printLoops(struct loop *loops,int nloops) {
  /*
    Utility function to print the content of the array of loop structs.
  */

  int i,j; //Counters

  printf("Loops:\nloop\tnbases\tnnicks\tnhel\tnphan\tbases\n");
  for (i = 0; i < nloops; i++) {
    printf("%d\t%d\t%d\t%d\t%d\t",
            i,loops[i].nbases,loops[i].nnicks,loops[i].nhelices,loops[i].nphantom);
    for (j = 0; j < loops[i].nbases - 1; j++) {
      printf("%d\t",loops[i].bases[j]);
    }
    printf("%d\n",loops[i].bases[loops[i].nbases-1]);
  }

}
/* ******************************************************************************** */


/* ******************************************************************************** */
void printHelices(struct helix *helices,int nhelices) {
  /*
    Utility function to print the content of the array of helix structs.
  */

  int i; //Counters

  printf("Helices:\nhelix\tnbases\tcorners\n");
  for (i = 0; i < nhelices; i++) {
    printf("%d\t%d\t%d\t%d\t%d\t%d\n",
           i,helices[i].nbases,helices[i].corners[0],helices[i].corners[1],
           helices[i].corners[2],helices[i].corners[3]);

  }

}
/* ******************************************************************************** */
