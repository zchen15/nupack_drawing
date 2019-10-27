/*
  Renders circular polymer graph.

  See ssdraw.c for comments.

  Justin Bois, 29 September, 2007
*/

#include "SecStructDrawHeader.h" // Header file for Concentrations

/* ******************************************************************************** */
void RenderCirclePolyGraph(struct base *bases, int nbases, double **StandardArrow,
                           double ArrowHeight, double ArrowRadius, int ColorStrands,
                           int white, int prob, int opacity, int showFreeEnergy, 
                           double freeEnergy, char *units, char *freeEnergyPreamble, 
                           int key, int thymine, int keyBuffer, int colorbarBuffer, 
                           int freeEnergyBuffer, int squareBuffer, int frameBuffer, 
                           char *OutputFile) {
                             

  int i,j,k; // counters
  int nStrandTypes=1; // Number of different sequences on circle diagram (initialized
                      // to 1 for now)
  int strandLengths[MAXSTRANDS]; // Strand lengths in structure
  int perm[MAXSTRANDS]; // The permutation
  int nPairs; // Total number of base pairs
  int nStrands; // Number of strands in complex
  int bigAngle; // Actually a dummy variable
  double arrowPos[2]; // The position of the arrowhead
  double arrowDir[2]; // The direction the arrowhead points
  double *segLengths; // The lengths of the segments the circle ensribes
  double SVGWidth = SVG_SIZE;  // Dimensions of SVG graphic
  double SVGHeight = SVG_SIZE; // Dimensions of SVG graphic
  double leftCorner, topCorner; // The top left corner coordinates of the view box
  double dim; // The dimension of the image of the circe plot itself
  double totalImageWidth, totalImageHeight; // Wifth and height of drawing with labels
  double extraTranslateX, extraTranslateY; // Extra distance to translate free energy
                                           // and colorbar/key.
  double vbw, vbh; // View box width and height
  double r; // The radius of the circle
  double b0 = INTERBASE_DISTANCE; // inter-base length
  double gap = GAP_FACT * INTERBASE_DISTANCE; // Gap size
  double fb0 = INTERBASE_DISTANCE*OVERHANG_FRAC; // overhang distance for strand ends
  double ArcWidth = BACKBONE_WIDTH_FRAC*INTERBASE_DISTANCE; // arc line thickness
  double BPWidth = BP_WIDTH_FRAC*INTERBASE_DISTANCE; // bp line thinkness
  double DotRadius = BASE_RADIUS_FRAC*INTERBASE_DISTANCE; // radius for base dots
  double gapAngle, b0Angle, fb0Angle; // Angles associated with lengths
  double arrowAngle, arrowEndAngle; // Angles associated with lengths
  double theta,theta0; // Angles along circle plot
  double TranslateDist[2]; // Distance to translate the images.
  double cbBuffer,kBuffer; // Buffer sizes for colorbar and key
  double fBuffer,sBuffer; // Buffer sizes for free energy and square buffering
  struct arc *arcs; // array of arc structs
  struct arrow *arrows; // array of arrow structs
  struct basedot *dots; // array of tick structs
  struct basepair *basepairs; // array of line structs


  // Go through bases and pick out convenient information
  nStrands = 0;
  nPairs = 0;
  strandLengths[nStrands] = 0;
  for (i = 0; i < nbases; i++) {
    (strandLengths[nStrands])++;
    if (bases[i].pair > -1) {
      nPairs++;
    }
    if (bases[i].nextConnection == 0) {
      perm[nStrands++] = bases[i].strandID;
      strandLengths[nStrands] = 0;
    }
  }
  nPairs /= 2;


  // Allocate memory
  arcs = (struct arc *) malloc(nStrands * sizeof(struct arc));
  arrows = (struct arrow *) malloc(nStrands * sizeof(struct arrow));
  dots = (struct basedot *) malloc(nbases * sizeof(struct basedot));
  basepairs = (struct basepair *) malloc(nPairs * sizeof(struct basepair));


  // Build the segment lengths
  segLengths = (double *) malloc((5*nStrands + nbases - nStrands) * sizeof(double));
  for (i = 0; i < nStrands; i++) {
    segLengths[i] = gap;
    segLengths[i+nStrands] = fb0;
    segLengths[i+2*nStrands] = fb0;
    segLengths[i+3*nStrands] = ArrowHeight;
    segLengths[i+4*nStrands] = (1-ARROW_FRAC*ArrowHeight);
  }
  for (i = 5*nStrands; i < 5*nStrands + nbases - nStrands; i++) {
    segLengths[i] = b0;
  }


  // Compute the radius of the circle diagram
  r = CircleRadius(&bigAngle,segLengths,5*nStrands+nbases-nStrands,TOLERANCE,
		   MAX_ITERATIONS);


  // Get the associated angles
  gapAngle = 2.0*asin(gap/2.0/r);
  b0Angle = 2.0*asin(b0/2.0/r);
  fb0Angle = 2.0*asin(fb0/2.0/r);
  arrowAngle = 2.0*asin(ArrowHeight/2.0/r);
  arrowEndAngle = 2.0*asin((1-ARROW_FRAC*ArrowHeight)/2.0/r);


  // Get the positions of the bases
  theta = -PI/2.0 - gapAngle/2.0 - fb0Angle;
  i = 0;
  for (j = 0; j < nStrands; j++) {
    for (k = 0; k < strandLengths[j]; k++) {
      bases[i].x[0] = r*cos(theta);
      bases[i].x[1] = r*sin(theta);
      theta -= b0Angle;
      i++;
    }
    theta -= (2*fb0Angle + arrowEndAngle + arrowAngle + gapAngle - b0Angle);
  }

 
  /* ***********  Enter data into structs for SVG rendering *************** */
  // First the coordinates for the arcs and arrowheads
  theta = -PI/2.0 - gapAngle/2.0;
  for (i = 0 ; i < nStrands; i++) {
    arcs[i].radius = r; // Arc radius

    theta0 = theta; // Arc start theta

    // Arc start
    arcs[i].x0[0] = r*cos(theta);
    arcs[i].x0[1] = r*sin(theta);
    
    // End of arc
    theta -= ((strandLengths[i]-1)*b0Angle + 2*fb0Angle + arrowEndAngle);
    arcs[i].x1[0] = r*cos(theta);
    arcs[i].x1[1] = r*sin(theta);
    if (theta0-theta > PI) {
      arcs[i].LargeArc = 1;
    }
    else {
      arcs[i].LargeArc = 0;
    }
    arcs[i].Sweep = 1;

    // Arrowhead, direction is the negative tangent to the circle
    arrowDir[0] = sin(theta);
    arrowDir[1] = -cos(theta);
    arrowPos[0] = r*cos(theta);
    arrowPos[1] = r*sin(theta);
    getArrow(arrows,i,StandardArrow,arrowPos,arrowDir);
  
    // Advance theta past the arrowhead and the gap
    theta -= (arrowAngle + gapAngle);
  }

  // Color the strands and arrows
  if (!ColorStrands || nStrandTypes == 1) { // Make them all default
    if (white) {
      for (i = 0; i < nStrands; i++) {
	arcs[i].RGB[0] = StrandColorWhite[0];
	arcs[i].RGB[1] = StrandColorWhite[1];
	arcs[i].RGB[2] = StrandColorWhite[2];
	arrows[i].RGB[0] = StrandColorWhite[0];
	arrows[i].RGB[1] = StrandColorWhite[1];
	arrows[i].RGB[2] = StrandColorWhite[2];
      }
    }
    else{
      for (i = 0; i < nStrands; i++) {
	arcs[i].RGB[0] = DefaultStrandColor[0];
	arcs[i].RGB[1] = DefaultStrandColor[1];
	arcs[i].RGB[2] = DefaultStrandColor[2];
	arrows[i].RGB[0] = DefaultStrandColor[0];
	arrows[i].RGB[1] = DefaultStrandColor[1];
	arrows[i].RGB[2] = DefaultStrandColor[2];
      }
    }
  }
  else {
    for (i = 0; i < nStrands; i++) {
      arcs[i].RGB[0] = StrandColors[perm[i]][0];
      arcs[i].RGB[1] = StrandColors[perm[i]][1];
      arcs[i].RGB[2] = StrandColors[perm[i]][2];
      arrows[i].RGB[0] = StrandColors[perm[i]][0];
      arrows[i].RGB[1] = StrandColors[perm[i]][1];
      arrows[i].RGB[2] = StrandColors[perm[i]][2];
    }
  }


  // Dots for the bases
  // Put bases in dot structs for writing
  for (i = 0; i < nbases; i++) {
    dots[i].x[0] = bases[i].x[0];
    dots[i].x[1] = bases[i].x[1];
  }
  if (prob && !opacity) { // Use RGB colors for probability
    for (i = 0; i < nbases; i++) {
      ColorMap(dots[i].RGB,bases[i].p);
    }
  }
  else { // Use sequence information for colors of bases
    for (i = 0; i < nbases; i++) {
      BaseColor(dots[i].RGB,bases[i].baseType,white);
    }
  }
  if (prob) {
    for (i = 0; i < nbases; i++) {
      dots[i].p = bases[i].p;
    }
  }
  else {
    for (i = 0; i < nbases; i++) {
      dots[i].p = 1.0;
    }
  }


  // Lines for base pairs
  j = 0;
  for (i = 0; i < nbases; i++) {
    if (bases[i].pair > -1 && i < bases[i].pair) {
      basepairs[j].x0[0] = bases[i].x[0];
      basepairs[j].x0[1] = bases[i].x[1];
      basepairs[j].x1[0] = bases[bases[i].pair].x[0];
      basepairs[j].x1[1] = bases[bases[i].pair].x[1];
      if (white) {
	basepairs[j].RGB[0] = BasePairColorWhite[0];
	basepairs[j].RGB[1] = BasePairColorWhite[1];
	basepairs[j].RGB[2] = BasePairColorWhite[2];
      }
      else {
	basepairs[j].RGB[0] = BasePairColor[0];
	basepairs[j].RGB[1] = BasePairColor[1];
	basepairs[j].RGB[2] = BasePairColor[2];
      }
      j++;
    }
  }
  /* ********************************************************************** */


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
  for (i = 0; i < nbases; i++) {
    dots[i].x[1] *= -1;
  }
  for (i = 0; i < nPairs; i++) {
    basepairs[i].x0[1] *= -1;
    basepairs[i].x1[1] *= -1;
  }

  // Translate (not necessary, but for some reason Illustrator likes it)
  TranslateDist[0] = r + DotRadius + EDGE_BUFFER_FRAC*(2*r+2*DotRadius);
  TranslateDist[1] = TranslateDist[0];

  // Calculate width and height of image
  dim = 2.0*TranslateDist[0];
  totalImageWidth = dim;
  totalImageHeight = dim;
  leftCorner = 0.0;
  topCorner = 0.0;

  // Calculate the buffers on the drawings
  cbBuffer = COLORBAR_OFFSET_FRAC*dim
    + dim*COLORBAR_HEIGHT_FRAC/COLORBAR_HEIGHT_TO_WIDTH*(1+COLORBAR_LABEL_OFFSET_FRAC)
    + COLORBAR_OVERSHOOT_FRAC*dim*COLORBAR_HEIGHT_FRAC*COLORBAR_LABEL_FONT_SIZE_FRAC;

  /*
    This is the old way.  See comments in utils.c.

  kBuffer =  BASE_KEY_OFFSET_FRAC*dim
      + dim*BASE_KEY_DOT_RADIUS_FRAC*(1 + BASE_KEY_OFFSET_FRAC)
      + BASE_KEY_OVERSHOOT_FRAC*max2(dim*BASE_KEY_DOT_RADIUS_FRAC,DotRadius)*BASE_KEY_FONT_SIZE_FRAC;
  */

  // New key buffer with better Feng Shui
  kBuffer =  BASE_KEY_OFFSET_FRAC*dim 
    + BASE_KEY_DOT_RADIUS_FRAC*dim*(2.0 + BASE_KEY_SEP_FRAC)
    + BASE_KEY_FONT_SIZE_FRAC_TOTAL*dim*(1 + BASE_KEY_OVERSHOOT_FRAC);

  fBuffer = dim*(FREE_ENERGY_BOTTOM_BUFFER_FRAC + FREE_ENERGY_TOP_BUFFER_FRAC
		 + FREE_ENERGY_HEIGHT_FRAC*FREE_ENERGY_FONT_SIZE_FRAC);

  if (colorbarBuffer && (!keyBuffer || cbBuffer >= kBuffer)) { // Include colorbar
    totalImageWidth += cbBuffer;
  }
  else if (keyBuffer && (!colorbarBuffer || cbBuffer < kBuffer)) { // Include base key
    totalImageWidth += kBuffer;
  }
  if (freeEnergyBuffer) { // free energy is shown
    totalImageHeight += fBuffer;
  }
  if (squareBuffer) { // Square buffer is max of all buffers on all sides
    sBuffer = max2(fBuffer,max2(cbBuffer,kBuffer));
    // All buffers are equal
    cbBuffer = sBuffer;
    kBuffer = sBuffer;
    fBuffer = sBuffer;
    totalImageHeight = dim + 2.0*sBuffer;
    totalImageWidth = totalImageHeight;
  }
  else { // Set the square buffer to zero otherwise
    sBuffer = 0.0;
  }
  vbh = max2(totalImageHeight,totalImageWidth);
  vbw = vbh;
  SVGHeight = vbh;
  SVGWidth = vbw;


  // Extra translation needed for buffers
  if (squareBuffer) {
    extraTranslateX = sBuffer;
    extraTranslateY = sBuffer;
  }
  else if (totalImageWidth > totalImageHeight) {
    extraTranslateY = (totalImageWidth - totalImageHeight)/2.0;
    extraTranslateX = 0.0;
  }
  else {
    extraTranslateX = (totalImageWidth - totalImageHeight)/2.0;
    extraTranslateY = 0.0;
  }


  // Translate so that all coordinates are positive.  This is unnecessary, but
  // for some reason it helps importing into Illustrator
  TranslateDist[0] += extraTranslateX;
  TranslateDist[1] += extraTranslateY;
  

  // Apply the translations to the drawing elements
  for (i = 0; i < nStrands; i++) {
    arcs[i].x0[0] += TranslateDist[0];
    arcs[i].x0[1] += TranslateDist[1];
    arcs[i].x1[0] += TranslateDist[0];
    arcs[i].x1[1] += TranslateDist[1];
    arrows[i].x0[0] += TranslateDist[0];
    arrows[i].x0[1] += TranslateDist[1];
    arrows[i].x1[0] += TranslateDist[0];
    arrows[i].x1[1] += TranslateDist[1];
    arrows[i].x2[0] += TranslateDist[0];
    arrows[i].x2[1] += TranslateDist[1];
    arrows[i].x3[0] += TranslateDist[0];
    arrows[i].x3[1] += TranslateDist[1];
    arrows[i].x4[0] += TranslateDist[0];
    arrows[i].x4[1] += TranslateDist[1];
    arrows[i].x5[0] += TranslateDist[0];
    arrows[i].x5[1] += TranslateDist[1];
  }
  for (i = 0; i < nbases; i++) {
    dots[i].x[0] += TranslateDist[0];
    dots[i].x[1] += TranslateDist[1];
  }
  for (i = 0; i < nPairs; i++) {
    basepairs[i].x0[0] += TranslateDist[0];
    basepairs[i].x0[1] += TranslateDist[1];
    basepairs[i].x1[0] += TranslateDist[0];
    basepairs[i].x1[1] += TranslateDist[1];
  }


  // Render the SVG
  SVGheader(OutputFile,SVGWidth,SVGHeight,leftCorner,topCorner,vbw,vbh);
  SVGbasepairs(OutputFile,basepairs,nPairs,BPWidth);
  SVGarcs(OutputFile,arcs,arrows,nStrands,ArcWidth);
  SVGdots(OutputFile,dots,bases, nbases,DotRadius,prob,opacity,0);
  if (prob && !opacity) {
    SVGColorBar(OutputFile,dim+COLORBAR_OFFSET_FRAC*dim+sBuffer,sBuffer,dim,white);
  }
  else if (key) { // Show the base coloring key
    SVGBaseKey(OutputFile,DotRadius,dim+BASE_KEY_OFFSET_FRAC*dim+extraTranslateX,
	       dim,thymine,white);
  }
  if (showFreeEnergy) { // print the free energy in lower left corner
    SVGEnergy(OutputFile,freeEnergy,units,freeEnergyPreamble,
	      dim*FREE_ENERGY_SIDE_BUFFER_FRAC + sBuffer,
	      dim*(1.0 + FREE_ENERGY_TOP_BUFFER_FRAC + 
		   FREE_ENERGY_HEIGHT_FRAC*FREE_ENERGY_FONT_SIZE_FRAC)+extraTranslateY,
	      dim,white);
  }
  SVGend(OutputFile);

  free(segLengths);
  free(arcs);
  free(dots);
  free(arrows);
  free(basepairs);

}
/* ******************************************************************************** */

