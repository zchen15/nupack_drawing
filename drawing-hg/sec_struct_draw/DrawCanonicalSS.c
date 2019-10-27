/*
 Renders canonical secondary structure drawing.
   
   See ssdraw.c for comments.
   
   Justin Bois, 29 September, 2007
   */

#include "SecStructDrawHeader.h" // Header file for Concentrations

/* ******************************************************************************** */
int DrawCanonicalSS(struct base *bases, int nbases, double **StandardArrow,
                    double ArrowHeight, double ArrowRadius, int ColorStrands,
                    int AllowOverlap, int white, int prob, int opacity, int drawBases,
                    int showFreeEnergy, double freeEnergy, char *units, 
                    char *freeEnergyPreamble, int key, int thymine, int coords, 
                    int material, int keyBuffer, int colorbarBuffer, 
                    int freeEnergyBuffer, int squareBuffer, int frameBuffer,
                    char *OutputFile, char *coordFile) {
                      
  /*
    Draws a canonical secondary structure.  Returns 1 if successful.
      Returns 0 if the structure cannot be drawn, depending on the
      AllowOverlap specification.
      */
  
  int i,j,h,l; // counters
  int multiTypes; // == 1 if there are multiple strand types in complex
  int npairs; // Total number of base pairs
  int nbackbone; // Number lines in backbone
  int nloops; // Total number of loops
  int nhelices; // Total number of helices
  int nnicks; // Number of nicks (equals number of strands)
  int nsegs; // Number of seqments in loop
  int *looppos; // looppos[i] = 1 if loop is already in position
  double x0[2]; // Place in space
  double u[2]; // Unit vector for helix
  double v[2]; // A vector
  double dv[2]; // Vector along which to translate point in space
  double dist; // Distance for translating
  double **R; // Rotation matrix
  double midpoint[2]; // Midpoint between two paired bases
  double denom; // Normalization we divide by to get unit vectors
  double minlim[2],maxlim[2]; // The extremities of the plot
  double arrowPos[2]; // The position of the arrowhead
  double arrowDir[2]; // The direction the arrowhead points
  double *seglengths; // The lengths of the segments the circle ensribes
  double SVGWidth = SVG_SIZE;  // Dimensions of SVG graphic
  double SVGHeight = SVG_SIZE; // Dimensions of SVG graphic
  double leftCorner, topCorner; // The top left corner coordinates of the view box
  double vbw, vbh; // View box width and height
  double dim; // Square dimension of view box minus colorbar, base key, free energy
  double imageWidth, imageHeight; /// Width and height of sec struc drawing
  double totalImageWidth, totalImageHeight; // Wifth and height of drawing with labels
  double extraTranslateX, extraTranslateY; // Extra distance to translate free energy
  // and colorbar/key.
  double b0; // inter-base length
  double bp0; // base pair distance
  double bs0; // base stack distance
  double n0; // interbase dist for nick
  double b0angle; // Angle in circle corresponding to length b0
  double bp0angle; // Angle in circle corresponding to length bp0
  double n0angle; // Angle in circle corresponding to length n0
  double fb0 = INTERBASE_DISTANCE*ARROW_OVERHANG_FRAC; // overhang distance
  double BackboneWidth = BACKBONE_WIDTH_FRAC*INTERBASE_DISTANCE; // arc line thickness
  double BPWidth = BP_WIDTH_FRAC*INTERBASE_DISTANCE; // bp line thinkness
  double arrowLineWidth = ARROW_LINE_FRAC*INTERBASE_DISTANCE; // arrow line thinkness
  double DotRadius = BASE_RADIUS_FRAC*INTERBASE_DISTANCE; // radius for base dots
  double theta; // Angles along circle plot
  double TranslateDist[2]; // Distance to translate the images.
  double cbBuffer,kBuffer; // Buffer sizes for colorbar and key
  double fBuffer,sBuffer; // Buffer sizes for free energy and square buffering
  FILE *fp; // File for writing coordinates
  struct helix *helices; // Contains information about helices
  struct loop *loops; // Contains information about loops
  struct basedot *dots; // array of dot structs
  struct basepair *basepairs; // array of structs for lines connecting paired bases
  struct backbone *backbones; // array of structs for lines connecting backbone bases
  struct arrowline *arrowlines; // Array of lines from last base to arrow head
  struct arrow *arrows; // Array of arrow structs containing arrowheads
  
  // Special case of a single base  (Two bases is special case also, dealt with below)
  if (nbases <= 1) {
    printf("Error: must have more than one base.  No SVG generated.\n");
    return 1;
  }
  
  
  // Determine physical parameters
  if (coords) {
    if (material == 1) { // RNA
      b0 = INTERBASE_DISTANCE_RNA; // inter-base length
      bp0 = INTERBASE_DISTANCE_RNA*BASE_PAIR_DIST_FRAC_RNA; // base pair distance
      bs0 = INTERBASE_DISTANCE_RNA*STACK_DIST_FRAC_RNA; // base stack distance
      n0 = INTERBASE_DISTANCE_RNA*NICK_PHANTOM_DIST_FRAC; // interbase dist for nick
    }
    else { // DNA
      b0 = INTERBASE_DISTANCE_DNA; // inter-base length
      bp0 = INTERBASE_DISTANCE_DNA*BASE_PAIR_DIST_FRAC_DNA; // base pair distance
      bs0 = INTERBASE_DISTANCE_DNA*STACK_DIST_FRAC_DNA; // base stack distance
      n0 = INTERBASE_DISTANCE_DNA*NICK_PHANTOM_DIST_FRAC; // interbase dist for nick
    }
  }
  else {
    b0 = INTERBASE_DISTANCE; // inter-base length
    bp0 = INTERBASE_DISTANCE*BASE_PAIR_DIST_FRAC; // base pair distance
    bs0 = INTERBASE_DISTANCE*STACK_DIST_FRAC; // base stack distance
    n0 = INTERBASE_DISTANCE*NICK_PHANTOM_DIST_FRAC; // interbase dist for nick
  }
  
  //  printf("b0 = %8.7f\nbp0 = %8.7f\nbs0 = %8.7f\nn0 = %8.7f\n",b0,bp0,bs0,n0);
  
  // Preliminary allocations
  R = (double **) malloc(2 * sizeof(double *));
  R[0] = (double *) malloc(2 * sizeof(double));
  R[1] = (double *) malloc(2 * sizeof(double));
  
  
  // Build the loops and helices structs.  Also update bases with loop information.
  getLoopsHelices(&loops,&helices,&nloops,&nhelices,bases,nbases);
  
  
  // Determine number of phantom bases 
  for (i = 0; i < nloops; i++) {
    if (loops[i].nnicks > 1) {
      printf("Error!  Disconnected structure.  Exiting....\n");
      exit(ERR_DISCONNECTED);
    }
    if (loops[i].nbases == 2) { // blunt helix
      loops[i].nphantom = 0;
    }
    else if (loops[i].nbases == 3) { // One base overhang
      loops[i].nphantom = loops[i].nnicks*N_PHANTOM_1;
    }
    else if (loops[i].nbases == 4) { // Two base overhang
      loops[i].nphantom = loops[i].nnicks*N_PHANTOM_2;
    }
    else {
      loops[i].nphantom = loops[i].nnicks*(N_PHANTOM_BASES_PER_NICK);
    }
  }
  
  // Find the radii of the circles for the loops
  for (i = 0; i < nloops; i++) {
    if (loops[i].nbases == 2) { // Blunt helix
      loops[i].radius = bp0/2.0;
      loops[i].bigAngle = 1;
    }
    else if (loops[i].nbases == 4 && loops[i].nnicks > 0  
             && ((bases[loops[i].bases[0]].pair == loops[i].bases[3]
                  && bases[loops[i].bases[1]].pair == loops[i].bases[2])
                 || (bases[loops[i].bases[0]].pair == loops[i].bases[1]
                     && bases[loops[i].bases[2]].pair == loops[i].bases[3]))) {
                       // Stacked bases in an exterior loop 
                       loops[i].radius = sqrt(bs0*bs0 + bp0*bp0) / 2.0;
                       loops[i].bigAngle = 0;
                     }
    else if (coords && loops[i].nbases <= 0 && loops[i].nnicks == 0) {
      // Hairpin of length 1 or 2 in real coordinates can't stretch
      loops[i].radius = bp0/2.0;
      loops[i].bigAngle = 1;
    }
    else {
      nsegs = loops[i].nbases + loops[i].nphantom;
      seglengths = (double *) malloc(nsegs * sizeof(double));
      for (j = 0; j < loops[i].nhelices; j++) {
        seglengths[j] = bp0;
      }
      for (j = loops[i].nhelices; j < loops[i].nhelices + loops[i].nphantom+1; j++) {
        seglengths[j] = n0;
      }
      for (j = loops[i].nhelices + loops[i].nphantom + 1; j < nsegs; j++) {
        seglengths[j] = b0;
      }
      
      loops[i].radius = CircleRadius(&loops[i].bigAngle,seglengths,nsegs,
                                     TOLERANCE,MAX_ITERATIONS);
      free(seglengths);
    }
  }
  
  
  // Make all helices, just ladders centered origin going in pos x-direction
  for (i = 0; i < nhelices; i++) {
    helices[i].x = (double **) malloc(helices[i].nbases * sizeof(double *));
    for (j = 0; j < helices[i].nbases; j++) {
      helices[i].x[j] = (double *) malloc(2 * sizeof(double));
    }
    helices[i].x[0][0] = -(helices[i].nbases/2-1)*bs0/2.0;
    helices[i].x[0][1] = bp0/2.0;
    helices[i].x[helices[i].nbases-1][0] = helices[i].x[0][0];
    helices[i].x[helices[i].nbases-1][1] = -bp0/2.0;
    for(j = 1; j < helices[i].nbases/2; j++) {
      helices[i].x[j][0] = helices[i].x[j-1][0] + bs0;;
      helices[i].x[j][1] = bp0/2.0;
      helices[i].x[helices[i].nbases-1-j][0] = helices[i].x[j][0];
      helices[i].x[helices[i].nbases-1-j][1] = -bp0/2.0;
    }
  }
  
  
  // Make loops, centered at origin and unrotated.
  for (i = 0; i < nloops; i++) {
    loops[i].x = (double **) malloc(loops[i].nbases * sizeof(double *));
    for (j = 0; j < loops[i].nbases; j++) {
      loops[i].x[j] = (double *) malloc(2 * sizeof(double));
    }
    loops[i].center[0] = 0.0;
    loops[i].center[1] = 0.0;
    if (loops[i].nbases == 2) { // A blunt helix
      loops[i].x[0][0] = -bp0/2.0;
      loops[i].x[0][1] = 0.0;
      loops[i].x[1][0] = bp0/2.0;
      loops[i].x[1][1] = 0.0;
    }
    else if (loops[i].nbases == 4 && loops[i].nnicks > 0  
             && ((bases[loops[i].bases[0]].pair == loops[i].bases[3]
                  && bases[loops[i].bases[1]].pair == loops[i].bases[2])
                 || (bases[loops[i].bases[0]].pair == loops[i].bases[1]
                     && bases[loops[i].bases[2]].pair == loops[i].bases[3]))) {
      // Stacked bases in an exterior loop 
      if (bases[loops[i].bases[0]].pair == loops[i].bases[3]
	  && bases[loops[i].bases[1]].pair == loops[i].bases[2]) { // case 1
	loops[i].x[0][0] = -bp0/2.0;
	loops[i].x[0][1] = -bs0/2.0;
	loops[i].x[1][0] = -bp0/2.0;
	loops[i].x[1][1] = bs0/2.0;
	loops[i].x[2][0] = bp0/2.0;
	loops[i].x[2][1] = bs0/2.0;
	loops[i].x[3][0] = bp0/2.0;
	loops[i].x[3][1] = -bs0/2.0;
      }
      else if (bases[loops[i].bases[0]].pair == loops[i].bases[1]
	       && bases[loops[i].bases[2]].pair == loops[i].bases[3]) { // case 2
	loops[i].x[0][0] = -bp0/2.0;
	loops[i].x[0][1] = bs0/2.0;
	loops[i].x[1][0] = bp0/2.0;
	loops[i].x[1][1] = bs0/2.0;
	loops[i].x[2][0] = bp0/2.0;
	loops[i].x[2][1] = -bs0/2.0;
	loops[i].x[3][0] = -bp0/2.0;
	loops[i].x[3][1] = -bs0/2.0;
      }
    }
    else{
      n0angle = 2.0*asin(n0/(2.0*loops[i].radius));
      b0angle = 2.0*asin(b0/(2.0*loops[i].radius));
      bp0angle = 2.0*asin(bp0/(2.0*loops[i].radius));
      theta = -PI/2.0;
      if (bases[loops[i].bases[0]].pair > -1) { // First base is paired
        theta -= bp0angle/2.0;
      }
      else { // This only happens for loop 0
        theta -= (loops[i].nphantom + 1)*n0angle/2.0;
      }
      for (j = 0; j < loops[i].nbases - 1; j++) {
        loops[i].x[j][0] = loops[i].radius*cos(theta);
        loops[i].x[j][1] = loops[i].radius*sin(theta);
        if (loops[i].bases[j+1] == loops[i].bases[j] + 1) {
          if (bases[loops[i].bases[j]].nextConnection) { // Connected base
            theta -= b0angle;
          }
          else { // nick
            theta -= (loops[i].nphantom + 1)*n0angle;
          }
        }
        else { // base pair
          theta -= bp0angle;
        }
      }
      loops[i].x[loops[i].nbases-1][0] = loops[i].radius*cos(theta);
      loops[i].x[loops[i].nbases-1][1] = loops[i].radius*sin(theta);
    }
  }  
  
  
  // looppos is an array where looppos[i] = 1 if loop i is already in position
  looppos = (int *) malloc(nloops * sizeof(int));
  for (i = 0; i < nloops; i++) {
    looppos[i] = 0;
  }
  
  
  // Position first loop and helix if first base is paired
  if (bases[loops[0].bases[0]].pair > -1) { // First base is paired
    if (loops[0].nbases > 2) { // Not just a blunt helix
      // Rotate loop so first helix points up
      midpoint[0] = (loops[0].x[0][0] + loops[0].x[1][0])/2.0;
      midpoint[1] = (loops[0].x[0][1] + loops[0].x[1][1])/2.0;
      denom = norm(midpoint,2);
      u[0] = midpoint[0]/denom;
      u[1] = midpoint[1]/denom;
      theta = PI/2.0 - atan2(u[1],u[0]);
      getRfromTheta(R,theta);
      for (j = 0; j < loops[0].nbases; j++) {
        v[0] = loops[0].x[j][0];
        v[1] = loops[0].x[j][1];
        MatrixVectorMult(loops[0].x[j],R,v,2);
      }
      
      // Put coordinates in bases
      for (j = 0; j < loops[0].nbases; j++) {
        bases[loops[0].bases[j]].x[0] = loops[0].x[j][0];
        bases[loops[0].bases[j]].x[1] = loops[0].x[j][1];
      }
      
      // Position for first helix
      x0[0] = loops[0].x[0][0];
      x0[1] = loops[0].x[0][1];
    }
    else { // Blunt helix
      // Position for first helix
      x0[0] = -bp0/2.0;
      x0[1] = 0.0;
      
      // Put coordinates in bases
      bases[0].x[0] = -bp0/2.0;
      bases[0].x[1] = 0.0;
      bases[bases[0].pair].x[0] = bp0/2.0;
      bases[bases[0].pair].x[0] = 0.0;
    }
    
    // Specify unit vector
    u[0] = 0.0;
    u[1] = 1.0;
    
    // Rotate helix
    theta = PI/2.0;
    getRfromTheta(R,theta);
    for (j = 0; j < helices[0].nbases; j++) {
      v[0] = helices[0].x[j][0];
      v[1] = helices[0].x[j][1];
      MatrixVectorMult(helices[0].x[j],R,v,2);
    }
    
    // Translate helix
    dv[0] = x0[0] - helices[0].x[0][0];
    dv[1] = x0[1] - helices[0].x[0][1];
    for (j = 0; j < helices[0].nbases; j++) {
      helices[0].x[j][0] += dv[0];
      helices[0].x[j][1] += dv[1];
    }
    
    looppos[0] = 1;
    i = helices[0].corners[2] + 1;
    if (helices[0].corners[0] == helices[0].corners[2]
        && helices[0].corners[2] < nbases - 1) { // lone pair
          l = bases[helices[0].corners[2]+1].loop;
        }
    else {
      l = bases[helices[0].corners[2]].loop;
    }
    h = 1;
  }
  else { // First base is not paired
    // Put coordinates in bases for loop 0
    for (j = 0; j < loops[0].nbases; j++) {
      bases[loops[0].bases[j]].x[0] = loops[0].x[j][0];
      bases[loops[0].bases[j]].x[1] = loops[0].x[j][1];
    }
    looppos[0] = 1;
    i = 0;
    l = 0;
    h = 0;
    theta = 0;
  }
  
  // Go around structure and stick loops and helices together
  while (i < nbases) {
    
    // Position loop if not done already
    if (!looppos[l]) {
      if (loops[l].bigAngle) {
        theta -= PI;
      }
      else {
        theta -= PI/2.0;
      }
      getRfromTheta(R,theta);
      
      // Rotate loop so it aligns along helix
      for (j = 0; j < loops[l].nbases; j++) {
        v[0] = loops[l].x[j][0];
        v[1] = loops[l].x[j][1];
        MatrixVectorMult(loops[l].x[j],R,v,2);
      }
      
      // Translate loop
      midpoint[0] = (helices[h-1].x[helices[h-1].nbases/2-1][0]
                     + helices[h-1].x[helices[h-1].nbases/2][0]) / 2.0;
      midpoint[1] = (helices[h-1].x[helices[h-1].nbases/2-1][1]
                     + helices[h-1].x[helices[h-1].nbases/2][1]) / 2.0;
      if (loops[l].bigAngle) {
        dist = -sqrt(pow(loops[l].radius,2) - pow(bp0,2)/4.0);
      }
      else {
        dist = sqrt(pow(loops[l].radius,2) - pow(bp0,2)/4.0);
      }
      dv[0] = midpoint[0] + u[0]*dist;
      dv[1] = midpoint[1] + u[1]*dist;      
      loops[l].center[0] = dv[0];
      loops[l].center[1] = dv[1];
      for (j = 0; j < loops[l].nbases; j++) {
        loops[l].x[j][0] += dv[0];
        loops[l].x[j][1] += dv[1];
      }
      looppos[l] = 1;
      
      
      // Put coordinates in bases
      for (j = 0; j < loops[l].nbases; j++) {
        bases[loops[l].bases[j]].x[0] = loops[l].x[j][0];
        bases[loops[l].bases[j]].x[1] = loops[l].x[j][1];
      }
    }
    
    // Traverse around loop until we get to a helix or last base
    while (i < nbases && bases[i].pair == -1) {
      i++;
    }
    
    if (i < nbases) { // Not to end yet
      if (i > bases[i].pair) { // Already did this helix, go to its 3' end
        i = helices[bases[i].helix].corners[1] + 1;
        if (i < nbases) {
          l = bases[i].loop;
        }
      }
      else {
        // Get unit vector for helix (outward normal from loop)
        midpoint[0] = (bases[bases[i].pair].x[0] + bases[i].x[0]) / 2.0;
        midpoint[1] = (bases[bases[i].pair].x[1] + bases[i].x[1]) / 2.0;
        u[0] = (midpoint[0]-loops[l].center[0]);
        u[1] = (midpoint[1]-loops[l].center[1]);
        denom = norm(u,2);
        u[0] /= denom;
        u[1] /= denom;
        theta = atan2(u[1],u[0]);
        getRfromTheta(R,theta);
        
        // Rotate helix so it aligns along unit vector
        for (j = 0; j < helices[h].nbases; j++) {
          v[0] = helices[h].x[j][0];
          v[1] = helices[h].x[j][1];
          MatrixVectorMult(helices[h].x[j],R,v,2);
        }
        
        // Translate helix
        dv[0] = bases[i].x[0] - helices[h].x[0][0];
        dv[1] = bases[i].x[1] - helices[h].x[0][1];
        for (j = 0; j < helices[h].nbases; j++) {
          helices[h].x[j][0] += dv[0];
          helices[h].x[j][1] += dv[1];
        }
        
        i = helices[h].corners[2] + 1;
        if (helices[h].corners[0] == helices[h].corners[2]
            && helices[h].corners[2] < nbases - 1) { // lone pair
              l = bases[helices[h].corners[2]+1].loop;
            }
        else {
          l = bases[helices[h].corners[2]].loop;
        }
        h++;
      }
    }      
  } 
  
  
  // Get lines for base pairs
  npairs = 0;
  for (i = 0; i < nhelices; i++) {
    npairs += helices[i].nbases/2;
  }
  
  dots = (struct basedot *) malloc(nbases * sizeof(struct basedot));
  // Transfer coordinates from helices to bases
  for (i = 0; i < nhelices; i++) {
    for (j = 0; j < helices[i].nbases; j++) {
      bases[helices[i].bases[j]].x[0] = helices[i].x[j][0];
      bases[helices[i].bases[j]].x[1] = helices[i].x[j][1];
    }
  }
  
  
  // Write coordinate file if necessary
  if (coords) {
    if ((fp = fopen(coordFile,"w")) == NULL) {
      printf("Error opening output file %s\n",OutputFile);
      printf("\nExiting....\n\n");
      exit(ERR_OUTPUTFILE);
    }
    
    for (i = 0; i < nbases; i++) {
      fprintf(fp,"%d\t%8.7f\t%8.7f\n",i+1,bases[i].x[0],bases[i].x[1]);
    }
    
    fclose(fp);
    
    // Free memory and return 0    
    for (i = 0; i < nloops; i++) {
      free(loops[i].bases);
      for (j = 0; j < loops[i].nbases; j++) {
        free(loops[i].x[j]);
      }
      free(loops[i].x);
    }
    free(loops);
    return 0;
  }
  
  
  // Check for clashes
  if (!AllowOverlap && detectClashes(bases,nbases)) {
    // For a clash, free memory and return 0    
    for (i = 0; i < nloops; i++) {
      free(loops[i].bases);
      for (j = 0; j < loops[i].nbases; j++) {
        free(loops[i].x[j]);
      }
      free(loops[i].x);
    }
    free(loops);
    return 0;
  }
  
  
  // Get the arrows
  nnicks = 0;
  for (i = 0; i < nloops; i++) {
    nnicks += loops[i].nnicks;
  }
  
  
  arrows = (struct arrow *) malloc(nnicks * sizeof(struct arrow));
  arrowlines = (struct arrowline *) malloc(nnicks * sizeof(arrowline));
  
  
  j = 0;
  for (i = 0; i < nbases; i++) {
    if (bases[i].nextConnection == 0) { // At 3' end of strand
      if (bases[i].pair > -1) { // 3'-most base is paired
        // Perpendicular to base pair
        arrowDir[0] = bases[i].x[1] - bases[bases[i].pair].x[1];
        arrowDir[1] = bases[bases[i].pair].x[0] - bases[i].x[0];
        denom = sqrt(dot(arrowDir,arrowDir,2));
        arrowDir[0] /= denom;
        arrowDir[1] /= denom;
      }
      else { // 3'-most base is unpaired
        // Arrowhead, direction is the negative tangent to the circle
        theta = atan2(bases[i].x[1] - loops[bases[i].loop].center[1],
                      bases[i].x[0] - loops[bases[i].loop].center[0]);
        arrowDir[0] = sin(theta);
        arrowDir[1] = -cos(theta);
      }
      arrowlines[j].x0[0] = bases[i].x[0];
      arrowlines[j].x1[0] = bases[i].x[0] + fb0*arrowDir[0];
      arrowlines[j].x0[1] = bases[i].x[1];
      arrowlines[j].x1[1] = bases[i].x[1] + fb0*arrowDir[1];
      arrowPos[0] = arrowlines[j].x1[0];
      arrowPos[1] = arrowlines[j].x1[1];
      getArrow(arrows,j,StandardArrow,arrowPos,arrowDir);
      j++;
    }
  }
  
  
  // Reflect all coordinates through x-axis
  for (i = 0; i < nbases; i++) {
    bases[i].x[1] *= -1;
  }
  for (i = 0; i < nnicks; i++) {
    arrows[i].x0[1] *= -1;
    arrows[i].x1[1] *= -1;
    arrows[i].x2[1] *= -1;
    arrows[i].x3[1] *= -1;
    arrows[i].x4[1] *= -1;
    arrows[i].x5[1] *= -1;
    arrowlines[i].x0[1] *= -1;
    arrowlines[i].x1[1] *= -1;
  }
  
  
  basepairs = (struct basepair *) malloc(npairs * sizeof(struct basepair));
  j = 0;
  for (i = 0; i < nbases; i++) {
    if (i < bases[i].pair) {
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
  
  
  // Get lines connecting bases (backbone)
  nbackbone = 0;
  for (i = 0; i < nbases - 1; i++) {
    if (bases[i].nextConnection) {
      nbackbone++;
    }
  }
  backbones = (struct backbone *) malloc(nbackbone * sizeof(struct backbone));
  j = 0;
  for (i = 0; i < nbases - 1; i++) {
    if (bases[i].nextConnection) {
      backbones[j].x0[0] = bases[i].x[0];
      backbones[j].x0[1] = bases[i].x[1];
      backbones[j].x1[0] = bases[i+1].x[0];
      backbones[j].x1[1] = bases[i+1].x[1];
      if (bases[i].pair > -1 && bases[i+1].pair > -1  
          && bases[i+1].pair == bases[i].pair - 1 
          && bases[bases[i].pair].prevConnection
          && bases[i].pair != i+1) {
            backbones[j].curved = 0;
          }
      else {
        backbones[j].curved = 1;
        if (i+1 == nbases-1 || 
            (bases[i+1].pair > -1 && helices[bases[i+1].helix].nbases == 2)) {
              backbones[j].radius = loops[bases[i].loop].radius;
            }
        else { 
          backbones[j].radius = loops[bases[i+1].loop].radius;
        }
      }
      j++;
    }
  }
  
  
  // Free memory we don't need any more
  free(R[0]);
  free(R[1]);
  free(R);
  for (i = 0; i < nhelices; i++) {
    free(helices[i].bases);
    for (j = 0; j < helices[i].nbases; j++) {
      free(helices[i].x[j]);
    }
    free(helices[i].x);
  }
  free(helices);
  free(looppos);
  
  
  // Put bases in dot structs for writing
  dots = (struct basedot *) malloc(nbases * sizeof(struct basedot));
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
  
  
  // Color the backbones
  // Check to see if we have multiple strand types
  multiTypes = 0;
  for (i = 1; i < nbases && multiTypes == 0; i++) {
    if (bases[i].strandID != bases[0].strandID) {
      multiTypes = 1;
    }
  }
  
  // Right now, we only use default strand color
  //  if (!ColorStrands || multiTypes == 0) { // Make them all default
  if (white) {
    for (i = 0; i < nbackbone; i++) {
      backbones[i].RGB[0] =  StrandColorWhite[0];
      backbones[i].RGB[1] =  StrandColorWhite[1];
      backbones[i].RGB[2] =  StrandColorWhite[2];
    }
    for (j = 0; j < nnicks; j++) {
      arrows[j].RGB[0] =  StrandColorWhite[0];
      arrows[j].RGB[1] =  StrandColorWhite[1];
      arrows[j].RGB[2] =  StrandColorWhite[2];
      arrowlines[j].RGB[0] =  StrandColorWhite[0];
      arrowlines[j].RGB[1] =  StrandColorWhite[1];
      arrowlines[j].RGB[2] =  StrandColorWhite[2];
    }
  }
  else{
    for (i = 0; i < nbackbone; i++) {
      backbones[i].RGB[0] =  DefaultStrandColor[0];
      backbones[i].RGB[1] =  DefaultStrandColor[1];
      backbones[i].RGB[2] =  DefaultStrandColor[2];
    }
    for (j = 0; j < nnicks; j++) {
      arrows[j].RGB[0] =  DefaultStrandColor[0];
      arrows[j].RGB[1] =  DefaultStrandColor[1];
      arrows[j].RGB[2] =  DefaultStrandColor[2];
      arrowlines[j].RGB[0] =  DefaultStrandColor[0];
      arrowlines[j].RGB[1] =  DefaultStrandColor[1];
      arrowlines[j].RGB[2] =  DefaultStrandColor[2];
    }
  }
  
  /*
    else {
      for (i = 0; i < nbackbone; i++) {
        backbones[i].RGB[0] =  StrandColors[bases[i].strandID][0];
        backbones[i].RGB[1] =  StrandColors[bases[i].strandID][1];
        backbones[i].RGB[2] =  StrandColors[bases[i].strandID][2];
      }
      for (j = 0; j < nnicks; j++) {
        arrows[j].RGB[0] =  StrandColors[bases[i].strandID][0];
        arrows[j].RGB[1] =  StrandColors[bases[i].strandID][1];
        arrows[j].RGB[2] =  StrandColors[bases[i].strandID][2];
        arrowlinesj].RGB[0] =  StrandColors[bases[i].strandID][0];
        arrowlines[j].RGB[1] =  StrandColors[bases[i].strandID][1];
        arrowlines[j].RGB[2] =  StrandColors[bases[i].strandID][2];
        }
        }
        */
  
  
  // Find the limits of the graphic
  minlim[0] = dots[0].x[0];
  minlim[1] = dots[0].x[1];
  maxlim[0] = dots[0].x[0];
  maxlim[1] = dots[0].x[1];
  
  for (i = 1; i < nbases; i++) {
    if (dots[i].x[0] - DotRadius < minlim[0]) {
      minlim[0] = dots[i].x[0] - DotRadius;
    }
    if (dots[i].x[0] + DotRadius > maxlim[0]) {
      maxlim[0] = dots[i].x[0] + DotRadius;
    }
    if (dots[i].x[1] - DotRadius < minlim[1]) {
      minlim[1] = dots[i].x[1] - DotRadius;
    }
    if (dots[i].x[1] + DotRadius > maxlim[1]) {
      maxlim[1] = dots[i].x[1] + DotRadius;
    }
  }
  for (i = 0; i < nnicks; i++) {
    if (arrows[i].x2[0] < minlim[0]) {
      minlim[0] = arrows[i].x2[0];
    }
    if (arrows[i].x2[0] > maxlim[0]) {
      maxlim[0] = arrows[i].x2[0];
    }
    if (arrows[i].x2[1] < minlim[1]) {
      minlim[1] = arrows[i].x2[1];
    }
    if (arrows[i].x2[1] > maxlim[1]) {
      maxlim[1] = arrows[i].x2[1];
    }
    if (arrows[i].x3[0] < minlim[0]) {
      minlim[0] = arrows[i].x3[0];
    }
    if (arrows[i].x3[0] > maxlim[0]) {
      maxlim[0] = arrows[i].x3[0];
    }
    if (arrows[i].x3[1] < minlim[1]) {
      minlim[1] = arrows[i].x3[1];
    }
    if (arrows[i].x3[1] > maxlim[1]) {
      maxlim[1] = arrows[i].x3[1];
    }
    if (arrows[i].x4[0] < minlim[0]) {
      minlim[0] = arrows[i].x4[0];
    }
    if (arrows[i].x4[0] > maxlim[0]) {
      maxlim[0] = arrows[i].x4[0];
    }
    if (arrows[i].x4[1] < minlim[1]) {
      minlim[1] = arrows[i].x4[1];
    }
    if (arrows[i].x4[1] > maxlim[1]) {
      maxlim[1] = arrows[i].x4[1];
    }
  }
  
  
  // Calculate width and height of image
  imageWidth = 2*EDGE_BUFFER_FRAC*(maxlim[0] - minlim[0]) + 2*DotRadius 
    + maxlim[0] - minlim[0];
  imageHeight = 2*EDGE_BUFFER_FRAC*(maxlim[1] - minlim[1]) + 2*DotRadius 
    + maxlim[1] - minlim[1];

  // Special case of two bases, adjust image height so we don't cut off the arc
  if (nbases == 2) {
    imageHeight -= 1.0*DotRadius;
  }
  
  // Specifications for SVG
  // We make the image square for display on the web site.
  // Right now, this is the only option.  We might wish to change this in the future
  // with a -square or a -website flag.  Once this is implemented, the -framebuffer
  // flag can also be activated.
  totalImageWidth = max2(imageWidth,imageHeight);
  totalImageHeight = totalImageWidth; // The image dimension is square
  dim = totalImageHeight; // Store the square dimension of the image in the view box
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
    totalImageHeight = max2(imageHeight,imageWidth) + 2.0*sBuffer;
    totalImageWidth = totalImageHeight;
  }
  else { // Set the square buffer to zero otherwise
    sBuffer = 0.0;
  }
  vbh = max2(totalImageHeight,totalImageWidth);
  vbw = vbh;
  SVGHeight = vbh;
  SVGWidth = vbw;
  
  // These translation distances of the images assume a square image.  This may change
  // in the future.
  if (imageWidth > imageHeight) {
    TranslateDist[0] = 0.0;
    TranslateDist[1] = (imageWidth-imageHeight)/2.0;
  }
  else {
    TranslateDist[0] = (imageHeight-imageWidth)/2.0;
    TranslateDist[1] = 0.0;
  }
  
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
  TranslateDist[0] += EDGE_BUFFER_FRAC*(maxlim[0] - minlim[0]) + DotRadius
    - minlim[0] + extraTranslateX;
  TranslateDist[1] += EDGE_BUFFER_FRAC*(maxlim[1] - minlim[1]) + DotRadius
    - minlim[1] + extraTranslateY;
  
  for (i = 0; i < nbases; i++) {
    dots[i].x[0] += TranslateDist[0];
    dots[i].x[1] += TranslateDist[1];
  }
  for (i = 0; i < npairs; i++) {
    basepairs[i].x0[0] += TranslateDist[0];
    basepairs[i].x0[1] += TranslateDist[1];
    basepairs[i].x1[0] += TranslateDist[0];
    basepairs[i].x1[1] += TranslateDist[1];
  }
  for (i = 0; i < nbackbone; i++) {
    backbones[i].x0[0] += TranslateDist[0];
    backbones[i].x0[1] += TranslateDist[1];
    backbones[i].x1[0] += TranslateDist[0];
    backbones[i].x1[1] += TranslateDist[1];
  }
  for (i = 0; i < nnicks; i++) {
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
    arrowlines[i].x0[0] += TranslateDist[0];
    arrowlines[i].x0[1] += TranslateDist[1];
    arrowlines[i].x1[0] += TranslateDist[0];
    arrowlines[i].x1[1] += TranslateDist[1];
  }
  
  
  // Render the SVG
  SVGheader(OutputFile,SVGWidth,SVGHeight,leftCorner,topCorner,vbw,vbh);
  SVGarrows(OutputFile,arrows,arrowlines,nnicks,arrowLineWidth);
  SVGbasepairs(OutputFile,basepairs,npairs,BPWidth);
  SVGbackbone(OutputFile,backbones,nbackbone,BackboneWidth);
  SVGdots(OutputFile,dots,bases,nbases,DotRadius,prob,opacity,drawBases);
  if (prob && !opacity) {
    SVGColorBar(OutputFile,dim+COLORBAR_OFFSET_FRAC*dim+extraTranslateX,
                extraTranslateY,dim,white);
  }
  else if (key) { // Show the base coloring key
    SVGBaseKey(OutputFile,DotRadius,dim+BASE_KEY_OFFSET_FRAC*dim+extraTranslateX,
               dim,thymine,white);
  }
  if (showFreeEnergy) { // print the free energy in lower left corner
    SVGEnergy(OutputFile,freeEnergy,units,freeEnergyPreamble,
              dim*FREE_ENERGY_SIDE_BUFFER_FRAC + extraTranslateX,
              dim*(1.0 + FREE_ENERGY_TOP_BUFFER_FRAC + 
                   FREE_ENERGY_HEIGHT_FRAC*FREE_ENERGY_FONT_SIZE_FRAC)+extraTranslateY,
              dim,white);
  }
  SVGend(OutputFile);
  
  
  // Free memory
  for (i = 0; i < nloops; i++) {
    free(loops[i].bases);
    for (j = 0; j < loops[i].nbases; j++) {
      free(loops[i].x[j]);
    }
    free(loops[i].x);
  }
  free(loops);
  
  free(dots);
  free(basepairs);
  free(backbones);
  free(arrows);
  free(arrowlines);
  
  return 1; // No overlaps
  
}

/* ******************************************************************************** */


/* ******************************************************************************** */
void getLoopsHelices(struct loop **loops, struct helix **helices, int *nloops, 
int *nhelices, struct base *bases, int nbases) {
  /*
   Builds the loops and helices structs from the input data and
     updates the loop information for the bases.
     
     helices and loops are allocated here, but are freed elsewhere.
     */

  int i,j,k; // Counters
  int *helixcount; // An array of counts of bases in each helix 

  // Allocate plenty of memory for helix count
  helixcount = (int *) malloc(nbases * sizeof(int));
  for (i = 0; i < nbases; i++) {
   helixcount[i] = 0;
  }


  // Get the helix labels
  if (bases[0].pair > -1) {
   bases[0].helix = 0;
   (*nhelices) = 0;
   (helixcount[bases[0].helix])++;
  }
  else {
   bases[0].helix = -1;
   (*nhelices) = -1;
  }
  for (i = 1; i < nbases; i++) {
   if (bases[i].pair == -1) {
     bases[i].helix = -1;
   }
   else if (bases[i].pair < i) {
     bases[i].helix = bases[bases[i].pair].helix;
     (helixcount[bases[i].helix])++;
   }
   else if (bases[i-1].pair == -1 // prev. base unpaired; new helix 
            || bases[i-1].pair != bases[i].pair + 1 // prev. base paired, bulge loop
            || bases[i].prevConnection == 0 // strand break on this side of helix
            || bases[bases[i].pair].nextConnection == 0) { //str. break on other side
              (*nhelices)++;
              bases[i].helix = (*nhelices);
              (helixcount[bases[i].helix])++;
            }
   else {
     bases[i].helix = bases[i-1].helix;
     (helixcount[bases[i].helix])++;
   }
  }
  (*nhelices)++;


  // Allocate memory for struct for helices
  (*helices) = (struct helix *) malloc((*nhelices) * sizeof(struct helix));
  for (i = 0; i < (*nhelices); i++) {
   (*helices)[i].bases = (int *) malloc(helixcount[i] * sizeof(int));
   (*helices)[i].nbases = helixcount[i];
  }


  // Put in information to helix struct
  i = 0;
  k = -1;
  for (i = 0; i < nbases; i++) {
   if (bases[i].helix > k) {
     k++;
     for (j = 0; j < (*helices)[bases[i].helix].nbases/2; j++) {
       (*helices)[k].bases[j] = i + j;
       (*helices)[k].bases[(*helices)[k].nbases-1-j] = bases[i+j].pair;
     }
     (*helices)[k].corners[0] = (*helices)[k].bases[0];
     (*helices)[k].corners[1] = (*helices)[k].bases[(*helices)[k].nbases-1];
     (*helices)[k].corners[2] = (*helices)[k].bases[(*helices)[k].nbases/2 - 1];
     (*helices)[k].corners[3] = (*helices)[k].bases[(*helices)[k].nbases/2];
   }
  }


  // Allocate memory for loops
  (*nloops) = (*nhelices) + 1;
  (*loops) = (struct loop *) malloc((*nloops) * sizeof(struct loop));

  // Get counts of number of bases in loops
  // First for the loop containing strand ends
  i = 0;
  (*loops)[0].nbases = 1;
  while (i != nbases-1) {
   if (bases[i].pair > i) {
     i = bases[i].pair;
   }
   else {
     i++;
   }
   ((*loops)[0].nbases)++;
  }


  // Now all other loops; they begin with 3' end of a helix
  for (j = 0; j < (*nhelices); j++) {
   i = (*helices)[j].corners[2] + 1;
   (*loops)[j+1].nbases = 1;
   while (i != bases[(*helices)[j].corners[2]].pair) {
     if (bases[i].pair > i) {
       i = bases[i].pair;
     }
     else {
       i++;
     }
     ((*loops)[j+1].nbases)++;
   }
   ((*loops)[j+1].nbases)++; // Have one more from ending base
  }


  // Fill in bases in loops
  for (j = 0; j < (*nloops); j++) {
   (*loops)[j].bases = (int *) malloc((*loops)[j].nbases * sizeof(int));
  }
  // First for the loop containing strand ends
  i = 0;
  k = 0;
  (*loops)[0].nnicks = 0;
  (*loops)[0].nhelices = 0;
  (*loops)[0].bases[k++] = i;
  while (i != nbases-1) {
   if (bases[i].prevConnection == 0 && bases[i].pair) {
     ((*loops)[0].nnicks)++;
   }
   if (bases[i].pair > i) {
     i = bases[i].pair;
     ((*loops)[0].nhelices)++;
   }
   else {
     i++;
   }
   (*loops)[0].bases[k++] = i;
   if ((*loops)[0].nbases == 2) { // Must have 1 nick and we never entered while loop
     (*loops)[0].nnicks = 1;
   }
  }


  // Now all other loops; they begin with 3' end of a helix
  for (j = 0; j < (*nhelices); j++) {
   i = (*helices)[j].corners[2] + 1;
   (*loops)[j+1].bases[0] = (*helices)[j].corners[2]; // First base
   k = 1;
   (*loops)[j+1].nnicks = 0;
   (*loops)[j+1].nhelices = 1; // Automatically has one from first helix
   (*loops)[j+1].bases[k++] = i;
   while (i != bases[(*helices)[j].corners[2]].pair) {
     if (bases[i].prevConnection == 0 && bases[i].pair != (*loops)[j+1].bases[k-2]) {
       ((*loops)[j+1].nnicks)++;
     }
     if (bases[i].pair > i) {
       i = bases[i].pair;
       ((*loops)[j+1].nhelices)++;
     }
     else {
       i++;
     }
     (*loops)[j+1].bases[k++] = i;
   }
   
   // Check to see if there's a nick right before last base in loop
   if (bases[bases[(*helices)[j].corners[2]].pair].prevConnection == 0) {
     ((*loops)[j+1].nnicks)++;
   }
   
   if ((*loops)[j+1].nbases == 2) { // Must have 1 nick and never entered while loop
     (*loops)[j+1].nnicks = 1;
   }
  }


  // Update bases with loop information (lone pair base belongs to 5' most loop
  // That's why we traverse loop backwards
  for (i = 0; i < nbases; i++) {
   bases[i].loop = -1;
  }
  for (i = (*nloops)-1; i >= 0; i--) {
   for (j = 0; j < (*loops)[i].nbases; j++) {
     bases[(*loops)[i].bases[j]].loop = i;
   }
  }

  free(helixcount);

}

/* ******************************************************************************** */


/* ******************************************************************************** */
int detectClashes(struct base *bases, int nbases) {
  /*
    Quick and dirty way of clash detection.  Loops through all dots,
    lines, and arrow and finds clashes.  Returns 1 is there is a clash
      and 0 otherwise.
      */
  
  int i,j; // Counters
  int nlines; // Number of lines in drawing.
  int clash = 0; // Whether or not there is a clash.
  double buffer = INTERBASE_DISTANCE*CLASH_BUFFER_FRAC; // clash buffer
  double **lines; // See description below.
  int **connLines; // See descriptions below.
  
  
  // First check base-base clashes, since it's easiest
  for (i = 0; i < nbases-1; i++) {
    for (j = i+1; j < nbases; j++) {
      if (pointpointdist(bases[i].x,bases[j].x,2) < buffer) { // Clash!
        return 1;
      }
    }
  }
  
  
  // Allocate plenty of memory for lines and connLines
  // The maximum number of lines is nbases*2.
  lines = (double **) malloc(nbases*2 * sizeof(double *));
  connLines = (int **) malloc(nbases*2 * sizeof(int *));
  
  // Go around structure and build array of lines
  // lines[i][0] and [i][1] is the x,y coords of first point in line i
  // lines[i][2] and [i][3] is the x,y coords of the second point in line i
  // connLines[i][0] and connLines[i][1] are the two bases line i connects
  nlines = 0;
  for (i = 0; i < nbases; i++) {
    if (bases[i].nextConnection) { // Backbone line
      lines[nlines] = (double *) malloc(4 * sizeof(double));
      connLines[nlines] = (int *) malloc(2 * sizeof(int));
      lines[nlines][0] = bases[i].x[0];
      lines[nlines][1] = bases[i].x[1];
      lines[nlines][2] = bases[i+1].x[0];
      lines[nlines][3] = bases[i+1].x[1];
      connLines[nlines][0] = i;
      connLines[nlines][1] = i+1;
      nlines++;
    }
    if (bases[i].pair > i) { // basepair line
      lines[nlines] = (double *) malloc(4 * sizeof(double));
      connLines[nlines] = (int *) malloc(2 * sizeof(int));
      lines[nlines][0] = bases[i].x[0];
      lines[nlines][1] = bases[i].x[1];
      lines[nlines][2] = bases[bases[i].pair].x[0];
      lines[nlines][3] = bases[bases[i].pair].x[1];
      connLines[nlines][0] = i;
      connLines[nlines][1] = bases[i].pair;
      nlines++;
    }
  }
  
  
  // Check point-line clashes
  clash = 0;
  for (i = 0; i < nbases && clash == 0; i++) {
    for (j = 0; j < nlines && clash == 0; j++) {
      if (i != connLines[j][0] && i != connLines[j][1]) {
        if (pointlinedist(bases[i].x,lines[j],2) < buffer ) {
          clash = 1;
        }
      }
    }
  }
  
  
  // Check line crossings
  for (i = 0; i < nlines-1 && clash == 0; i++) {
    for (j = i+1; j < nlines && clash == 0; j++) {
      // Only check if they do not share a point
      if (connLines[i][0] != connLines[j][0] && connLines[i][1] != connLines[j][0] &&
          connLines[i][0] != connLines[j][1] && connLines[i][1] != connLines[j][1]) {
            if (linecross2d(lines[i],lines[j])) {
              clash = 1;
            }
          }
    }
  }
  
  for (i = 0; i < nlines; i++) {
    free(lines[i]);
    free(connLines[i]);
  }
  free(lines);
  free(connLines);
  
  return clash;
  
}
/* ******************************************************************************** */
