#! /usr/bin/python

# Python script for generating dot plots from outputs of complexes and
# concentrations for use with NUPACK.
#
# In running this code, it is called from the command line as follows:
# DotPlot.py [-mfe] -[target] [-complex ComplexID] [-order PermID] [-grid]
# [-svgfile OutputFile] [-title TitleString] [-labels LabelFile]
# [-datafile dataOutputFile] [-longheader] [-energy free_energy]
# [-freeenergy] [-energypreamble preamble_text] -[units units_text] prefix.suffix
#
# -mfe is chosen if the lower triangle is to show the MFE structure.
# -target is chosen if the target structure is overlayed.  The file name
#    for the target is prefix.target or prefix.ocx-target.  The format
#    of the .target file is identical to that of the .mfe files in NUPACK
#    (without the line containing the value of the free energy) and that
#    of the .ocx-target file is identical to that of the .ocx-mfe files
#    in NUPACK (also without the line containing the value of the free
#    energy).   NOTE: PSEUDOKNOTS ARE NOT ALLOWED FOR DOT-PAREN STRUCTURES!
#    IF THE DOT-PAREN STRUCTURE HAS PSEUDOKNOTS, THERE WILL BE AN ERROR.
# -mfetarget
#  Same as -target, except the mfe file is used as a target
# -complex flag gives the ID of the complex we're plotting for
# -order flag gives the permutation ID of the ordered complex
# -grid flag indicates that grid lines are drawn in plot
# -title flag gives the title for the plot.  If none is chosen,
#    the default is not to have a title.
# -svgfile gives the name of the file to write the output SVG.  If none is
#    chosen, the default is prefix-dotplot.svg.
# -labels gives the name of the file containing the strand labels.  Line
#    1 is the strand label for the first strand, the second line for the
#    second, and so on.  Number of entries must equal number of strands.
#    Right now, there is no error checking for this.
# -datafile gives the name of the file to which the data used to
#    generate the dot plot are written.
# -longheader is chosen for a verbose header.
# -energy: Argument of flag is printed on plot as free energy in %.2f format.
#          Overrides the -freeenergy flag
# -freeenergy: No argument, active flag means free energy is extracted from input file
#              and printed on dot plot in %.2f format.
# -energypreamble: Argument of flag is printed before the free energy, followed
#                  immediately by a single whitespace and then the free energy.
#                  Example: -energypreamble "Free energy (log Q):"
# -units: A string giving units of free energy, default is kcal/mol.
# -freeenergybuffer: Makes whitespace on the drawing of the same size the free
#                    energy text would be on the bottom of the drawing
# -fontscale NUMBER  Scale fonts from the default by this factor
#
# prefix.suffix is the name of the input file.  Acceptable values for suffix
# are:
#
# .ocx-ppairs (-complex and -order flags required)
# .ocx-epairs (-complex and -order flags required, -mfe flag ignored))
# .cx-epairs (-complex flag required, -mfe and -order flag ignored)
# .fpairs (-complex, -order, -mfe, -energy, -freeeenergy flags all ignored)
# .ppairs (-complex, -order, -mfe flags all ignored)
#
# The file prefix.cx must always be present.
# prefix.ocx-key must be present for suffix = .ocx-ppairs
# prefix.ocx-mfe must be present when the -mfe flag is chosen and
#   suffix = .ocx-ppairs.
#
# For now, we assume no sorting of the entries in the output files
# from complexes and concentrations.  We therefore have an inefficient
# way of finding the permutation of interest.  If sorting is assumed,
# there are more efficient ways to pick out the permutation of
# interest.  The plots are still generated quickly despite this.
#
# There isn't much error checking here.  All input from other NUPACK
# code is assumed to be consistent and correct.
#
#
# Right now we don't translate dot-paren MFE (or design) structures to
# fit into the dot plot framework for complexes where the strands are
# treated as distinguishable. Therefore, the only option for plotting
# an MFE in the lower triangle is for a specific ordered complex when
# the suffix is ppairs.
#
# Justin Bois, bois AT caltech DOT edu
# Caltech
# 25 October, 2006
#
# Updated 24 November, 2006
#
# Updated 24 December, 2006
#
# Updated 29 September, 2007
#
# Updated 18 October, 2007


# Import a modules we'll use
import sys
argv = sys.argv
exit = sys.exit
import re
WhiteSpaceSearch = re.compile('\s+')
import math
log10 = math.log10
floor = math.floor
ceil = math.ceil
sqrt = math.sqrt
import time
strftime = time.strftime
import ColorMap
getRGB = ColorMap.getRGB
import getInputData
getInputDataComplexes = getInputData.getInputDataComplexes
getInputDataNuPackLib = getInputData.getInputDataNuPackLib
getInputDataTarget = getInputData.getInputDataTarget
import writeDataFile
writeDataFile = writeDataFile.writeDataFile


# ######## CONSTANTS DESCRIBING DIMENSIONS, ETC OF PLOT ######
Font = 'sans-serif'

MakeSquare = True # Make the overall SVG square
IllustratorCompliant = 1 # Make Illlustrator compliant (overrides SVGHeight and Width)

DiscreteScaling = 0 # Color and box size scaled discrete, specified by minP, below
BoxScaling = 1 # Whether or not to scale the size of the boxes
AreaScaling = 1  # = 1 is scaling of box size is by area and 0 if by side length

red,green,blue = getRGB(1.0)
OtherColorContinuous = [0,0,0] # Color of squares of structure in lower triangle (black)

# Options for discrete scaling
MediumBigRatio = 0.5 # Ratio of area of medium square to big square
SmallBigRatio = 0.25 # Ratio of area of small square to big square
BigColor = [255,0,0] # RGB (red)
MediumColor = [0,100,0] # RGB (green)
SmallColor = [0,0,255] # RGB (blue)
OtherColorDiscrete = [0,0,0] # Color of squares of structure in lower triangle (black)
TargetColor = [0,0,0] # Color for target structure (black)
TargetLineWidthFrac = 0.15  # Width of line in targets as fraction of full base width
#BaseOpacityOverTarget = 0.5 # Opacity of base over target
noX = 1

# Cutoffs for size of boxes (taken from paradigms paper)
# This is used if DiscreteScaling = 1.  MUST HAVE ONLY 3 ELEMENTS
# If have continuous scaling, the last entry is the cutoff probability
# for being shown on the plot
minP = [0.75,0.5,0.005]

# Font scale factor
FontScale = 1.0

# Find the fontscale option
if '-fontscale' in argv:
    FontScale = float(argv[argv.index('-fontscale') + 1])



# Options for the colorbar
ColorBarFrac = 1.0  # Fraction of total height of plot taken b colorbar
ColorBarWidth = 20   # Width of the colorbar
NColorBarLabels = 6 # Number of colorbar labels (6 give [0.0 , 0.2 , 0.4 , 0.6 , 0.8 , 1.0]
ColorBarLabelOffset = 6  # Distance from colorbar of label
ColorBarLabelFontSize = 24*FontScale  # Font size of colorbar labels
ColorBarLabelColor = 'black'  # Color of test for colorbar labels
NColors = 256 # Number of levels in RGB graphics

# Works well when there is a colormap
WidthToHeightRatio = 1.06  # old value is 1.17
SVGHeight = 500 # Height in pixels.  This works for all plots.
vbh = 1045 # height of view box (this works well)
# vbw and SVGWidth are set off the height and such that everything fits
d = 800.0 # Width and height of dot plot matrix
TopMarginFull = 145.0 # Distance from top of document to top of dot matrix
TopMarginFrac00 = 1.0/12.0 # Fraction of margin used when there is no title and no strand labels
TopMarginFrac10 = 2.0/3.0 # Fraction of margin used when there is a title but no strand labels
TopMarginFrac01 = 3.0/8.0 # Fraction of margin used when there is no title but there are strand labels

LeftMargin = 102.5 # Distance from left of document to left edge of dot matrix (excludes axis labels)
RightMargin = 50.0 # Distance from edge of colorbar to edge of document (excludes colorbar labels)

AxisWidth = 2.0 # Width of lines of axes
TickWidth = 1.0 # Width of tick lines
GridWidth = 0.75 # Width of grid lines
BreakWidth = 2.0 # Width of line marking strand breaks
DiagWidth = 2.0 # Width of the diagonal line in the dot plot
AxisColor = 'gray' # Color of axes
TickColor = 'gray' # Color of ticks
GridColor = 'gray' # Color of grid lines
BreakColor = 'gray' # Color of break lines
DiagColor = 'gray' # Color of diagonal line
GridDashLength = 5.0 # Length of a dash in gridlines
GridGapLength = 5.0 # Length of a gap in the dash array for grid lines
TickLengthFrac = 0.015 # Fraction of the length of axes that ticks are
StrandBreakGap = 1 # Strand break gap in units of size of a full square

LastTick = 1  # Whether or not to include a tick at the strand length if withing LastTickFrac
LastTickFrac = 0.5  # The maximal fraction toward the end of the strand that the last tick will be
AdjustStrandLabelFontSize = 1 # Adjust strand label font size instead of truncating strand label
AllStrandLabelsSameSize = 1 # Works well for web with cap on strand label text
StrandLabelFontSizeStandard = 36*FontScale  # Font size of the strand labels
AxisLabelFontSize = 36*FontScale  # Font size of the axis labels
AxisLabelOffset = 12*FontScale  # How much the axis label is offset from the tick labels
StrandLabelOffset = 6*FontScale  # How must the strand label is offset from the axis
TickLabelFontSize = 24*FontScale  # Font size for tick labels
TickLabelOffset = 6  # How much the tick labels are offset from the axis
AxisLabelColor = 'black'  # Font color for axis labels
StrandLabelColor = 'black'  # Font color for strand labels
TickLabelColor = 'black'  # Font color for tick labels
TitleColor = 'black'  # Font color for title
TitleFontSize = 36*FontScale  # Font size for title
TitleOffset = 36*FontScale # How much title is offset
FreeEnergyFontSize = 10 # font size for free energy display
FreeEnergyFontSizeFrac = 0.035*FontScale # Free energy font size as fraction of total image height
                               # This is currently active, and FreeEnergyFontSize is ignored.
                               # The hack that employs this can be found by searching for
                               # FreeEnergyFontSizeFrac in the text
FreeEnergyFontSizeHack = 0.9 # This is a hack factor to get font size to match MFE
                              # drawing font size.  Actual font size typically ends up being about 37.
FreeEnergyColor = 'black' # Font color for free energy
FreeEnergyVertOffset = 64 # How much the free energy is offset from bottom of plot
FreeEnergyHorizOffset = 128 # How much free energy is offset from edge of plot


# Strand2S is the minimal fraction of bases a strand must have to be
# labeled as "Strand X" as opposed to "S1".  The setting of 0.175 works
# well for the setting for other parameters we have chosen here.
Strand2S = 0.175

# CharWidth is the approximate width of a character
# This is used to determine how many letters of a strand label
# we can write out.  The value below seems to work well.
CharWidth = 0.5*StrandLabelFontSizeStandard

LowerTriangle = 1  # Do we fill the lower triangle?  1 for yes and 0 for no
                   # For Method = 2 and Asymmetric = 0, the only way to fill
                   # it is to mirror the upper.

CutUnattached = 1 # Cut out strands that are not in complex

Asymmetric = 0 # If using .ocx-epairs files, compute probabilities in rows, meaning
               # the matrix is asymmetric

BasePairConcs = 1  # If using permAvg files, compute fraction of strands that
                   # have base pair

# The Tick Interval is now defined as:
#  N = 0-10, TickInterval = 1
#  N = 10-50, TickInterval = 5
#  N = 51-100, TickInterval = 10
#  N = 101-200, TickInterval = 20
#  N = 901-1000, TickInterval = 100
#  N = 1001-2000, TickInterval = 200
#   ...and so on

# Exit codes
INVALID_INPUT_FILE_ERROR = 2
FAILED_TO_FIND_PERMUTATION_ERROR = 3
FAILED_TO_FIND_INITIAL_CONCENTRATIONS = 4
NO_COMPLEX_ID_ERROR = 5
NO_PERM_ID_ERROR = 6

# ##### Fudge factors for vertical alignment of text
# These should vanish when we finally figure out how to align these things
# This works well for TickLabelFontSize = 24
TickHorizFudgeFact = 2.4/3.0*TickLabelFontSize
TickVertFudgeFact = 1.0/3.0*TickLabelFontSize

# This works well for ColorBarLabelFontSize = 24
ColorBarVertFudgeFact = 1.0/3.0*ColorBarLabelFontSize

# These work for StrandLabelFontSizeStandard = 36
StrandHorizFudgeFact = 1.0/3.0*StrandLabelFontSizeStandard
StrandVertFudgeFact = -7.0/3.0*StrandLabelFontSizeStandard

# These work for AxisLabelFontSize = 36
AxisHorizFudgeFact = 2.1/3.0*AxisLabelFontSize
AxisVertFudgeFact = 0.0/3.0*AxisLabelFontSize
#############################################################



# ######### GET DATA FROM COMMAND LINE INPUT ##################
# Get the name of the input file
InputFile = argv[-1]

# Get the suffix for the input file to see what kind of dot plot we're doing
InputFileList = InputFile.split('.')
InputFileSuffix = InputFileList[-1]
prefix = InputFile[0:(len(InputFile) - len(InputFileSuffix) - 1)]
cxFile = prefix + '.cx'
ocxFile = prefix + '.ocx'
permKeyFile = prefix + '.ocx-key'
if InputFileSuffix == 'ocx-ppairs':
    Method = 1  # Treat strands as if distinguishable
    NoPermID = 0   # had perm IDs
    Asymmetric = 0  # Asymmetry only allowed when Method == 2 and 3
elif InputFileSuffix == 'ocx-epairs':
    Method = 2  # strand indistinguishable, have PermIDs
    NoPermID = 0
elif InputFileSuffix == 'cx-epairs':
    Method = 2  # indistinguishable strands
    NoPermID = 1 # no perm IDs
elif InputFileSuffix == 'fpairs':
    Method = 3  # Do fraction paired style
    NoPermID = 1  # Perm ID not relevant
    Asymmetric = 1 # Must be asymmetric
elif InputFileSuffix == 'ppairs':
    Method = 4  # Input file didn't come from Complexes, no perm or complex IDs
    Asymmetric = 0 #  Must by symmetric
    NoPermID = 1  # Not relevant, there's no perm or complex ID.
# We currently do not allow dot plots for .mfe files because that's pointless
else: # Can't yet do these types of files
    print 'Sorry, we cannot yet do this type of dot plot.\nExiting...\n\n'
    exit(INVALID_INPUT_FILE_ERROR)

# Find the fontscale option

LegendLabelIndex = -1
if '-legend_label' in argv:
    LegendLabelIndex = argv.index('-legend_label') + 1
elif '--legend_label' in argv:
    LegendLabelIndex = argv.index('--legend_label') + 1

if LegendLabelIndex == -1:
    LegendLabel = "Equilibrium probability"
else:
    LegendLabel = argv[LegendLabelIndex]

# Find the output file
OutFileIndex = -1  # Initialize in case the option isn't chosen
if '-svgfile' in argv:
   OutFileIndex  = argv.index('-svgfile') + 1
elif '--svgfile' in argv:
    OutFileIndex = argv.index('--svgfile') + 1
elif '-s' in argv:
    OutFileIndex = argv.index('-s') + 1

if OutFileIndex == -1:
    OutputFile = prefix + '-dotplot.svg'
else:
    OutputFile = argv[OutFileIndex]


# Find the file with strand labels
LabelFileIndex = -1  # Initialize in case the option isn't chosen
if '-labels' in argv:
   LabelFileIndex  = argv.index('-labels') + 1
elif '--labels' in argv:
    LabelFileIndex = argv.index('--labels') + 1
elif '-l' in argv:
    LabelFileIndex = argv.index('-l') + 1

if LabelFileIndex == -1:
    UseLabels = 0
else:
    LabelFile = argv[LabelFileIndex]
    UseLabels = 1


# Whether or not to include a grid
ShowGrid = 0
if ('-grid' in argv) or ('--grid' in argv) or ('-g' in argv):
    ShowGrid = 1

# Whether or not we do long header
LongHeader = 0
if ('-longheader' in argv) or ('--longheader' in argv) or ('-h' in argv):
    LongHeader = 1


# Find the complex ID
if Method == 1 or Method == 2:
    if '-complex' in argv:
        CompIndex = argv.index('-complex') + 1
    elif '--complex' in argv:
        CompIndex = argv.index('--complex') + 1
    elif '-c' in argv:
        CompIndex = argv.index('-c') + 1
    else:
        print 'Error: No complex ID!'
        exit(NO_COMPLEX_ID_ERROR)
    ComplexID = int(argv[CompIndex])
else:
    ComplexID = -1  # Dummy complex ID


# Find the permuation ID
if NoPermID == 0:
    if '-order' in argv:
        PermIndex = argv.index('-order') + 1
    elif '--order' in argv:
        PermIndex = argv.index('--order') + 1
    elif '-o' in argv:
        PermIndex = argv.index('-o') + 1
    else:
        print 'Error: No order ID!'
        exit(NO_PERM_ID_ERROR)
    PermID = int(argv[PermIndex])
else:
    PermID = -1  # Dummy perm ID


# Check to see if we use the MFE structure
UseMFE = 0
MFEFile = 'dummy_string'
if ('-mfe' in argv) or ('--mfe' in argv) or ('-m' in argv):
    UseMFE = 1
    if Method == 4:
        MFEFile = prefix + '.mfe'
    else:
        MFEFile = prefix + '.ocx-mfe'


# Check to see if we use a target structure
UseTarget = 0  # grandfathered in, next time use True or False
UseMfeTarget = False
TargetFile = 'dummy_string'
if ('-target' in argv) or ('--target' in argv):
    UseTarget = 1
    if Method == 4:
        TargetFile = prefix + '.target'
    else:
        TargetFile = prefix + '.ocx-target'
    if Method == 2 or Method == 3:
        print '\nCannot have a target for chosen format of input file.  Ignoring target.\n'
        UseTarget = 0

if ('-mfetarget' in argv) or ('--mfetarget' in argv):
    UseTarget = 1
    UseMfeTarget =True
    if Method == 4:
        TargetFile = prefix + '.mfe'
    else:
        TargetFile = prefix + '.ocx-mfe'
    if Method == 2 or Method == 3:
        print '\nCannot have a target for chosen format of input file.  Ignoring target.\n'
        UseTarget = 0

# Check for title
UseTitle = 0
if '-title' in argv:
    TitleStrIndex = argv.index('-title') + 1
    UseTitle = 1
elif '--title' in argv:
    TitleStrIndex = argv.index('--title') + 1
    UseTitle = 1
elif '-t' in argv:
    TitleStrIndex = argv.index('-t') + 1
    UseTitle = 1
if  UseTitle:
    TitleString = argv[TitleStrIndex]
else:
    TitleString = ''  # Dummy place holder


# Check for free energy buffering (square overrides)
UseFreeEnergyBuffer = 0
if '-freeenergybuffer' in argv:
    UseFreeEnergyBuffer = 1
elif '--freeenergybuffer' in argv:
    UseFreeEnergyBuffer = 1
elif '-b' in argv:
    UseFreeEnergyBuffer = 1


# Check for free energy (buffer if necessary)
ShowFreeEnergy = 0
FreeEnergyInFile = 0
if Method != 3:
    if '-freeenergy' in argv:
        FreeEnergyInFile = 1
        ShowFreeEnergy = 1
    elif '--freeenergy' in argv:
        FreeEnergyInFile = 1
        ShowFreeEnergy = 1
    elif '-f' in argv:
        FreeEnergyInFile = 1
        ShowFreeEnergy = 1


    if '-energy' in argv:
        FreeEnergyIndex  = argv.index('-energy') + 1
        ShowFreeEnergy = 1
        FreeEnergyInFile = 0
    elif '--energy' in argv:
        FreeEnergyIndex = argv.index('--energy') + 1
        ShowFreeEnergy = 1
        FreeEnergyInFile = 0
    elif '-e' in argv:
        FreeEnergyIndex = argv.index('-e') + 1
        ShowFreeEnergy = 1
        FreeEnergyInFile = 0

    if ShowFreeEnergy:
        UseFreeEnergyBuffer = 1

        # Get the free energy off the command line if not in the file
        if not(FreeEnergyInFile):
            FreeEnergy = float(argv[FreeEnergyIndex])

        # Check for units of free energy
        UnitsString = 'kcal/mol' # default
        FreeEnergyUnitsIndex = -1
        if '-units' in argv:
            FreeEnergyUnitsIndex  = argv.index('-units') + 1
        elif '--units' in argv:
            FreeEnergyUnitsIndex = argv.index('--units') + 1
        elif '-u' in argv:
            FreeEnergyUnitsIndex = argv.index('-u') + 1

        if FreeEnergyUnitsIndex != -1:
            UnitsString = argv[FreeEnergyUnitsIndex]

        # Check for free energy preamble
        UseFreeEnergyPreamble = 0
        if '-energypreamble' in argv:
            FreeEnergyPreambleIndex  = argv.index('-energypreamble') + 1
            UseFreeEnergyPreamble = 1
        elif '--energypreamble' in argv:
            FreeEnergyPreambleIndex = argv.index('--energypreamble') + 1
            UseFreeEnergyPreamble = 1
        elif '-p' in argv:
            FreeEnergyPreambleIndex = argv.index('-p') + 1
            UseFreeEnergyPreamble = 1

        if UseFreeEnergyPreamble:
            PreambleString = argv[FreeEnergyPreambleIndex]


# Find whether or not we do a data file
DataFileIndex = -1  # Initialize in case the option isn't chosen
if '-datafile' in argv:
   DataFileIndex  = argv.index('-datafile') + 1
elif '--datafile' in argv:
    DataFileIndex = argv.index('--datafile') + 1
elif '-d' in argv:
    DataFileIndex = argv.index('-d') + 1

if DataFileIndex == -1 or Method == 2:
    OutputDataFile = 0
else:
    OutputDataFile = 1
    DataFile = argv[DataFileIndex]
if Method == 2 and DataFileIndex != -1:
    print 'No data file written.  Do not yet have capability for files of type .%s' % InputFileSuffix


# Check to see if -mfe flag is to be ignored
if UseMFE and (Asymmetric or Method == 2 or Method == 3):
    print 'Ignoring -mfe flag.  Not relevant for given input file.'
    UseMFE = 0

# #############################################################


# Read in strand labels
if  UseLabels:
    LabelStrings = []
    f = getInputData.delay_open(LabelFile,'r')
    line = f.readline()
    # Blow through comments
    while (line[0] == '%' or line[0] == '\n') and line != '':
        line = f.readline()

    # Read in strand labels
    while line != '':
        # Chop off newline and add to label list
        if line[-1] == '\n':
            LabelString = line[0:-1]
        LabelStrings.append(LabelString)
        line = f.readline()
    f.close()


# Read in and store comments form input file
if LongHeader:
    CommentList = []
    f = getInputData.delay_open(InputFile,'r')
    line = f.readline()
    while (line[0] == '%' or line[0] == '\n') and line != '':
        CommentList.append(line)
        line = f.readline()
    f.close()


# Make error messages for getInputData
Errors = [FAILED_TO_FIND_PERMUTATION_ERROR,FAILED_TO_FIND_INITIAL_CONCENTRATIONS]


# Reads in the input data
# Make parameters to sent to getInputData
Params = [Method,NoPermID,Asymmetric,LowerTriangle,CutUnattached,UseMFE,FreeEnergyInFile]

if Method < 4:
    PairProb,UnPairProb,Sequences,PairList,N,K,permKey,Mirror,FileFreeEnergy = \
               getInputDataComplexes(InputFile,cxFile,ocxFile,permKeyFile,MFEFile,ComplexID,PermID,Params,Errors)
else:
    PairProb,UnPairProb,Sequences,PairList,N,K,permKey,Mirror,FileFreeEnergy = \
              getInputDataNuPackLib(InputFile,MFEFile,LowerTriangle,UseMFE,FreeEnergyInFile,Errors)

if UseTarget or UseMfeTarget:
    TargetList,TargetUnpaired = \
              getInputDataTarget(TargetFile, permKey, Sequences, N, K,
              ComplexID, PermID, Params, Errors, UseMfeTarget)


# Set the free energy equal to the file free energy
if FreeEnergyInFile:
    FreeEnergy = FileFreeEnergy


if ShowFreeEnergy:
    # Build the free energy string
    if UseFreeEnergyPreamble:
        FreeEnergyString = PreambleString + ' '
    else:
        FreeEnergyString = ''

    FreeEnergyString += '%.2f %s' % (FreeEnergy,UnitsString)



# At this point, the only data we have is in the form of:
# PairProb, UnPairProb, Sequences, PairList, TargetList, TargetUnpaired, N, permKey,
# FreeEnergyString
# Everything else is parameters for plotting

# ################ WRITE DATA TO OUTPUT FILE ##################
if OutputDataFile:
    writeDataFile(Sequences,PairProb,UnPairProb,permKey,Method,UseTitle,TitleString,DataFile)
# #############################################################


###############################################################
# AT THIS POINT, ALL THE DATA ARE ORGANIZED AND READY FOR     #
# USE IN GENERATING THE DOT PLOT.  NO CODE BELOW THIS POINT   #
# MESSES WITH ANY PROBABILITIES, PARSING OF FILES, ETC.       #
#                                                             #
#         HERE BEGINS THE SCRIPT FOR SVG GENERATION           #
###############################################################


# Get the number of strands
Nstrands = len(permKey)


# Get the appropriate top margin
if UseTitle == 0 and Nstrands > 1:  # No title, strand labels
    TopMargin = TopMarginFull * TopMarginFrac01
else: # Title and strand labels
    TopMargin = TopMarginFull

# COMMENT BELOW OUT TO INCLUDE LABELS FOR SINGLE STRANDED CASE
if UseTitle == 0 and Nstrands == 1: # No strand labels, no title
    TopMargin = TopMarginFull * TopMarginFrac00
elif UseTitle == 1 and Nstrands == 1: # No strand labels, title
    TopMargin = TopMarginFull * TopMarginFrac10
elif UseTitle == 0 and Nstrands > 1:  # No title, strand labels
    TopMargin = TopMarginFull * TopMarginFrac01
else: # Title and strand labels
    TopMargin = TopMarginFull


# The coordinates for the corner of the box (the axes) containing the plot
BoxCornerX = LeftMargin
BoxCornerY = TopMargin


# BaseWidth is width a base occupies in dot plot
BaseWidth = 1.0*(d - 2*(Nstrands-1)*AxisWidth)/(N + (Nstrands-1)*(StrandBreakGap))
TargetLineWidth = TargetLineWidthFrac * BaseWidth


# Convert size ratios from area to side lengths (more convenient)
if AreaScaling == 1:
    MediumBigRatio = sqrt(MediumBigRatio)
    SmallBigRatio = sqrt(SmallBigRatio)


# ################## GET DOT POSITIONS #######################
# Note that in the case where we're asymmetrical, the upper and lower
# triangle are both taken care of by this part of the code.
# We append each entry in PairProb to include the x and y coordinates of the
# lower left corner of the square on the dot plot and its dimension (length
# of a side)
# After this is done, one entry x in PairProb has:
# x[0] = base i of the pair i,j
# x[1] = base j of base j in the pair with j > i
# x[2] = probability that i and j are paired
# x[3] = number of strand breaks to the left of base i
# x[4] = number of strand breaks to the left of base j
# x[5] = left corner of the square on the dotplot
# x[6] = right corner of the square on the dotplot
# x[7] = the length of a side of the square
# x[8] = the color of the square in rgb

if DiscreteScaling == 0:
    for x in PairProb:
        SquareDim = BaseWidth
        if BoxScaling == 1:
            if AreaScaling == 1:
                SquareDim *= sqrt(x[2])
            else:
                SquareDim *= x[2]
        LeftCornerX = BoxCornerX + x[1]*BaseWidth + x[4]*(StrandBreakGap*BaseWidth + 2*AxisWidth) \
                      + (BaseWidth-SquareDim)/2.0
        LeftCornerY = BoxCornerY + x[0]*BaseWidth + x[3]*(StrandBreakGap*BaseWidth + 2*AxisWidth) \
                      + (BaseWidth-SquareDim)/2.0
        red,green,blue = getRGB(x[2])
        x.append(LeftCornerX)
        x.append(LeftCornerY)
        x.append(SquareDim)
        x.append([red,green,blue])
else: # Discrete scaling
    for x in PairProb:
        SquareDim = BaseWidth
        if x[2] >= minP[0]: # Biggest square
            Color = BigColor
        elif x[2] >= minP[1]:  # Medium square
            Color = MediumColor
            if BoxScaling == 1:
                SquareDim *= MediumBigRatio
        elif x[2] >= minP[2]:  # Small square
            Color = SmallColor
            if BoxScaling == 1:
                SquareDim *= SmallBigRatio
        # Otherwise, probability is too small

        LeftCornerX = BoxCornerX + x[1]*BaseWidth + x[4]*(StrandBreakGap*BaseWidth + 2*AxisWidth) \
                      + (BaseWidth-SquareDim)/2.0
        LeftCornerY = BoxCornerY + x[0]*BaseWidth + x[3]*(StrandBreakGap*BaseWidth + 2*AxisWidth) \
                      + (BaseWidth-SquareDim)/2.0
        x.append(LeftCornerX)
        x.append(LeftCornerY)
        x.append(SquareDim)
        x.append(Color)
# #############################################################


# ######## GET DOT POSITIONS FOR TARGET #######################
if UseTarget:
    for x in TargetList:
        SquareDim = BaseWidth
        LeftCornerX = BoxCornerX + x[1]*BaseWidth + x[3]*(StrandBreakGap*BaseWidth + 2*AxisWidth) \
                      + (BaseWidth-SquareDim)/2.0
        LeftCornerY = BoxCornerY + x[0]*BaseWidth + x[2]*(StrandBreakGap*BaseWidth + 2*AxisWidth) \
                      + (BaseWidth-SquareDim)/2.0
        x.append(LeftCornerX)
        x.append(LeftCornerY)
# #############################################################


# ######### GET DOT POS. FOR LOWER TRI. STRUCT. ###############
if Asymmetric == 0 and LowerTriangle == 1 and Mirror == 0:
    # We append each entry in PairList to include the x and y coordinates of the
    # lower left corner of the square on the dot plot and its dimension (length
    # of a side).  All squares in other structure are maximal size.
    # When completed, entry x in PairList is:
    # x[0] = base i of the pair i,j
    # x[1] = base j of base j in the pair with j > i
    # x[2] = dummy place holder to maintain consistency with PairProb
    # x[3] = number of strand breaks to the left of base i
    # x[4] = number of strand breaks to the left of base j
    # x[5] = left corner of the square on the dotplot
    # x[6] = right corner of the square on the dotplot
    # x[7] = the length of a side of the square
    # x[8] = the color of the square
    for x in PairList:
        LeftCornerX = BoxCornerX + x[1]*BaseWidth + x[4]*(StrandBreakGap*BaseWidth + 2*AxisWidth)
        LeftCornerY = BoxCornerY + x[0]*BaseWidth + x[3]*(StrandBreakGap*BaseWidth + 2*AxisWidth)
        SquareDim = BaseWidth
        x.append(LeftCornerX)
        x.append(LeftCornerY)
        x.append(SquareDim)
        if DiscreteScaling == 0:
            x.append(OtherColorContinuous)
        else:
            x.append(OtherColorDiscrete)
###############################################################


# ############# GET DOT POS. FOR LOWER TRI. MIRROR ############
# We append each entry in PairList to include the x and y coordinates of the
# lower left corner of the square on the dot plot and its dimension (length
# of a side)
# After this is done, one entry x in PairList has:
# x[0] = base i of the pair i,j
# x[1] = base j of base j in the pair with j > i
# x[2] = probability that i and j are paired
# x[3] = number of strand breaks to the left of base i
# x[4] = number of strand breaks to the left of base j
# x[5] = left corner of the square on the dotplot
# x[6] = right corner of the square on the dotplot
# x[7] = the length of a side of the square
# x[8] = the color of the square
if LowerTriangle == 1 and Mirror == 1:
    if DiscreteScaling == 0:
        for x in PairList:
            SquareDim = BaseWidth
            if BoxScaling == 1:
                if AreaScaling == 1:
                    SquareDim *= sqrt(x[2])
                else:
                    SquareDim *= x[2]
            LeftCornerX = BoxCornerX + x[1]*BaseWidth + x[4]*(StrandBreakGap*BaseWidth + 2*AxisWidth) \
                          + (BaseWidth-SquareDim)/2.0
            LeftCornerY = BoxCornerY + x[0]*BaseWidth + x[3]*(StrandBreakGap*BaseWidth + 2*AxisWidth) \
                          + (BaseWidth-SquareDim)/2.0
            red,green,blue = getRGB(x[2])
            x.append(LeftCornerX)
            x.append(LeftCornerY)
            x.append(SquareDim)
            x.append([red,green,blue])
    else: # Discrete scaling
        for x in PairList:
            SquareDim = BaseWidth
            if x[2] >= minP[0]: # Biggest square
                Color = BigColor
            elif x[2] >= minP[1]:  # Medium square
                Color = MediumColor
                if BoxScaling == 1:
                    SquareDim *= MediumBigRatio
            elif x[2] >= minP[2]:  # Small square
                Color = SmallColor
                if BoxScaling == 1:
                    SquareDim *= SmallBigRatio
            # Otherwise, probability is too small

            LeftCornerX = BoxCornerX + x[1]*BaseWidth + x[4]*(StrandBreakGap*BaseWidth + 2*AxisWidth) \
                          + (BaseWidth-SquareDim)/2.0
            LeftCornerY = BoxCornerY + x[0]*BaseWidth + x[3]*(StrandBreakGap*BaseWidth + 2*AxisWidth) \
                          + (BaseWidth-SquareDim)/2.0
            x.append(LeftCornerX)
            x.append(LeftCornerY)
            x.append(SquareDim)
            x.append(Color)
# #############################################################


# ############### GET DOTS FOR UNPAIRED PROBS #################
# Where the left corner should start on unpaired dots
LeftXUnPair = BoxCornerX + d + AxisWidth + BaseWidth + AxisWidth

if DiscreteScaling == 0:
    for x in UnPairProb:
        SquareDim = BaseWidth
        if SquareDim>40:
          SquareDim==40
        if BoxScaling == 1:
            if AreaScaling == 1:
                SquareDim *= sqrt(x[1])
            else:
                SquareDim *= x[1]
        LeftCornerX = LeftXUnPair + (BaseWidth-SquareDim)/2.0
        LeftCornerY = BoxCornerY + x[0]*BaseWidth + x[2]*(2*AxisWidth + StrandBreakGap*BaseWidth) \
                      + (BaseWidth-SquareDim)/2.0
        red,green,blue = getRGB(x[1])
        x.append(LeftCornerX)
        x.append(LeftCornerY)
        x.append(SquareDim)
        x.append([red,green,blue])
else: # Discrete scaling
    for x in UnPairProb:
        SquareDim = BaseWidth
        if SquareDim>40:
          SquareDim==40
        if x[1] >= minP[0]: # Biggest square
            Color = BigColor
        elif x[1] >= minP[1]:  # Medium square
            Color = MediumColor
            if BoxScaling == 1:
                SquareDim *= MediumBigRatio
        elif x[1] >= minP[2]:  # Small square
            Color = SmallColor
            if BoxScaling == 1:
                SquareDim *= SmallBigRatio
        # Otherwise, probability is too small

        LeftCornerX = LeftXUnPair + (BaseWidth-SquareDim)/2.0
        LeftCornerY = BoxCornerY + x[0]*BaseWidth + x[2]*(2*AxisWidth + StrandBreakGap*BaseWidth) \
                      + (BaseWidth-SquareDim)/2.0
        x.append(LeftCornerX)
        x.append(LeftCornerY)
        x.append(SquareDim)
        x.append(Color)
# #############################################################


# ########## GET DOTS FOR UNPAIRED TARGET #####################
if UseTarget:
    for x in TargetUnpaired:
        SquareDim = BaseWidth
        LeftCornerX = LeftXUnPair + (BaseWidth-SquareDim)/2.0
        LeftCornerY = BoxCornerY + x[0]*BaseWidth + x[1]*(2*AxisWidth + StrandBreakGap*BaseWidth) \
                      + (BaseWidth-SquareDim)/2.0
        x.append(LeftCornerX)
        x.append(LeftCornerY)
# #############################################################


# ########## DETERMINE IF PAIR IS PART OF TARGET ##############
if UseTarget:
    for x in PairProb:
        InTarget = 0
        NpairsInTarget = len(TargetList)
        i = 0
        while InTarget == 0 and i < NpairsInTarget:
            if x[0] == TargetList[i][0] and x[1] == TargetList[i][1]:
                InTarget = 1
            i += 1
        x.append(InTarget)
    for x in UnPairProb:
        InTarget = 0
        NUnpairedInTarget = len(TargetUnpaired)
        i = 0
        while InTarget == 0 and i < NUnpairedInTarget:
            if x[0] == TargetUnpaired[i][0]:
                InTarget = 1
            i += 1
        x.append(InTarget)
# #############################################################


# ############### SET UP AXES FOR PLOTTING ####################
# AxisData[i] has [upper left x-coord, upper left corner y-coord, x-dimension, y-dimension]
# for the given box in the plot
AxisData = []
CornerY = BoxCornerY - AxisWidth/2.0
for x in permKey:
    CornerX = BoxCornerX - AxisWidth/2.0
    for y in permKey:
        BoxDimX = BaseWidth*Sequences[y-1][2] + AxisWidth
        BoxDimY = BaseWidth*Sequences[x-1][2] + AxisWidth
        AxisData.append([CornerX,CornerY,BoxDimX,BoxDimY])
        CornerX += BaseWidth*Sequences[y-1][2] + 2*AxisWidth + StrandBreakGap*BaseWidth
    CornerY += BaseWidth*Sequences[x-1][2] + 2*AxisWidth + StrandBreakGap*BaseWidth
###############################################################


# ############### SET UP DIAGONAL LINE ########################
if Asymmetric == 0 or K == 1: # Only include if matrix is symmetric
    # DiagData contains the x,y start and end points for diagonal lines
    DiagData = []
    x1 = BoxCornerX - AxisWidth/2.0
    y1 = BoxCornerY - AxisWidth/2.0
    for x in permKey:
        x2 = x1 + BaseWidth*Sequences[x-1][2] + AxisWidth
        y2 = y1 + BaseWidth*Sequences[x-1][2] + AxisWidth
        DiagData.append([x1,y1,x2,y2])
        x1 = x1 + BaseWidth*Sequences[x-1][2] + 2*AxisWidth + StrandBreakGap*BaseWidth
        y1 = y1 + BaseWidth*Sequences[x-1][2] + 2*AxisWidth + StrandBreakGap*BaseWidth
###############################################################


# ############ UNPAIRED PROBABILITY BOUNDING BOX ##############
UnPairedBoxData = []
x1 = BoxCornerX + d + AxisWidth + BaseWidth + AxisWidth/2.0
y1 = BoxCornerY - AxisWidth/2.0
for x in permKey:
    BoxDimX = BaseWidth + AxisWidth
    BoxDimY = BaseWidth*Sequences[x-1][2] + AxisWidth
    UnPairedBoxData.append([x1,y1,BoxDimX,BoxDimY])
    y1 = y1 + BaseWidth*Sequences[x-1][2] + 2*AxisWidth + StrandBreakGap*BaseWidth
###############################################################


# ################# SET UP TICKS FOR PLOT #####################
# The Tick Interval is now defined as:
#  N = 0-10, TickInterval = 1
#  N = 11-20, TickInterval = 2
#  N = 11-50, TickInterval = 5
#  N = 51-100, TickInterval = 10
#  N = 101-200, TickInterval = 20
#  N = 901-1000, TickInterval = 100
#  N = 1001-2000, TickInterval = 200
#   ...and so on
virtualN = d / BaseWidth
if virtualN <= 10:
    TickInterval = 1
elif virtualN <= 20:
    TickInterval = 2
elif virtualN <= 50:
    TickInterval = 5
elif virtualN <= 100:
    TickInterval = 10
else:
    pow10 = floor(log10(virtualN) + 0.01)
    TickInterval = ceil(virtualN/10.0**pow10) * 10**(pow10-1)


# Each entry in Ticks has [x1,y1,x2,y2], the coords for tick lines
Ticks = []
XTickLabels = []
YTickLabels = []


# First do ticks on horizontal axes
CornerX = BoxCornerX
x1 = CornerX - BaseWidth/2.0
yLabel = BoxCornerY + d + AxisWidth + TickLabelOffset + TickHorizFudgeFact
ybottom1 = BoxCornerY + d
ybottom2 = ybottom1 - TickLengthFrac*d
ytop1 = BoxCornerY
ytop2 = BoxCornerY + TickLengthFrac*d
for x in permKey:
    index = TickInterval
    while index <= Sequences[x-1][2]:
        x1 += TickInterval*BaseWidth
        Ticks.append([x1,ybottom1,x1,ybottom2])
        Ticks.append([x1,ytop1,x1,ytop2])
        XTickLabels.append([x1,yLabel,index])
        index += TickInterval
    if LastTick == 1:
        index -= TickInterval  # Bring us back to where we were
        if Sequences[x-1][2] - index > LastTickFrac*TickInterval: # include last tick
            x1 += (Sequences[x-1][2] - index)*BaseWidth
            Ticks.append([x1,ybottom1,x1,ybottom2])
            Ticks.append([x1,ytop1,x1,ytop2])
            XTickLabels.append([x1,yLabel,Sequences[x-1][2]])
    CornerX += BaseWidth*Sequences[x-1][2] + 2*AxisWidth + StrandBreakGap*BaseWidth
    x1 = CornerX - BaseWidth/2.0


# Ticks on vertical axes
CornerY = BoxCornerY
y1 = CornerY - BaseWidth/2.0
xLabel = BoxCornerX - AxisWidth - TickLabelOffset
xleft1 = BoxCornerX
xleft2 = xleft1 + TickLengthFrac*d
xright1 = BoxCornerX + d
xright2 = xright1 - TickLengthFrac*d
for x in permKey:
    index = TickInterval
    while index <= Sequences[x-1][2]:
        y1 += TickInterval*BaseWidth
        Ticks.append([xleft1,y1,xleft2,y1])
        Ticks.append([xright1,y1,xright2,y1])
        YTickLabels.append([xLabel,y1 + TickVertFudgeFact,index])
        index += TickInterval
    if LastTick == 1:
        index -= TickInterval  # Bring us back to where we were
        if Sequences[x-1][2] - index > LastTickFrac*TickInterval: # include last tick
            y1 += (Sequences[x-1][2] - index)*BaseWidth
            Ticks.append([xleft1,y1,xleft2,y1])
            Ticks.append([xright1,y1,xright2,y1])
            YTickLabels.append([xLabel,y1 + TickVertFudgeFact,Sequences[x-1][2]])
    CornerY += BaseWidth*Sequences[x-1][2] + 2*AxisWidth + StrandBreakGap*BaseWidth
    y1 = CornerY - BaseWidth/2.0
###############################################################


# ################ SET UP GRIDLINES FOR PLOT ##################
# Each entry in GridLines has [x1,y1,x2,y2], the coords for start and end of gridline
GridLines = []
if ShowGrid:
    # First do horizontal grid lines
    CornerX = BoxCornerX

    HorizStart = BoxCornerX
    VertStart = BoxCornerY

    for x in permKey:
        HorizEnd = HorizStart + BaseWidth*Sequences[x-1][2]
        for y in permKey:
            VertEnd = VertStart + BaseWidth*Sequences[y-1][2]

            # Vertical grid lines
            x1 = HorizStart - BaseWidth/2.0
            index = TickInterval
            while index <= Sequences[x-1][2]:
                x1 += TickInterval*BaseWidth
                GridLines.append([x1,VertStart,x1,VertEnd])
                index += TickInterval
            if LastTick == 1:
                index -= TickInterval  # Bring us back to where we were
                if Sequences[x-1][2] - index > LastTickFrac*TickInterval: # include last tick
                    x1 += (Sequences[x-1][2] - index)*BaseWidth
                    GridLines.append([x1,VertStart,x1,VertEnd])

            # Horizontal grid lines
            y1 = VertStart - BaseWidth/2.0
            index = TickInterval
            while index <= Sequences[y-1][2]:
                y1 += TickInterval*BaseWidth
                GridLines.append([HorizStart,y1,HorizEnd,y1])
                index += TickInterval
            if LastTick == 1:
                index -= TickInterval  # Bring us back to where we were
                if Sequences[y-1][2] - index > LastTickFrac*TickInterval: # include last tick
                    y1 += (Sequences[y-1][2] - index)*BaseWidth
                    GridLines.append([HorizStart,y1,HorizEnd,y1])

            VertStart = VertEnd + 2*AxisWidth + StrandBreakGap*BaseWidth
        VertStart = BoxCornerY
        HorizStart = HorizEnd + 2*AxisWidth + StrandBreakGap*BaseWidth
###############################################################



# ################# SET UP STRAND LABELS ######################
if Nstrands > 1:  # If we have more than one strand -- CHANGE TO >= IF WE ALWAYS WANT LABELS
    StrandLabelFontSize = StrandLabelFontSizeStandard # Start with standard size
    XStrandLabels = []  # List of coordinates and strings for strand labels on y-axis
    YStrandLabels = []  # List of coordinates and strings for strand labels on x-axis
    xhoriz0 = BoxCornerX - AxisWidth/2.0
    yvert0 = BoxCornerY - AxisWidth/2.0
    xvert = BoxCornerX + d + AxisWidth + BaseWidth + AxisWidth + BaseWidth + AxisWidth \
            + StrandLabelFontSizeStandard + StrandLabelOffset + StrandVertFudgeFact/4.0
    yhoriz = BoxCornerY - AxisWidth - StrandLabelOffset - StrandHorizFudgeFact


    # Set up array for counts of strand types in permKey
    permKeyCount = []
    StrandCount = []
    for i in range(max(permKey)):
        permKeyCount.append(0)
        StrandCount.append(0)
    for x in permKey:
        permKeyCount[x-1] += 1


    for x in permKey:
        xhoriz = xhoriz0  + Sequences[x-1][2]*BaseWidth/2.0 + AxisWidth
        yvert = yvert0 + Sequences[x-1][2]*BaseWidth/2.0 + AxisWidth

        StrandCount[x-1] += 1
        if permKeyCount[x-1] == 1:
            StrandNumberStr = str(x)
        else:
            StrandNumberStr = str(x) + '.' + str(StrandCount[x-1])

        if UseLabels:
            if AdjustStrandLabelFontSize:
                Nchars = floor(Sequences[x-1][2]*BaseWidth/CharWidth)
                StrandStr = LabelStrings[x-1]
                if AllStrandLabelsSameSize:
                    newStrandLabelFontSize = StrandLabelFontSizeStandard*Nchars/float(len(StrandStr))
                    if newStrandLabelFontSize < StrandLabelFontSize:
                        StrandLabelFontSize = newStrandLabelFontSize
                else:
                    if Nchars < len(StrandStr):
                        StrandLabelFontSize = StrandLabelFontSizeStandard*Nchars/float(len(StrandStr))
            else:
                Nchars = int(floor(Sequences[x-1][2]*BaseWidth/CharWidth))
                if Nchars < len(LabelStrings[x-1]):
                    StrandStr = LabelStrings[x-1][0:Nchars]
                else:
                    StrandStr = LabelStrings[x-1]
        else:
            if Sequences[x-1][2] > Strand2S*virtualN: # Enough room to write "Strand x"
                StrandStr = 'Strand ' + StrandNumberStr
            else:
                StrandStr = 'S' + StrandNumberStr

        XStrandLabels.append([xhoriz,yhoriz,StrandStr,StrandLabelFontSize])
        YStrandLabels.append([xvert,yvert,StrandStr,StrandLabelFontSize])
        xhoriz0 = xhoriz + Sequences[x-1][2]*BaseWidth/2 + AxisWidth/2.0 + StrandBreakGap*BaseWidth \
                  + AxisWidth/2.0
        yvert0 = yvert + Sequences[x-1][2]*BaseWidth/2 + AxisWidth/2.0 + StrandBreakGap*BaseWidth \
                 + AxisWidth/2.0
# #############################################################


# ################# SET UP AXIS LABELS ########################
# Horizontal Label
x1 = BoxCornerX + d/2.0
y1 = BoxCornerY + d + AxisWidth + TickLabelOffset + TickHorizFudgeFact + TickLabelFontSize \
     + AxisLabelOffset + AxisHorizFudgeFact
XAxisLabel = [x1,y1]


# Vertical Label
x1 = BoxCornerX - AxisWidth - TickLabelOffset - TickHorizFudgeFact - TickLabelFontSize \
     - AxisLabelOffset
y1 = BoxCornerY + d/2.0
YAxisLabel=[x1,y1]

x2 = BoxCornerX + d + AxisWidth + TickLabelOffset + TickHorizFudgeFact + \
     AxisLabelOffset + 3*AxisLabelFontSize
if Nstrands>1:
  x2+=2*AxisLabelFontSize
RightYAxisLabel=[x2,y1]

# #############################################################


# ################# SET UP FREE ENERGY DISPLAY ################
if ShowFreeEnergy:
    XFreeEnergy = BoxCornerX + FreeEnergyHorizOffset
    YFreeEnergy = XAxisLabel[1] + FreeEnergyVertOffset + FreeEnergyFontSize
# #############################################################


# #################### BUILD COLORBAR ########################
if Nstrands>1:
  ColorBarOffset = 120  # How far the colorbar is offset from strand labels
else:
  ColorBarOffset = 60
if DiscreteScaling == 0:
    ColorBarHeight = d*ColorBarFrac
    x1 = BoxCornerX + d + AxisWidth + BaseWidth + AxisWidth + BaseWidth + AxisWidth \
         + StrandLabelFontSizeStandard + StrandLabelOffset + StrandVertFudgeFact + ColorBarOffset
    y1 = BoxCornerY - AxisWidth/2.0 + (1.0-ColorBarFrac)/2.0*d
    ColorBar = []
    RightYAxisLabel=[x1+80,RightYAxisLabel[1]]
    for x in range(NColors):
        P = float(x)/NColors
        red,green,blue = getRGB(1.0-P)
        ColorBar.append([x1,y1+P*ColorBarHeight,ColorBarWidth,d/float(NColors),[red,green,blue]])
# #############################################################


# ################ SET UP COLORBAR LABELS #####################
if DiscreteScaling == 0:
    # entry x in CBLabelData has:
    # x[0] = x-position of label
    # x[1] = y-position of label
    # x[2] = number for label
    CBLabelData = []
    x1 = BoxCornerX + d + AxisWidth + BaseWidth + AxisWidth + BaseWidth + AxisWidth \
         + StrandLabelFontSizeStandard + StrandLabelOffset + StrandVertFudgeFact + ColorBarOffset \
         + ColorBarWidth + ColorBarLabelOffset
    y1 = BoxCornerY - AxisWidth/2.0 + (1.0-ColorBarFrac)/2.0*d + ColorBarVertFudgeFact
    for i in range(NColorBarLabels):
        y = y1 + i*ColorBarHeight/float(NColorBarLabels-1.0)
        CBLabelData.append([x1,y,1.0-i*1.0/float(NColorBarLabels-1.0)])
# #############################################################


# ########### BUILD BACKGROUND COLORATIONS ###################
# BGData[i] has [upper left x-coord, upper left corner y-coord, x-dimension, y-dimension]
# for the given box in the plot
if BoxScaling == 0:
    BGData = []
    CornerY = BoxCornerY
    for x in permKey:
        CornerX = BoxCornerX
        for y in permKey:
            BoxDimX = BaseWidth*Sequences[y-1][2]
            BoxDimY = BaseWidth*Sequences[x-1][2]
            BGData.append([CornerX,CornerY,BoxDimX,BoxDimY])
            CornerX += BaseWidth*Sequences[y-1][2] + 2*AxisWidth + StrandBreakGap*BaseWidth
        CornerY += BaseWidth*Sequences[x-1][2] + 2*AxisWidth + StrandBreakGap*BaseWidth
# #############################################################


# ########### FIND DIMENSIONS OF VIEWBOX AND DOCUMENT #########
if DiscreteScaling == 0: # Uses a colorbar
    if Nstrands == 1:
      vbw = RightYAxisLabel[0] + ColorBarWidth + RightMargin/2
    else:
      vbw = RightYAxisLabel[0] + ColorBarWidth + RightMargin/2+TickLabelFontSize
else: # No color bar
    if Nstrands == 1:   # CHANGE THIS TO Nstrands < 1 if we want to have labels even for just one strand
        vbw = BoxCornerX + d + AxisWidth + BaseWidth + AxisWidth + BaseWidth + AxisWidth \
              + TickLabelOffset + RightMargin
    else:
        vbw = BoxCornerX + d + AxisWidth + BaseWidth + AxisWidth + BaseWidth + \
              AxisWidth +StrandLabelFontSizeStandard + StrandLabelOffset + \
              StrandVertFudgeFact + RightMargin +2*TickLabelFontSize

# Adjust height depending on presence of titles and strand labels
# UNCOMMENT THIS IF ALL PLOTS HAVE STRAND LABELS, EVEN FOR Nstrands = 1
#if UseTitle == 0:
#    vbh -= (TitleOffset + TitleFontSize + StrandHorizFudgeFact)

# COMMENT OUT IF STATEMENTS BELOW IF WE WANT STRAND LABELS EVEN WHEN WE ONLY HAVE ONE STRAND
if Nstrands == 1 and UseTitle:
    vbh -= (StrandLabelOffset + StrandLabelFontSizeStandard)
elif Nstrands > 1 and UseTitle == 0:
    vbh -= (TitleOffset + TitleFontSize + StrandHorizFudgeFact)
elif Nstrands == 1 and UseTitle == 0:
    vbh -= (StrandLabelOffset + StrandHorizFudgeFact + StrandLabelFontSizeStandard + TitleOffset + TitleFontSize)

if UseFreeEnergyBuffer:
    vbh += FreeEnergyVertOffset + FreeEnergyFontSize
SVGWidth = float(vbw)/vbh*SVGHeight
# #############################################################


# ##################### SET UP TITLE ##########################
if UseTitle:
    XTitle = vbw/2.0;
    YTitle = BoxCornerY - AxisWidth
    if Nstrands > 1:  # REVERT TO Nstrands >= 1 is we have strand labels for single-strand case
        YTitle -= (StrandLabelOffset + StrandHorizFudgeFact + StrandLabelFontSizeStandard)
    YTitle -= TitleOffset

    # Fix the title string to have degree symbols (replace degC and degF)
    TitleString = TitleString.replace('degC','&#x2103;')
# #############################################################


# ########## MAKE SQUARE AND ILLUSTRATOR COMPLIANT ############
if MakeSquare:
    if vbh > vbw:  # height is greater than width
        TranslateDist = (vbh-vbw)/2.0
        if BoxScaling == 0 and DiscreteScaling == 0:
            for i in range(len(BGData)):
                BGData[i][0] += TranslateDist
        if DiscreteScaling == 0:
            for i in range(len(ColorBar)):
                ColorBar[i][0] += TranslateDist
            for i in range(len(CBLabelData)):
                CBLabelData[i][0] += TranslateDist
        for i in range(len(AxisData)):
            AxisData[i][0] += TranslateDist
        for i in range(len(UnPairedBoxData)):
            UnPairedBoxData[i][0] += TranslateDist
        if Asymmetric == 0 or K == 1:
            for i in range(len(DiagData)):
                DiagData[i][0] += TranslateDist
                DiagData[i][2] += TranslateDist
        for i in range(len(Ticks)):
            Ticks[i][0] += TranslateDist
            Ticks[i][2] += TranslateDist
        if ShowGrid:
            for i in range(len(GridLines)):
                GridLines[i][0] += TranslateDist
                GridLines[i][2] += TranslateDist
        for i in range(len(XTickLabels)):
            XTickLabels[i][0] += TranslateDist
        for i in range(len(YTickLabels)):
            YTickLabels[i][0] += TranslateDist
        if Nstrands > 1:
            for i in range(len(XStrandLabels)):
                XStrandLabels[i][0] += TranslateDist
            for i in range(len(YStrandLabels)):
                YStrandLabels[i][0] += TranslateDist
        if UseTitle:
            XTitle += TranslateDist
        XAxisLabel[0] += TranslateDist
        YAxisLabel[0] += TranslateDist
        RightYAxisLabel[0] += TranslateDist
        if ShowFreeEnergy:
            XFreeEnergy += TranslateDist
        for i in range(len(PairProb)):
            PairProb[i][5] += TranslateDist
        if LowerTriangle == 1 and Asymmetric == 0:
            for i in range(len(PairList)):
                PairList[i][5] += TranslateDist
        for i in range(len(UnPairProb)):
            UnPairProb[i][3] += TranslateDist
        if UseTarget:
            for i in range(len(TargetList)):
                TargetList[i][4] += TranslateDist
            for i in range(len(TargetUnpaired)):
                TargetUnpaired[i][2] += TranslateDist

        vbw = vbh
    else:  # width is greater than height
        TranslateDist = (vbw-vbh)/2.0
        if BoxScaling == 0 and DiscreteScaling == 0:
            for i in range(len(BGData)):
                BGData[i][1] += TranslateDist
        if DiscreteScaling == 0:
            for i in range(len(ColorBar)):
                ColorBar[i][1] += TranslateDist
            for i in range(len(CBLabelData)):
                CBLabelData[i][1] += TranslateDist
        for i in range(len(AxisData)):
            AxisData[i][1] += TranslateDist
        for i in range(len(UnPairedBoxData)):
            UnPairedBoxData[i][1] += TranslateDist
        if Asymmetric == 0 or K == 1:
            for i in range(len(DiagData)):
                DiagData[i][1] += TranslateDist
                DiagData[i][3] += TranslateDist
        for i in range(len(Ticks)):
            Ticks[i][1] += TranslateDist
            Ticks[i][3] += TranslateDist
        if ShowGrid:
            for i in range(len(GridLines)):
                GridLines[i][1] += TranslateDist
                GridLines[i][3] += TranslateDist
        for i in range(len(XTickLabels)):
            XTickLabels[i][1] += TranslateDist
        for i in range(len(YTickLabels)):
            YTickLabels[i][1] += TranslateDist
        if Nstrands > 1:
            for i in range(len(XStrandLabels)):
                XStrandLabels[i][1] += TranslateDist
            for i in range(len(YStrandLabels)):
                YStrandLabels[i][1] += TranslateDist
        if UseTitle:
            YTitle += TranslateDist
        XAxisLabel[1] += TranslateDist
        YAxisLabel[1] += TranslateDist
        RightYAxisLabel[0] += TranslateDist
        if ShowFreeEnergy:
            YFreeEnergy += TranslateDist
        for i in range(len(PairProb)):
            PairProb[i][6] += TranslateDist
        if LowerTriangle == 1 and Asymmetric == 0:
            for i in range(len(PairList)):
                PairList[i][6] += TranslateDist
        for i in range(len(UnPairProb)):
            UnPairProb[i][4] += TranslateDist
        if UseTarget:
            for i in range(len(TargetList)):
                TargetList[i][5] += TranslateDist
            for i in range(len(TargetUnpaired)):
                TargetUnpaired[i][3] += TranslateDist

        vbh = vbw
    SVGWidth = SVGHeight

if IllustratorCompliant:
    SVGWidth = vbw
    SVGHeight = vbh

if RightYAxisLabel[0]>vbw-TitleFontSize:
  RightYAxisLabel[0]=vbw-TitleFontSize
# #############################################################


# ################ SET THE FREE ENERGY FONT SIZE ##############
# This is a hack
FreeEnergyFontSize = FreeEnergyFontSizeFrac*vbh*FreeEnergyFontSizeHack
# #############################################################


# ############### WRITE OUTPUT TO SVG FILE ####################
# Open the output file
f = open(OutputFile,'w')


# Print the header
HeadStr =  '<?xml version="1.0" standalone="no"?>\n'
f.write(HeadStr)
HeadStr = '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"\n'
f.write(HeadStr)
HeadStr = '  "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n'
f.write(HeadStr)
HeadStr = '<svg width="%8.7fpx" height="%8.7fpx" viewBox="0 0 %8.7f %8.7f"\n' % (SVGWidth,SVGHeight,vbw,vbh)
f.write(HeadStr)
HeadStr = '    xmlns="http://www.w3.org/2000/svg" version="1.1" id="dotplot_top">\n'
f.write(HeadStr)
HeadStr = '  <desc>NUPACK DOT PLOT</desc>\n\n'
f.write(HeadStr)


# Write the comments from the input file in the SVG file
if LongHeader:
    CommentStr = '  <!-- ***************************************************************** \n'
    f.write(CommentStr)
    CommentStr = '       The command used to generate this graphic was: \n'
    f.write(CommentStr)
    CommentStr = '       '
    f.write(CommentStr)
    for x in argv:
        if ' ' in x:
            f.write('"%s" ' % x)
        else:
            f.write('%s ' % x)
    f.write('\n\n')
    CommentStr = '       This graphic was generated on %s PST\n\n' % strftime('%c')
    f.write(CommentStr)
    CommentStr = '       Following is the header from the input file %s\n\n' % InputFile
    f.write(CommentStr)
    for x in CommentList:
        f.write('       %s' % x)
    f.write('       ***************************************************************** -->\n\n')
else:
    CommentStr = '  <!-- This graphic was generated by nupack.org on %s PST -->\n\n' % strftime('%c')
    f.write(CommentStr)



# Put in the background color
if BoxScaling == 0 and DiscreteScaling == 0:
    red,green,blue = getRGB(0)
    BGcolor = [red,green,blue]
    CommentStr = '  <!-- Begin Backgroud Coloring --> \n'
    f.write(CommentStr)
    for x in BGData:
        BGStr = '  <rect x="%8.7f" y="%8.7f" ' % (x[0],x[1])
        f.write(BGStr)
        BGStr = 'width="%8.7f" height="%8.7f"\n' % (x[2],x[3])
        f.write(BGStr)
        BGStr = '    fill="rgb(%d,%d,%d)" stroke="none"/>\n' % (BGcolor[0],BGcolor[1],BGcolor[2])
        f.write(BGStr)
    CommentStr = '  <!-- End Backgroud Coloring --> \n\n'
    f.write(CommentStr)


# Draw the colorbar
if DiscreteScaling == 0:
    CommentStr = '  <!-- Begin Colorbar --> \n'
    f.write(CommentStr)
    for x in ColorBar:
        CBStr = '  <rect x="%8.7f" y="%8.7f" ' % (x[0],x[1])
        f.write(CBStr)
        CBStr = 'width="%8.7f" height="%8.7f"\n' % (x[2],x[3])
        f.write(CBStr)
        CBStr = '    fill="rgb(%d,%d,%d)" stroke="none"/>\n' % (x[4][0],x[4][1],x[4][2])
        f.write(CBStr)
    CommentStr = '  <!-- End Colorbar --> \n\n'
    f.write(CommentStr)


# Draw the colorbar labels
if DiscreteScaling == 0:
    CommentStr = '  <!-- Begin Colorbar Labels --> \n'
    f.write(CommentStr)
    for x in CBLabelData:
        CBLabelStr = '  <text x="%8.7f" y="%8.7f" ' % (x[0],x[1])
        f.write(CBLabelStr)
        CBLabelStr = 'font-size="%g" font-family="%s" fill="%s" ' % (ColorBarLabelFontSize,Font,ColorBarLabelColor)
        f.write(CBLabelStr)
        CBLabelStr = 'text-anchor="start">%2.1f</text>\n' % x[2]
        f.write(CBLabelStr)
    CommentStr = '  <!-- End Colorbar Labels --> \n\n'
    f.write(CommentStr)


# Draw the axes
CommentStr = '  <!-- Begin Axes --> \n'
f.write(CommentStr)
for x in AxisData:
    AxesStr = '  <rect x="%8.7f" y="%8.7f" ' % (x[0],x[1])
    f.write(AxesStr)
    AxesStr = 'width="%8.7f" height="%8.7f"\n' % (x[2],x[3])
    f.write(AxesStr)
    AxesStr = '    fill="none" stroke="%s" stroke-width="%8.7f"/>\n' % (AxisColor,AxisWidth)
    f.write(AxesStr)
CommentStr = '  <!-- End Axes --> \n\n'
f.write(CommentStr)


# Draw bounding box for unpaired
CommentStr = '  <!-- Begin Bounding Box for Unpaired --> \n'
f.write(CommentStr)
for x in UnPairedBoxData:
    AxesStr = '  <rect x="%8.7f" y="%8.7f" ' % (x[0],x[1])
    f.write(AxesStr)
    AxesStr = 'width="%8.7f" height="%8.7f"\n' % (x[2],x[3])
    f.write(AxesStr)
    AxesStr = '    fill="none" stroke="%s" stroke-width="%8.7f"/>\n' % (AxisColor,AxisWidth)
    f.write(AxesStr)
CommentStr = '  <!-- End Bounding Box for Unpaired --> \n\n'
f.write(CommentStr)


# Draw the diagonal
if Asymmetric == 0 or K == 1:
    CommentStr = '  <!-- Begin Diagonal Line --> \n'
    f.write(CommentStr)
    for x in DiagData:
        DiagStr = '  <path d="M %8.7f %8.7f L %8.7f %8.7f"\n' % (x[0],x[1],x[2],x[3])
        f.write(DiagStr)
        DiagStr = '    stroke="%s" stroke-width="%8.7f" />\n' % (DiagColor,DiagWidth)
        f.write(DiagStr)
    CommentStr = '  <!-- End Diagonal Line --> \n\n'
    f.write(CommentStr)


# Draw the ticks
CommentStr = '  <!-- Begin Tick Marks --> \n'
f.write(CommentStr)
for x in Ticks:
    TickStr = '  <path d="M %8.7f %8.7f L %8.7f %8.7f"\n' % (x[0],x[1],x[2],x[3])
    f.write(TickStr)
    TickStr = '    stroke="%s" stroke-width="%8.7f" />\n' % (TickColor,TickWidth)
    f.write(TickStr)
CommentStr = '  <!-- End Tick Marks --> \n\n'
f.write(CommentStr)


# Draw the grid lines if desired
if ShowGrid:
    CommentStr = '  <!-- Begin Grid Lines --> \n'
    f.write(CommentStr)
    for x in GridLines:
        GridStr = '  <path d="M %8.7f %8.7f L %8.7f %8.7f"\n' % (x[0],x[1],x[2],x[3])
        f.write(GridStr)
        GridStr = '    stroke="%s" stroke-width="%8.7f"\n' % (GridColor,GridWidth)
        f.write(GridStr)
        GridStr = '    stroke-dasharray="%8.7f,%8.7f" />\n' % (GridDashLength,GridGapLength)
        f.write(GridStr)
    CommentStr = '  <!-- End Grid Lines --> \n\n'
    f.write(CommentStr)


# Draw tick labels
CommentStr = '  <!-- Begin X-Tick Labels --> \n'
f.write(CommentStr)
for x in XTickLabels:
    TickLabelStr = '  <text x="%8.7f" y="%8.7f" ' % (x[0],x[1])
    f.write(TickLabelStr)
    TickLabelStr = 'font-size="%g" font-family="%s" fill="%s" ' % (TickLabelFontSize,Font,TickLabelColor)
    f.write(TickLabelStr)
    TickLabelStr = 'text-anchor="middle">%d</text>\n' % x[2]
    f.write(TickLabelStr)
CommentStr = '  <!-- End X-Tick Labels --> \n\n'
f.write(CommentStr)

CommentStr = '  <!-- Begin Y-Tick Labels --> \n'
f.write(CommentStr)
for x in YTickLabels:
    TickLabelStr = '  <text x="%8.7f" y="%8.7f" ' % (x[0],x[1])
    f.write(TickLabelStr)
    TickLabelStr = 'font-size="%g" font-family="%s" fill="%s" ' % (TickLabelFontSize,Font,TickLabelColor)
    f.write(TickLabelStr)
    TickLabelStr = 'alignment-baseline="text-after-edge" '
    f.write(TickLabelStr)
    TickLabelStr = 'text-anchor="end">%d</text>\n' % x[2]
    f.write(TickLabelStr)
CommentStr = '  <!-- End Y-Tick Labels --> \n\n'
f.write(CommentStr)


# Draw the strand labels
if Nstrands > 1: # we have more than one strand, CHANGE TO >= 1 IF WE ALWAYS WANT LABELS
    CommentStr = '  <!-- Begin X-Strand Labels --> \n'
    f.write(CommentStr)
    for x in XStrandLabels:
        StrandLabelStr = '  <text x="%8.7f" y="%8.7f" ' % (x[0],x[1])
        f.write(StrandLabelStr)
        StrandLabelStr = 'font-size="%g" font-family="%s" fill="%s" ' % (x[3],Font,StrandLabelColor)
        f.write(StrandLabelStr)
        StrandLabelStr = 'text-anchor="middle">%s</text>\n' % x[2]
        f.write(StrandLabelStr)
    CommentStr = '  <!-- End X-Strand Labels --> \n\n'
    f.write(CommentStr)
    CommentStr = '  <!-- Begin Y-Strand Labels --> \n'
    f.write(CommentStr)
    for x in YStrandLabels:
        StrandLabelStr = '  <g transform="rotate(90,%8.7f,%8.7f)">\n' % (x[0],x[1])
        f.write(StrandLabelStr)
        StrandLabelStr = '    <text x="%8.7f" y="%8.7f" ' % (x[0],x[1])
        f.write(StrandLabelStr)
        StrandLabelStr = 'font-size="%g" font-family="%s" fill="%s" ' % (x[3],Font,StrandLabelColor)
        f.write(StrandLabelStr)
        StrandLabelStr = 'alignment-baseline="text-after-edge" '
        f.write(StrandLabelStr)
        StrandLabelStr = 'text-anchor="middle">%s</text>\n' % x[2]
        f.write(StrandLabelStr)
        StrandLabelStr = '  </g>\n'
        f.write(StrandLabelStr)
    CommentStr = '  <!-- End Y-Strand Labels --> \n\n'
    f.write(CommentStr)


# Display the free energy
if ShowFreeEnergy:
    CommentStr = '  <!-- Begin Free Energy --> \n'
    f.write(CommentStr)
    TitleStr = '  <text x="%8.7f" y="%8.7f" ' % (XAxisLabel[0],YFreeEnergy)
    f.write(TitleStr)
    TitleStr = 'font-size="%g" font-family="%s" fill="%s" ' % (FreeEnergyFontSize,Font,FreeEnergyColor)
    f.write(TitleStr)
    TitleStr = 'text-anchor="middle">%s</text>\n' % FreeEnergyString
    f.write(TitleStr)
    CommentStr = '  <!-- End Free Energy --> \n\n'
    f.write(CommentStr)


# Put in the title
if UseTitle:
    CommentStr = '  <!-- Begin Title --> \n'
    f.write(CommentStr)
    TitleStr = '  <text x="%8.7f" y="%8.7f" ' % (XAxisLabel[0],YTitle)
    f.write(TitleStr)
    TitleStr = 'font-size="%g" font-family="%s" fill="%s" ' % (TitleFontSize,Font,TitleColor)
    f.write(TitleStr)
    TitleStr = 'text-anchor="middle">%s</text>\n' % TitleString
    f.write(TitleStr)
    CommentStr = '  <!-- End Title --> \n\n'
    f.write(CommentStr)


# Draw Axis labels
# X-axis
CommentStr = '  <!-- Begin X-Axis Label --> \n'
f.write(CommentStr)
AxisLabelStr = '  <text x="%8.7f" y="%8.7f" ' % (XAxisLabel[0],XAxisLabel[1])
f.write(AxisLabelStr)
AxisLabelStr = 'font-size="%g" font-family="%s" fill="%s" ' % (AxisLabelFontSize,Font,AxisLabelColor)
f.write(AxisLabelStr)
AxisLabelStr = 'text-anchor="middle">Base index</text>\n'
f.write(AxisLabelStr)
CommentStr = '  <!-- End X-Axis Label --> \n\n'
f.write(CommentStr)

# Y-axis
CommentStr = '  <!-- Begin Y-Axis Label --> \n'
f.write(CommentStr)
AxisLabelStr = '  <g transform="rotate(-90,%8.7f,%8.7f)">\n' % (YAxisLabel[0],YAxisLabel[1])
f.write(AxisLabelStr)
AxisLabelStr = '    <text x="%8.7f" y="%8.7f" ' % (YAxisLabel[0],YAxisLabel[1])
f.write(AxisLabelStr)
AxisLabelStr = 'font-size="%g" font-family="%s" fill="%s" ' % (AxisLabelFontSize,Font,AxisLabelColor)
f.write(AxisLabelStr)
AxisLabelStr = 'alignment-baseline="text-after-edge" '
f.write(AxisLabelStr)
AxisLabelStr = 'text-anchor="middle">Base index</text>\n'
f.write(AxisLabelStr)
AxisLabelStr = '  </g>\n'
f.write(AxisLabelStr)
CommentStr = '  <!-- End Y-Axis Label --> \n\n'
f.write(CommentStr)

# Right Y-axis
CommentStr = '  <!-- Begin Y-Axis Label --> \n'
f.write(CommentStr)
AxisLabelStr = '  <g transform="rotate(90,%8.7f,%8.7f)">\n' % (RightYAxisLabel[0],RightYAxisLabel[1])
f.write(AxisLabelStr)
AxisLabelStr = '    <text x="%8.7f" y="%8.7f" ' % (RightYAxisLabel[0],RightYAxisLabel[1])
f.write(AxisLabelStr)
AxisLabelStr = 'font-size="%g" font-family="%s" fill="%s" ' % (AxisLabelFontSize,Font,AxisLabelColor)
f.write(AxisLabelStr)
AxisLabelStr = 'alignment-baseline="text-after-edge" '
f.write(AxisLabelStr)
AxisLabelStr = 'text-anchor="middle">%s</text>\n' % (LegendLabel)
f.write(AxisLabelStr)
AxisLabelStr = '  </g>\n'
f.write(AxisLabelStr)
CommentStr = '  <!-- End Y-Axis Label --> \n\n'
f.write(CommentStr)



# Draw the squares (dots)
CommentStr = '  <!-- Begin BP Prob Dots --> \n'
f.write(CommentStr)
# Old code with transparency
if noX:
    for x in PairProb:
        if x[2] >= minP[2]:
            DotStr = '  <rect x="%8.7f" y="%8.7f" ' % (x[5],x[6])
            f.write(DotStr)
            DotStr = 'width="%8.7f" height="%8.7f"\n' % (x[7],x[7])
            f.write(DotStr)
            DotStr = '    fill="rgb(%d,%d,%d)" stroke="none"' % (x[8][0],x[8][1],x[8][2])
            f.write(DotStr)
#            if x[-1] == 1:
#                DotStr = ' opacity="%.2f"' % BaseOpacityOverTarget
#                f.write(DotStr)
            DotStr = '/>\n'
            f.write(DotStr)
else :
    for x in PairProb:
        if x[2] >= minP[2]:
            DotStr = '  <rect x="%8.7f" y="%8.7f" ' % (x[5],x[6])
            f.write(DotStr)
            DotStr = 'width="%8.7f" height="%8.7f"\n' % (x[7],x[7])
            f.write(DotStr)
            DotStr = '    fill="rgb(%d,%d,%d)" stroke="none"/>\n' % (x[8][0],x[8][1],x[8][2])
            f.write(DotStr)
CommentStr = '  <!-- End BP Prob Dots --> \n\n'
f.write(CommentStr)


# Draw the squares in the lower triangle
if LowerTriangle == 1 and Asymmetric == 0:
    CommentStr = '  <!-- Begin Lower Triangle Dots --> \n'
    f.write(CommentStr)
    for x in PairList:
        if x[2] >= minP[2]:
            DotStr = '  <rect x="%8.7f" y="%8.7f" ' % (x[5],x[6])
            f.write(DotStr)
            DotStr = 'width="%8.7f" height="%8.7f"\n' % (x[7],x[7])
            f.write(DotStr)
            DotStr = '    fill="rgb(%d,%d,%d)" stroke="none"/>\n' % (x[8][0],x[8][1],x[8][2])
            f.write(DotStr)
    CommentStr = '  <!-- End Lower Triangle Dots --> \n\n'
    f.write(CommentStr)


# Draw the squares for unpaired probabilities
CommentStr = '  <!-- Begin Unpaired Prob Dots --> \n'
f.write(CommentStr)
for x in UnPairProb:
    if x[1] >= minP[2]:
        DotStr = '  <rect x="%8.7f" y="%8.7f" ' % (x[3],x[4])
        f.write(DotStr)
        DotStr = 'width="%8.7f" height="%8.7f"\n' % (x[5],x[5])
        f.write(DotStr)
        DotStr = '    fill="rgb(%d,%d,%d)" stroke="none"/>\n' % (x[6][0],x[6][1],x[6][2])
        f.write(DotStr)
CommentStr = '  <!-- End Unpaired Prob Dots --> \n\n'
f.write(CommentStr)


# Draw squares for target
if UseTarget:
    corr = TargetLineWidth/2.0  # Correction to position
    CommentStr = '  <!-- Begin Target Pairs --> \n'
    f.write(CommentStr)

    if noX:
        circleRadius = BaseWidth * sqrt(2) / 2.0 - TargetLineWidth/2.0
        for x in TargetList:
            #            TargetStr = '  <rect x="%8.7f" y="%8.7f" ' % (x[4],x[5])
            #            f.write(TargetStr)
            #            TargetStr = 'width="%8.7f" height="%8.7f"\n' % (BaseWidth,BaseWidth)
            #            f.write(TargetStr)
            #            TargetStr = '    fill="rgb(%d,%d,%d)" stroke="none"/>\n' % (TargetColor[0],TargetColor[1],TargetColor[2])
            #            f.write(TargetStr)
            TargStr = '  <circle cx="%8.7f" cy="%8.7f" r = "%8.7f"\n' % (x[4]+BaseWidth/2.0,x[5]+BaseWidth/2.0,circleRadius)
            f.write(TargStr)
            TargStr = '    fill="none" stroke="rgb(%d,%d,%d)" stroke-width="%8.7f"/>\n' % (TargetColor[0],TargetColor[1],TargetColor[2],TargetLineWidth)
            f.write(TargStr)
    else:
        for x in TargetList:
            TargStr = '  <path d="M %8.7f %8.7f\n' % (x[4]+corr,x[5]+corr)
            f.write(TargStr)
            TargStr = '           L %8.7f %8.7f\n' % (x[4]+BaseWidth-corr,x[5]+corr)
            f.write(TargStr)
            TargStr = '           L %8.7f %8.7f\n' % (x[4]+BaseWidth-corr,x[5]+BaseWidth-corr)
            f.write(TargStr)
            TargStr = '           L %8.7f %8.7f\n' % (x[4]+corr,x[5]+BaseWidth-corr)
            f.write(TargStr)
            TargStr = '           L %8.7f %8.7f\n' % (x[4]+corr,x[5]+corr)
            f.write(TargStr)
            TargStr = '           M %8.7f %8.7f\n' % (x[4]+corr,x[5]+corr)
            f.write(TargStr)
            TargStr = '           L %8.7f %8.7f\n' % (x[4]+BaseWidth-corr,x[5]+BaseWidth-corr)
            f.write(TargStr)
            TargStr = '           M %8.7f %8.7f\n' % (x[4]+BaseWidth-corr,x[5]+corr)
            f.write(TargStr)
            TargStr = '           L %8.7f %8.7f"\n' % (x[4]+corr,x[5]+BaseWidth-corr)
            f.write(TargStr)
            TargStr = '     stroke="rgb(%d,%d,%d)" fill = "none"\n' % (TargetColor[0],TargetColor[1],TargetColor[2])
            f.write(TargStr)
            TargStr = '     stroke-width="%8.7f"/>\n' % TargetLineWidth
            f.write(TargStr)
    CommentStr = '  <!-- End Target Pairs --> \n\n'
    f.write(CommentStr)


# Draw unpaired squares for target
if UseTarget:
    corr = TargetLineWidth/2.0  # Correction to position
    CommentStr = '  <!-- Begin Target Unpaired --> \n'
    f.write(CommentStr)
    if noX:
        circleRadius = BaseWidth /2.0
        for x in TargetUnpaired:
            TargStr = '  <circle cx="%8.7f" cy="%8.7f" r = "%8.7f"\n' % (x[2]+BaseWidth/2.0,x[3]+BaseWidth/2.0,circleRadius)
            f.write(TargStr)
            TargStr = '    fill="none" stroke="rgb(%d,%d,%d)" stroke-width="%8.7f"/>\n' % (TargetColor[0],TargetColor[1],TargetColor[2],TargetLineWidth)
            f.write(TargStr)
    else:
        for x in TargetUnpaired:
            TargStr = '  <path d="M %8.7f %8.7f\n' % (x[2]+corr,x[3]+corr)
            f.write(TargStr)
            TargStr = '           L %8.7f %8.7f\n' % (x[2]+BaseWidth-corr,x[3]+corr)
            f.write(TargStr)
            TargStr = '           L %8.7f %8.7f\n' % (x[2]+BaseWidth-corr,x[3]+BaseWidth-corr)
            f.write(TargStr)
            TargStr = '           L %8.7f %8.7f\n' % (x[2]+corr,x[3]+BaseWidth-corr)
            f.write(TargStr)
            TargStr = '           L %8.7f %8.7f\n' % (x[2]+corr,x[3]+corr)
            f.write(TargStr)
            TargStr = '           M %8.7f %8.7f\n' % (x[2]+corr,x[3]+corr)
            f.write(TargStr)
            TargStr = '           L %8.7f %8.7f\n' % (x[2]+BaseWidth-corr,x[3]+BaseWidth-corr)
            f.write(TargStr)
            TargStr = '           M %8.7f %8.7f\n' % (x[2]+BaseWidth-corr,x[3]+corr)
            f.write(TargStr)
            TargStr = '           L %8.7f %8.7f"\n' % (x[2]+corr,x[3]+BaseWidth-corr)
            f.write(TargStr)
            TargStr = '     stroke="rgb(%d,%d,%d)" fill = "none"\n' % (TargetColor[0],TargetColor[1],TargetColor[2])
            f.write(TargStr)
            TargStr = '     stroke-width="%8.7f"/>\n' % TargetLineWidth
            f.write(TargStr)
    CommentStr = '  <!-- End Target Unpaired --> \n\n'
    f.write(CommentStr)


# End the SVG file
EndStr = '</svg>\n'
f.write(EndStr)


# Close output
f.close()
###############################################################


# Exit successfully
exit(0)

