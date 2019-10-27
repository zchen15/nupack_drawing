#! /usr/bin/python

# Python script for generating SVG plots for melt curves
#
# Usage: MeltPlot.py [-datafile dataouputfile] [-title titleString] [-white] inputfile outputfile
#
#  input file must contain only two columns of data.  The left is the
#  temperature and the right is the fraction of unpaired bases at
#  equilibrium.  The right column must have values between 0 and 1.
#  outputfile is the name of the SVG file to which the plot is
#  written.  If -datafile is chosen, the raw data for the melt curve
#  is written to the file dataoutputfile.
#  -white flag makes all axes and lines white

# Import a modules we'll use
import sys
argv = sys.argv
exit = sys.exit
import re
WhiteSpaceSearch = re.compile('\s+')
import time
strftime = time.strftime
import math
floor = math.floor


# ######## CONSTANTS DESCRIBING DIMENSIONS, ETC OF PLOT ######
Font = 'sans-serif'

# Delimiters for data plot
Delimiters = [',',' ','\t']
DelimitIndex = 2 # = 0 for comma, 1 for space, and 2 for tab

# Whether or not to have grid lines
GridLines = 1

# Works well when there is a colormap
WidthToHeightRatio = 1.17
SVGHeight = 500 # Height in pixels.  This works for all plots.
MakeSquare = 1 # This makes the dimensions of the SVG graphic square
IllustratorCompliant = 1 # This makes the plot compliant for importing into Illustrator
# vbw and SVGWidth are set off the height and such that everything fits
YLength = 800.0 # Height of the plot
YXRatio = 1.2 # Ratio of y to x axis lengths
TopMargin = 30.0 # Distance from top of document to top of plot
BottomMargin = 30.0 # Buffer on bottom of plot
LeftMargin = 136 # Distance from left of document to left edge of plot (excludes axis labels)
RightMargin = 30 # Buffer on right size of document
StandardNTicks = 10 # Target number of ticks (actual can vary around this value)
OverhangFrac = 0.05 # Fraction of axis that goes beyond last tick
fUnpairedMinCutoff = 0.0 # Minimal value of fUnpaired before we make 0.0 the minimal f value
fUnpairedMaxCutoff = 1.0 # Maximal value of fUnpaired before we make 1.0 the maximal f value
Enforce0to1 = 1 # = 1 if we force the y-axis to go from zero to 1
    
AxisWidth = 2.0 # Width of lines of axes
TickWidth = 1.0 # Width of tick lines
TickLengthFrac = 0.015 # Fraction of the length of axes that ticks are

LastTick = 1  # Whether or not to include a tick at the strand length if within LastTickFrac
LastTickFrac = 0.5  # The maximal fraction toward the end of the strand that the last tick will be
AxisLabelFontSize = 30  # Font size of the axis labels
XAxisLabelOffset = 12  # How much the axis label is offset from the tick labels
YAxisLabelOffset = 36  # How much the axis label is offset from the tick labels
TitleOffset = 36 # How much to offset the title of the plot
TitleFontSize = 36  # Font size for title
TickLabelFontSize = 24  # Font size for tick labels
TickLabelOffset = 6  # How much the tick labels are offset from the axis
GridWidth = 1

DotRadiusFrac = 0.01 # Dot radius as fraction of size of plot

# Define colors (this changes if -white is chosen)
DataPointColor = 'black'  # Color of data points
AxisColor = 'black' # Color of axes
TickColor = 'black' # Color of ticks
TitleColor = 'black'  # Font color for title
AxisLabelColor = 'black'  # Font color for axis labels
TickLabelColor = 'black'  # Font color for tick labels
GridColor = 'gray'

# Exit codes
INVALID_INPUT_FILE_ERROR = 2
FAILED_TO_FIND_PERMUTATION_ERROR = 3
FAILED_TO_FIND_INITIAL_CONCENTRATIONS = 4

# ##### Fudge factors for vertical alignment of text
# These should vanish when we finally figure out how to align these things
# This works well for TickLabelFontSize = 24
TickHorizFudgeFact = 2.4/3.0*TickLabelFontSize
TickVertFudgeFact = 1.0/3.0*TickLabelFontSize

# These work for AxisLabelFontSize = 36
AxisHorizFudgeFact = 2.1/3.0*AxisLabelFontSize
AxisVertFudgeFact = 0.0/3.0*AxisLabelFontSize
#############################################################


# ######### GET DATA FROM COMMAND LINE INPUT ##################
InputFile = argv[-2]
OutputFile = argv[-1]

# Find whether or not we do a data file
DataFileIndex = -1  # Initialize in case the option isn't chosen
if '-datafile' in argv:
   DataFileIndex  = argv.index('-datafile') + 1
elif '--datafile' in argv:
    DataFileIndex = argv.index('--datafile') + 1
elif '-d' in argv:
    DataFileIndex = argv.index('-d') + 1

if DataFileIndex == -1:
    WriteDataFile = 0
else:
    WriteDataFile = 1
    DataFile = argv[DataFileIndex]

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

AllWhite = 0
if '-white' in argv:
    AllWhite = 1
elif '--white' in argv:
    AllWhite = 1
elif '-w' in argv:
    AllWhite = 1
if AllWhite:
   DataPointColor = 'white'  # Color of data points
   AxisColor = 'white' # Color of axes
   TickColor = 'white' # Color of ticks
   TitleColor = 'white'  # Font color for title
   AxisLabelColor = 'white'  # Font color for axis labels
   TickLabelColor = 'white'  # Font color for tick labels
   GridColor = 'white'
###############################################################


# ################### READ INPUT DATA #########################
T = []
fUnpaired = []

# Open the input file for reading
f = open(InputFile,'r')
    
# Read pertinent lines that have data
line = f.readline()
while line != '':
    if line[0] != '%' and line[0] != '\n':
        LineData = WhiteSpaceSearch.split(line)
        T.append(float(LineData[0]))
        fUnpaired.append(float(LineData[1]))
    line = f.readline()
f.close()
# #############################################################


nPoints = len(T)


# ############### WRITE DATA FILE #############################
if WriteDataFile:
    f = open(DataFile,'w')
    f.write('% This is a data file of a melt curve generated by nupack.org\n')
    f.write('%% on %s PST\n' % strftime('%c'))
    f.write('% The left column contains the temperatures in degrees Celsius.\n')
    f.write('% The right column contains the corresponding fractions of bases that\n')
    f.write('% are unpaired at equilibrium.\n')
    for i in range(nPoints):
        f.write('%.1f%s%8.7f\n' % (T[i],Delimiters[DelimitIndex],fUnpaired[i]))
    f.close()
# #############################################################


# Set up the top margin depending on the title
if UseTitle:
   TopMargin += TitleOffset + TitleFontSize

# The coordinates for the corner of the box (the axes) containing the plot
BoxCornerX = LeftMargin
BoxCornerY = TopMargin
XLength = YXRatio*YLength
DotRadius = DotRadiusFrac*min(XLength,YLength)

# ############### SET UP AXES FOR PLOTTING ####################
# AxisData has [upper left x-coord, upper left corner y-coord, x-dimension, y-dimension]
# for the axes in the plot
AxisData = []
CornerY = BoxCornerY - AxisWidth/2.0
CornerX = BoxCornerX - AxisWidth/2.0
BoxDimX = XLength + AxisWidth
BoxDimY = YLength + AxisWidth
AxisData = [CornerX,CornerY,BoxDimX,BoxDimY]
###############################################################


# ################ SET UP TICKS FOR PLOT ######################
# Min and max values for the plot
fUnpairedMin = min(fUnpaired)
fUnpairedMax = max(fUnpaired)
Tmin = min(T)
Tmax = max(T)
Tinterval = (Tmax - Tmin) / (len(T)-1.0)

# Get beginning and end values for plot
T0 = Tmin - OverhangFrac*(Tmax-Tmin)
Tf = Tmax + OverhangFrac*(Tmax-Tmin)

if Enforce0to1:
    fUnpaired0 = 0.0
    fUnpairedf = 1.0
else:    
    if fUnpairedMin <= fUnpairedMinCutoff:
        fUnpaired0 = 0.0
    else:
        fUnpaired0 = fUnpairedMin - OverhangFrac*(fUnpairedMax - fUnpairedMin)

    if fUnpairedMax >= fUnpairedMaxCutoff:
        fUnpairedf = 1.0
    else:
        fUnpairedf = fUnpairedMax + OverhangFrac*(fUnpairedMax - fUnpairedMax)

# The tick intervals
XInterval = (floor(nPoints/StandardNTicks) + 1) * Tinterval
YInterval = (fUnpairedMax - fUnpairedMin) / (StandardNTicks - 1.0)


# Each entry in Ticks has [x1,y1,x2,y2], the coords for tick lines
XTickLabels = []
YTickLabels = []


Ticks = []
Grid = []
# First do vertical grid lines
Temp = Tmin
x1 = BoxCornerX + (Temp - T0)/(Tf - T0)*XLength
yLabel = BoxCornerY + YLength + AxisWidth + TickLabelOffset + TickHorizFudgeFact
ybottom1 = BoxCornerY + YLength
ybottom2 = ybottom1 - TickLengthFrac*min(XLength,YLength)
ytop1 = BoxCornerY
ytop2 = BoxCornerY + TickLengthFrac*min(XLength,YLength)
while Temp < Tf:
    Ticks.append([x1,ybottom1,x1,ybottom2])
    Ticks.append([x1,ytop1,x1,ytop2])
    Grid.append([x1,BoxCornerY+YLength,x1,BoxCornerY])
    XTickLabels.append([x1,yLabel,Temp])
    Temp += XInterval
    x1 = BoxCornerX + (Temp - T0)/(Tf - T0)*XLength


# Do vertical ticks
if Enforce0to1:
    fUp = 0.0
    y1 = BoxCornerY + (fUnpairedf - fUp)/(fUnpairedf - fUnpaired0)*YLength
    xLabel = BoxCornerX - AxisWidth - TickLabelOffset
    xleft1 = BoxCornerX
    xleft2 = xleft1 + TickLengthFrac*min(XLength,YLength)
    xright1 = BoxCornerX + XLength
    xright2 = xright1 - TickLengthFrac*min(XLength,YLength)
    while fUp <= 1.0:
        if fUp >= 0.0001 and fUp <= 0.999: # No ticks at 0 and 1
            Ticks.append([xleft1,y1,xleft2,y1])
            Ticks.append([xright1,y1,xright2,y1])
            Grid.append([BoxCornerX,y1,BoxCornerX+XLength,y1])
        YTickLabels.append([xLabel,y1+TickVertFudgeFact,fUp])
        fUp += 0.1
        y1 = BoxCornerY + (fUnpairedf - fUp)/(fUnpairedf - fUnpaired0)*YLength
else:
    fUp = fUnpairedMin
    y1 = BoxCornerY + (fUnpairedf - fUp)/(fUnpairedf - fUnpaired0)*YLength + TickVertFudgeFact
    xLabel = BoxCornerX - AxisWidth - TickLabelOffset
    xleft1 = BoxCornerX
    xleft2 = xleft1 + TickLengthFrac*min(XLength,YLength)
    xright1 = BoxCornerX + XLength
    xright2 = xright1 - TickLengthFrac*min(XLength,YLength)
    while fUp < fUnpairedf:
        Ticks.append([xleft1,y1,xleft2,y1])
        Ticks.append([xright1,y1,xright2,y1])
        Grid.append([BoxCornerX,y1,BoxCornerX+XLength,y1])
        YTickLabels.append([xLabel,y1,fUp])
        fUp += YInterval
        y1 = BoxCornerY + (fUnpairedf - fUp)/(fUnpairedf - fUnpaired0)*YLength
# #############################################################


# ################ SET UP DATA FOR PLOT #######################
Data = []
for i in range(nPoints):
    x = BoxCornerX + (T[i] - T0)/(Tf - T0)*XLength
    y = BoxCornerY + (fUnpairedf - fUnpaired[i])/(fUnpairedf - fUnpaired0)*YLength
    Data.append([x,y])
# #############################################################


# ################# SET UP AXIS LABELS ########################
# Horizontal Label
x1 = BoxCornerX + XLength/2.0
y1 = BoxCornerY + YLength + AxisWidth + TickLabelOffset + TickHorizFudgeFact + TickLabelFontSize \
     + XAxisLabelOffset + AxisHorizFudgeFact
XAxisLabel = [x1,y1]


# Vertical Label
x1 = BoxCornerX - AxisWidth - TickLabelOffset - TickHorizFudgeFact - TickLabelFontSize \
     - YAxisLabelOffset
y1 = BoxCornerY + YLength/2.0
YAxisLabel=[x1,y1]
# #############################################################


# ##################### SET UP TITLE ##########################
if UseTitle:
    XTitle = BoxCornerX + BoxDimX/2.0;
    YTitle = BoxCornerY - AxisWidth
    YTitle -= TitleOffset

    # Fix the title string to have degree symbols (replace degC and degF)
    TitleString = TitleString.replace('degC','&#x2103;')
# #############################################################


# ########### FIND DIMENSIONS OF VIEWBOX AND DOCUMENT #########
vbh = BoxCornerY + YLength + AxisWidth + TickLabelOffset + TickHorizFudgeFact + TickLabelFontSize \
     + XAxisLabelOffset + AxisHorizFudgeFact + BottomMargin
vbw = BoxCornerX + XLength + AxisWidth + RightMargin
if MakeSquare:
    if vbh > vbw:
        TranslateDist = (vbh-vbw)/2.0
        AxisData[0] += TranslateDist
        if GridLines:
            for i in range(len(Grid)):
                Grid[i][0] += TranslateDist
                Grid[i][2] += TranslateDist
        for i in range(len(Ticks)):
            Ticks[i][0] += TranslateDist
            Ticks[i][2] += TranslateDist
        for i in range(len(XTickLabels)):
            XTickLabels[i][0] += TranslateDist
        for i in range(len(YTickLabels)):
            YTickLabels[i][0] += TranslateDist
        XAxisLabel[0] += TranslateDist
        YAxisLabel[0] += TranslateDist
        if UseTitle:
           XTitle += TranslateDist
        for i in range(len(Data)):
            Data[i][0] += TranslateDist
        vbw = vbh
    else: # vbh < vbw
        TranslateDist = (vbw-vbh)/2.0
        AxisData[1] += TranslateDist
        if GridLines:
            for i in range(len(Grid)):
                Grid[i][1] += TranslateDist
                Grid[i][3] += TranslateDist
        for i in range(len(Ticks)):
            Ticks[i][1] += TranslateDist
            Ticks[i][3] += TranslateDist
        for i in range(len(XTickLabels)):
            XTickLabels[i][1] += TranslateDist
        for i in range(len(YTickLabels)):
            YTickLabels[i][1] += TranslateDist
        XAxisLabel[1] += TranslateDist
        YAxisLabel[1] += TranslateDist
        if UseTitle:
           YTitle += TranslateDist
        for i in range(len(Data)):
            Data[i][1] += TranslateDist
        vbh = vbw

if IllustratorCompliant:
    SVGWidth = vbw
    SVGHeight = vbh
else:
    SVGWidth = float(vbw)/vbh*SVGHeight
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
HeadStr = '    xmlns="http://www.w3.org/2000/svg" version="1.1">\n'
f.write(HeadStr)
HeadStr = '  <desc>NUPACK MELT CURVE</desc>\n\n'
f.write(HeadStr)
CommentStr = '  <!-- This graphic was generated by nupack.org on %s PST -->\n\n' % strftime('%c')
f.write(CommentStr)

# Draw the axes
CommentStr = '  <!-- Begin Axes --> \n'
f.write(CommentStr)
AxesStr = '  <rect x="%8.7f" y="%8.7f" ' % (AxisData[0],AxisData[1])
f.write(AxesStr)
AxesStr = 'width="%8.7f" height="%8.7f"\n' % (AxisData[2],AxisData[3])
f.write(AxesStr)
AxesStr = '    fill="none" stroke="%s" stroke-width="%8.7f"/>\n' % (AxisColor,AxisWidth)
f.write(AxesStr)
CommentStr = '  <!-- End Axes --> \n\n'
f.write(CommentStr)


# Draw the grid lines
if GridLines == 1:
    CommentStr = '  <!-- Begin Grid Lines --> \n'
    f.write(CommentStr)
    for x in Grid:
        GridStr = '  <path d="M %8.7f %8.7f L %8.7f %8.7f"\n' % (x[0],x[1],x[2],x[3])
        f.write(GridStr)
        GridStr = '    stroke="%s" stroke-width="%8.7f" \n' % (GridColor,GridWidth)
        f.write(GridStr)
        GridStr = '    stroke-dasharray="2,4" />\n'
        f.write(GridStr)
    CommentStr = '  <!-- End Grid Lines --> \n\n'
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


# Draw tick labels
CommentStr = '  <!-- Begin X-Tick Labels --> \n'
f.write(CommentStr)
for x in XTickLabels:
    TickLabelStr = '  <text x="%8.7f" y="%8.7f" ' % (x[0],x[1])
    f.write(TickLabelStr)
    TickLabelStr = 'font-size="%g" font-family="%s" fill="%s" ' % (TickLabelFontSize,Font,TickLabelColor)
    f.write(TickLabelStr)
    TickLabelStr = 'text-anchor="middle">%.1f</text>\n' % x[2]
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
    TickLabelStr = 'text-anchor="end">%.2f</text>\n' % x[2]
    f.write(TickLabelStr)
CommentStr = '  <!-- End Y-Tick Labels --> \n\n'
f.write(CommentStr)

# Put in the title
if UseTitle:
    CommentStr = '  <!-- Begin Title --> \n'
    f.write(CommentStr)
    TitleStr = '  <text x="%8.7f" y="%8.7f" ' % (XTitle,YTitle)
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
AxisLabelStr = 'text-anchor="middle">Temperature (C)</text>\n'
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
AxisLabelStr = 'text-anchor="middle">Fraction of bases unpaired</text>\n'
f.write(AxisLabelStr)
AxisLabelStr = '  </g>\n'
f.write(AxisLabelStr)
CommentStr = '  <!-- End Y-Axis Label --> \n\n'
f.write(CommentStr)

# Draw data points
CommentStr = '  <!-- Begin X-Y Data --> \n'
f.write(CommentStr)
for x in Data:
    DataStr = '  <circle cx="%8.7f" cy="%8.7f" r="%8.7f"\n' % (x[0],x[1],DotRadius)
    f.write(DataStr)
    DataStr = '        fill="%s" stroke="none" />' % DataPointColor
    f.write(DataStr)
CommentStr = '  <!-- End X-Y Data --> \n'
f.write(CommentStr)

# End the SVG file
EndStr = '</svg>'
f.write(EndStr)


# Close output
f.close()
###############################################################


# Exit successfully
exit(0)

