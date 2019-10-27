#! /usr/bin/python

# Python script to take input from command line and process .ppairs and
# files to make x-y data for plotting a melt curve.
#
# For the example of a melt (output from with a melt range from 10 to 30
# deg. C in 10 degree increments for three different strands, the input
# files are:
#   ./10.0/prefix_0.ppairs
#   ./10.0/prefix_1.ppairs
#   ./10.0/prefix_2.ppairs
#   ./20.0/prefix_0.ppairs
#   ./20.0/prefix_1.ppairs
#   ./20.0/prefix_2.ppairs
#   ./30.0/prefix_0.ppairs
#   ./30.0/prefix_1.ppairs
#   ./30.0/prefix_2.ppairs
#
# The usage is:
#  MeltParserMany.py prefix nstrands tempStart tempIncrement tempEnd outputfile
#
#  prefix specifies the name of the files, e.g. prefix_0.ppairs.
#    The files must be of the format as outputted
#    by the programs pairs and concentrations, respectively.
#  tempStart is the starting temperature for the melt profile
#  tempIncrement is the temperature increment for the melt profile
#  tempEnd is the ending temperature
#  outputfile is the name of the output file to which the data are
#   written.  This file will in turn be read by MeltPlot.Py to generate
#   the plot.
#
# January, 2007, Justin Bois

# Import the modules we'll use
import sys
argv = sys.argv
exit = sys.exit
import re
WhiteSpaceSearch = re.compile('\s+')


# ######### GET DATA FROM COMMAND LINE INPUT ##################
prefix = argv[1]
nstrands = int(argv[2])
tempStart = float(argv[3])
tempIncrement = float(argv[4])
tempEnd = float(argv[5])
outputFile = argv[6]

Tlist = []
temp = tempStart
while temp <= tempEnd + 0.0000001: # Add small number to avoid rounding error
    Tlist.append(temp)
    temp += tempIncrement

fUnpaired = [0.0]*len(Tlist)
nBases = [0]*nstrands  # Number of bases in each strand

for i in range(len(Tlist)):
    for j in range(nstrands):
        inputFile = './%.1f/%s_%d.ppairs' % (Tlist[i],prefix,j)
        f = open(inputFile)
        
        # Blow through comment lines until we get to the total number of bases
        line = f.readline()
        while not(line[0].isdigit()) and line != '':
            line = f.readline()

        if line == '':
            print 'Error: No records found in input file %s.' % inputFile
            exit(1)

        # The line we're at now is the total number of bases
        LineList = WhiteSpaceSearch.split(line)
        nBases[j] = int(LineList[0]);

        # Now read records and find unpaired bases
        line = f.readline()
        while line != '' and line[0].isdigit():
            LineList = WhiteSpaceSearch.split(line)
            if int(LineList[1]) == nBases[j] + 1:
                fUnpaired[i] += float(LineList[2])

            line = f.readline()
        f.close()

        nBasesTotal = sum(nBases)

    fUnpaired[i] /= float(nBasesTotal)

print nBasesTotal

# Make input file
f = open(outputFile,'w')
for i in range(len(Tlist)):
    f.write('%.1f\t%8.7f\n' % (Tlist[i],fUnpaired[i]))
f.close()


exit(0)

