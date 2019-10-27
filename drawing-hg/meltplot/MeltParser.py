#! /usr/bin/python

# Python script to take input from command line and process .pairs and
# .fpairs files to make x-y data for plotting a melt curve.
#
# For the example of a multi stranded melt (output from
# concentrations) with a melt range from 10 to 40 deg. C in 10 degree
# increments, the input files are:
#   ./10.0/prefix.fpairs
#   ./20.0/prefix.fpairs
#   ./30.0/prefix.fpairs
#   ./40.0/prefix.fpairs
#
# The usage is:
#  MeltParser.py filename tempStart tempIncrement tempEnd outputfile
#
#  filename is the name of the files, e.g. prefix.ppairs or
#    prefix.fpairs.  The files must be of the format as outputted
#    by the programs pairs and concentrations, respectively.
#  tempStart is the starting temperature for the melt profile
#  tempIncrement is the temperature increment for the melt profile
#  tempEnd is the ending temperature
#  outputfile is the name of the output file to which the data are
#   written.  This file will in turn be read by MeltPlot.Py to generate
#    the plot.
#
# January, 2007, Justin Bois

# Import the modules we'll use
import sys
argv = sys.argv
exit = sys.exit
import re
WhiteSpaceSearch = re.compile('\s+')


# ######### GET DATA FROM COMMAND LINE INPUT ##################
filename = argv[1]
tempStart = float(argv[2])
tempIncrement = float(argv[3])
tempEnd = float(argv[4])
outputFile = argv[5]

if filename[-6:] == 'ppairs':
    ppairsMode = 1
else:
    ppairsMode = 0

Tlist = []
temp = tempStart
while temp <= tempEnd + 0.0000001: # Add small number to avoid rounding error
    Tlist.append(temp)
    temp += tempIncrement

fUnpaired = [0.0]*len(Tlist)

if ppairsMode:
    for k in range(len(Tlist)):
        inputFile = './%.1f/%s' % (Tlist[k],filename)
        f = open(inputFile)

        # Now blow through lines to number of bases
        line = f.readline()
        while not(line[0].isdigit()) and line != '':
            line = f.readline()
        if line == '':
            print 'Error: No records found in input file %s.' % inputFile
            exit(1)

        # The line we're at now is the total number of bases
        LineList = WhiteSpaceSearch.split(line)
        nBases = int(LineList[0])

        # Now read records and find unpaired bases
        line = f.readline()
        while line != '' and line[0].isdigit():
            LineList = WhiteSpaceSearch.split(line)
            if int(LineList[1]) == nBases + 1:
                fUnpaired[k] += float(LineList[2])
            line = f.readline()
        f.close()

        fUnpaired[k] /= nBases

else: # fpairs file
    for k in range(len(Tlist)):
        inputFile = './%.1f/%s' % (Tlist[k],filename)
        f = open(inputFile)

        # Blow through comment lines until we get to inital monomer concentrations
        line = f.readline()
        while line[0:17] != '% Initial monomer' and line != '':
            line = f.readline()
        if line == '':
            print 'Error: No initial monomer concentrations found in %s.' % inputFile
            exit(1)

        line = f.readline()
        # We are now at the first initial monomer concentration
        x0 = []
        while (len(line) > 10 and line[4].isdigit()):
            LineList = WhiteSpaceSearch.split(line)
            x0.append(float(LineList[2]))
            line = f.readline()

        # We have passed the initial monomer concentrations.  Now read to sequences
        Nstrands = len(x0)
        while line[0:13] != '% id sequence' and line != '':
            line = f.readline()
        if line == '':
            print 'Error: No sequences found in %s.' % inputFile
            exit(1)

        line = f.readline()
        # We are now at the first sequence
        sequences = []
        for i in range(Nstrands):
            LineList = WhiteSpaceSearch.split(line)
            sequences.append(LineList[2])
            line = f.readline()

        # Now blow through lines to number of bases
        while not(line[0].isdigit()) and line != '':
            line = f.readline()
        if line == '':
            print 'Error: No records found in input file %s.' % inputFile
            exit(1)

        # The line we're at now is the total number of bases
        LineList = WhiteSpaceSearch.split(line)
        nBases = int(LineList[0])

        # Set up variables
        fup = [0.0]*Nstrands
        StrandKey = []
        for i in range(Nstrands):
            for j in range(len(sequences[i])):
                StrandKey.append(i)

        # Now read records and find unpaired bases
        line = f.readline()
        while line != '' and line[0].isdigit():
            LineList = WhiteSpaceSearch.split(line)
            if int(LineList[1]) == nBases + 1:
                fup[StrandKey[int(LineList[0])-1]] += float(LineList[2])
            line = f.readline()
        f.close()

        for i in range(Nstrands):
            # fup[i] /= len(sequences[i])
            fup[i] *= x0[i]
            x0[i] *= len(sequences[i])

        fUnpaired[k] = sum(fup)/sum(x0)



# Make input file
f = open(outputFile,'w')
for i in range(len(Tlist)):
    f.write('%.1f\t%8.7f\n' % (Tlist[i],fUnpaired[i]))
f.close()


exit(0)

