#! /usr/bin/python

# Python script to make a single .fpairs file from multiple .ppairs files,
# each for one strand
#
# Usage: makefpairsFromppairs.py prefix nstrands fpairsFilePrefix
#
# Makes an .fpairs file and a dummy .cx file for use with DotPlot.py
# The dummy .cx file just has comment lines with number of strands and sequences.
#
#  prefix: The input files to be processed are called prefix_0.ppairs,
#    prefix_1.ppairs, and so on.  There are a total of nstrands of them.
#  nstrands: The number of files from which data are extracted
#  fpairsFilePrefix: The name of the output fpairs file.  In almost all
#    circumstances, this should be the same as prefix.
#
# NOTE: No temperature data are written to the fpairs file.

# Import a modules we'll use
import sys
argv = sys.argv
exit = sys.exit
import re
WhiteSpaceSearch = re.compile('\s+')


# Get parameters from command line
fpairsFilePrefix = argv[-1]
nstrands = int(argv[-2])
prefix = argv[-3]

fpairsFile = fpairsFilePrefix + '.fpairs'
cxFile = fpairsFilePrefix + '.cx'

sequences = []
currentLen = 0
PairProbs = []
for strandID in range(nstrands):
    ppairsFile = '%s_%d.ppairs' % (prefix,strandID)
    f = open(ppairsFile,'r')

    # Blow through comment lines until we get to the sequence
    line = f.readline()
    while len(line) > 10 and line[0:10] != '% Sequence' and line != '':
        line = f.readline()

    if line == '':
        print 'Error!  No sequence found in input file %s!' % ppairsFile
        exit(1)

    # Pull out the sequence
    LineList = WhiteSpaceSearch.split(line)
    sequences.append(LineList[2])

    # Read on until the sequence length
    while line != '' and not(line[0].isdigit()):
        line = f.readline()

    if line == '':
        print 'Error!  No data found in input file %s!' % ppairsFile
        exit(1)

    # This is the number of bases in the sequence
    LineList = WhiteSpaceSearch.split(line)
    seqlen = int(LineList[0])

    # Now go through the data
    line = f.readline()
    while line != '' and line[0].isdigit():
        LineList = WhiteSpaceSearch.split(line)
        i = int(LineList[0]) + currentLen
        j = int(LineList[1])
        p = float(LineList[2])
        if j == seqlen + 1:
            j = -1
        else:
            j += currentLen
            PairProbs.append([j,i,p])
        PairProbs.append([i,j,p])
        line = f.readline()

    f.close()
    currentLen += seqlen
    
# Go through and fix unpaired probs
for i in range(len(PairProbs)):
    if PairProbs[i][1] == -1:
        PairProbs[i][1] = currentLen + 1

# Write dummy .cx file
f = open(cxFile,'w')
f.write('% This is an dummy .cx file generated from multiple .ppairs files.\n')
f.write('%% Number of strands: %d\n' % nstrands)
f.write('% id sequence\n')
# Write sequences
for i in range(nstrands):
    f.write('%%  %d %s\n' % (i+1,sequences[i]))
f.write('%\n')
f.close()


# Write fpairs file
f = open(fpairsFile,'w')
f.write('% This is an .fpairs file generated from multiple .ppairs files.\n')
# Put in dummy concentrations
f.write('% Initial monomer concentrations:\n')
for i in range(nstrands):
    f.write('%%   %d: %8.7e Molar\n' % (i+1,1e-6))
f.write('%\n')
f.write('%% Number of strands: %d\n' % nstrands)
f.write('% id sequence\n')

# Write sequences
for i in range(nstrands):
    f.write('%%  %d %s\n' % (i+1,sequences[i]))
f.write('%\n')

# Write number of bases
f.write('%d\n' % currentLen)

# Write pair probabilities
for x in PairProbs:
    f.write('%d\t%d\t%8.7e\n' % (x[0],x[1],x[2]))
f.close()

exit(0)
    
