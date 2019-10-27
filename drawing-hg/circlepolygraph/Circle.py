#! /usr/bin/python

# Python script to take input from command line and process .ocx-mfe file
# for circle polymer graph rendering.


# Import a modules we'll use
import sys
argv = sys.argv
exit = sys.exit
import os
system = os.system
import re
WhiteSpaceSearch = re.compile('\s+')
import string
digits = string.digits


# ######### GET DATA FROM COMMAND LINE INPUT ##################
# Get the name of the input file
InputFile = argv[-1]

# Get the suffix for the input file to see what kind of dot plot we're doing
InputFileList = InputFile.split('.')
InputFileSuffix = InputFileList[-1]
prefix = InputFile[0:(len(InputFile) - len(InputFileSuffix) - 1)]
if InputFileSuffix == 'ocx-mfe':
    Method = 1  # Output from complexes
    permKeyFile = prefix + '.ocx-key'
else:
    Method = 2  # Output from mfe ONLY WORKS FOR A SINGLE STRAND!  (This is a hack)


# Find the output file
OutFileIndex = -1  # Initialize in case the option isn't chosen
if '-svgfile' in argv:
   OutFileIndex  = argv.index('-svgfile') + 1
elif '--svgfile' in argv:
    OutFileIndex = argv.index('--svgfile') + 1
elif '-s' in argv:
    OutFileIndex = argv.index('-s') + 1

if OutFileIndex == -1:
    OutputFile = prefix + '-circle.svg'
else:
    OutputFile = argv[OutFileIndex]


# Find the complex ID
if Method == 1:
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


# Find the permuation ID
if Method == 1:
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


# Find the executable name
ExecutableIndex = -1  # Initialize in case the option isn't chosen
if '-executable' in argv:
   ExecutableIndex  = argv.index('-executable') + 1
elif '--executable' in argv:
    ExecutableIndex = argv.index('--executable') + 1
elif '-e' in argv:
    ExecutableIndex = argv.index('-e') + 1

if ExecutableIndex == -1:
    ExecutableCommand = 'circlepolygraph'
else:
    ExecutableCommand = argv[ExecutableIndex]
# #############################################################


# Input and output files for circlepolygraph
CirclePolyGraphInputFile = 'temp.in'
CirclePolyGraphOutputFile = 'temp-mfe_polygraph.svg'


# Open file with design target or MFE structure
f = open(InputFile,'r')

if Method == 1:
    # Blow through comments and newlines until we find the number of strands
    line = f.readline()
    while line[0:19] != '% Number of strands':
        line = f.readline()


    # Number of strands is the fifth element in the split line
    LineData = WhiteSpaceSearch.split(line)
    K = int(LineData[4])


    # Next line is "% id sequence"
    line = f.readline()


    # Next K lines have sequences
    # Sequences[i] is [sequence ID, sequence, length of sequence]
    Sequences = []
    for i in range(K):
        line = f.readline()
        LineData = WhiteSpaceSearch.split(line)
        Sequences.append([int(LineData[1]),LineData[2],len(LineData[2])])


    # Search for complex and perm id we're after
    TargetString = '% complex' + str(ComplexID) + '-order' + str(PermID)        
    TargetStringLen = len(TargetString)
    line = f.readline()
    while line[0:TargetStringLen] != TargetString and line != '':
        line = f.readline()

    if line == '': # Didn't find record
        print '\nComplex/order not found in input file.\nExiting...\n'
        exit(1)

    line = f.readline()  # This is the total number of bases

else:  # Method == 2
    K = 1; # Only one strand
    
    # Blow through comments and newlines until we find the sequence
    line = f.readline()
    while line[0:11] != '% Sequence:':
        line = f.readline()

    # This is the sequence information
    Sequences = []
    LineData = WhiteSpaceSearch.split(line)
    Sequences.append([1,LineData[2],len(LineData[2])])


    # Blow through rest of comments until we get to the record
    while line[0] == '%' or line[0] == '\n' or line[0] == '\0':
        line = f.readline()


# We are now on the line with the total number of bases.
line = f.readline()  # This is the free energy
line = f.readline()  # This is the dot paren structure
line = f.readline()  # This is the first base pair
# Now we keep reading the bp information as long as the first character is a number
Pairs = []
while line != '' and (line[0] in digits):
    LineData = WhiteSpaceSearch.split(line)
    Pairs.append([int(LineData[0]) , int(LineData[1])])
    line = f.readline()
        
# Close MFE or design target input file
f.close()
    

if Method == 1:
    # Get the permKey file for this complex
    f = open(permKeyFile,'r')


    # Blow through comments and newlines
    line = f.readline()
    while line[0] == '%' or line[0] == '\0' or line[0] == '\n':
        line = f.readline()


    # Search for complex and perm id we're after
    FoundLine = 0
    while FoundLine == 0 and line != '':
        LineList = WhiteSpaceSearch.split(line)
        Comp = int(LineList[0])
        Perm = int(LineList[1])
        if Comp == ComplexID and Perm == PermID:
            FoundLine = 1
        line = f.readline()

        
    # Close permKey file
    f.close()


    # Line list usually has empty character as last entry from new line; delete it
    if LineList[-1] == '':
        del LineList[-1]


    # First two entries are complex and perm IDs, so cut them out
    permKey = LineList[2:len(LineList)]


    # Convert to integers
    for i in range(len(permKey)):
        permKey[i] = int(permKey[i])
else:
    permKey = [1]


# Make input file
f = open(CirclePolyGraphInputFile,'w')
f.write('%d\n' % K)
for x in Sequences:
    f.write('%s\n' % x[1])
for x in permKey:
    f.write('%d\t' % x)
f.write('\n')
for x in Pairs:
    f.write('%d\t%d\n' % (x[0],x[1]))
f.close()


# Run the C-code
cmd = ExecutableCommand + ' temp'
system(cmd)
cmd = 'mv ' + CirclePolyGraphOutputFile + ' ' + OutputFile
system(cmd)
cmd = 'rm -f ' + CirclePolyGraphInputFile
system(cmd);

exit(0)

