#! /usr/bin/python

# Python script to take input from command line and process .ocx-mfe file
# for secondary structure rendering.
#
# Creates an input file for secondary structure drawing.  The input
# file has the following content on successive lines
#   Number of strands in complex (K)
#   The sequences present in the complex (one per line)
#   The tab-delimited ordering for the complex (permKey)
#   The pairs in the secondary structure in pair list format.
#


# Import the modules we'll use
import sys
argv = sys.argv
exit = sys.exit
import os
system = os.system
import re
WhiteSpaceSearch = re.compile('\s+')
PlusSearch = re.compile('\+')
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
    permKeyFile = prefix + '.ocx-key'

    # Find the complex ID
    if '-complex' in argv:
        CompIndex = argv.index('-complex') + 1
    elif '--complex' in argv:
        CompIndex = argv.index('--complex') + 1
    elif '-c' in argv:
        CompIndex = argv.index('-c') + 1
    else:
        print 'Error: No complex ID!'
        exit(1)
    ComplexID = int(argv[CompIndex])


    # Find the permuation ID
    if '-order' in argv:
        PermIndex = argv.index('-order') + 1
    elif '--order' in argv:
        PermIndex = argv.index('--order') + 1
    elif '-o' in argv:
        PermIndex = argv.index('-o') + 1
    else:
        print 'Error: No order ID!'
        exit(1)
    PermID = int(argv[PermIndex])

elif InputFileSuffix != 'mfe':
    print 'Error: Invalid input file.'
    exit(1)


FileNameIndex = -1
if '-filename' in argv:
    FileNameIndex = argv.index('-filename') + 1
elif '--filename' in argv:
    PermIndex = argv.index('--filename') + 1
elif '-o' in argv:
    PermIndex = argv.index('-f') + 1

if FileNameIndex == -1:
    FileName = prefix + '-ssdraw2d-input.in'
else:
    FileName = argv[FileNameIndex]

    
    

if InputFileSuffix == 'ocx-mfe':
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



    # Open file with design target or MFE structure
    f = open(InputFile,'r')

    # Blow through comments and newlines until we find the number of strands
    line = f.readline()
    while line[0:19] != '% Number of strands' and line != '':
        line = f.readline()
    if line == '':
        print 'Error in input file.  No number of strands line!'
        exit(0)


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
        Sequences.append(LineData[2])
        

    # Search for complex and perm id we're after
    TargetString = '% complex' + str(ComplexID) + '-order' + str(PermID)        
    TargetStringLen = len(TargetString)
    line = f.readline()
    while line[0:TargetStringLen] != TargetString and line != '':
        line = f.readline()

    if line == '': # Didn't find record
        print '\nComplex/order not found in input file.\nExiting...\n'
        exit(1)


else: # .mfe input file
    # Open file with design target or MFE structure
    f = open(InputFile,'r')

    # Blow through comments and newlines until we find the number of strands
    line = f.readline()
    while line[0:12] != '% Sequence: ' and line != '':
        line = f.readline()
    if line == '':
        print 'Error in input file.  No number of strands line!'
        exit(0)

    # This line has the sequences
    LineList = WhiteSpaceSearch.split(line)
    Sequences = PlusSearch.split(LineList[2])

    # Number of strands is the fifth element in the split line
    K = len(Sequences)

    # Keep reading until we get a row of % signs
    while line[0:5] != '% %%%' and line != '':
        line = f.readline()
    if line == '':
        print 'Error in input file.  No number of strands line!'
        exit(0)

    # Generate perm key
    permKey = []
    for x in Sequences:
        if x in permKey:
            permKey.append(permKey.index(x))
        elif len(permKey) == 0:
            permKey.append(1)
        else:
            permKey.append(max(permKey)+1)
    

# Now read in the pair list
line = f.readline()  # This is the total number of bases
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


# Make input file
f = open(FileName,'w')
f.write('%d\n' % K)
for x in Sequences:
    f.write('%s\n' % x)
for x in permKey:
    f.write('%d\t' % x)
f.write('\n')
for x in Pairs:
    f.write('%d\t%d\n' % (x[0],x[1]))
f.close()


#*******
# Joe added, first time ever touching python, be wary!
if '-ppairs' in argv or '--ppairs' in argv or '-p' in argv:
	Probs = {}
	PairAssoc = {}
	extrafilename = ''
	if InputFileSuffix == 'ocx-mfe':
		extrafilename = '_'  + str(ComplexID) + '_' + str(PermID)
		f = open(prefix + '.ocx-ppairs','r')
		TargetString = '% complex' + str(ComplexID) + '-order' + str(PermID)        
		TargetStringLen = len(TargetString)
		line = f.readline()
		while line[0:TargetStringLen] != TargetString and line != '':
			line = f.readline()

		if line == '': # Didn't find record
			print '\nComplex/order not found in input file.\nExiting...\n'
			exit(1)

		
	else:
		f = open(prefix + '.ppairs','r')
		line = f.readline()
	
	while line[0] == '%' or line.strip() == '':
		line = f.readline()
	
	length = int(line) #length
	line = f.readline()
	while line.strip() != '' and (line[0] in digits):
		LineData = WhiteSpaceSearch.split(line)
		Probs[str(int(LineData[0])) + ","+ str(int(LineData[1]))] = float(LineData[2])
		Probs[str(int(LineData[1])) + ","+ str(int(LineData[0]))] = float(LineData[2])

		line = f.readline()
	f.close()	
	PairsName = prefix + extrafilename + '.bp'
	f = open(PairsName,'w')
	for x in Pairs:
		PairAssoc[x[0]] = x[1]
		PairAssoc[x[1]] = x[0]
		#f.write('%f\n' % Probs[str(x[0])+","+str(x[1])])
	for i in range(1,length+1):
		if not PairAssoc.has_key(i):
			PairAssoc[i] = length + 1
	#print Probs
	#print Pairs
	#print PairAssoc
	for x in PairAssoc.keys():
		#if x < PairAssoc[x] or PairAssoc[x] == length + 1:
		#print str(x) + "," + str(PairAssoc[x])
		f.write('%f\n' % (Probs[str(x) + "," + str(PairAssoc[x])]))
	f.close()

exit(0)

