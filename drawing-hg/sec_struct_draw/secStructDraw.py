#! /usr/bin/python

# Python script to take input from command line and process .ocx-mfe
# or .mfe file for secondary structure rendering.
#
# Creates input files for secondary structure drawing.  The input
# file has the following content on successive lines
#   Number of strands in complex (K)
#   The sequences present in the complex (one per line)
#   The tab-delimited ordering for the complex (permKey)
#   The pairs in the secondary structure in pair list format.
#
#
# If necessary, a .prob files is created, which has the probabilities
# that each base is in the state it's in in the drawing.
#
# These files are then used in a call to ssdraw2d, which generates the
# SVG graphic.
#
#
# The usage is as follows:
#  secStructDraw.py [option flags] prefix.suffix
#
#  Acceptable values for suffix are .in (format as per ssdraw2d
#  standard, described above), .mfe (format per standard NUPACK), or
#  .ocx-mfe (also per standard NUPACK).  By default (unless -svgfile
#  flag is chosen), the SVG graphic is output to prefix.svg.  See
#  flags below for details of usage.
#
#  -complex: REQUIRED if the input file is .ocx-mfe.  The argument of this
#            flag is the complex ID to be drawn.
#  -order: REQUIRED if the input file is .ocx-mfe.  The argument of this
#            flag is the ID of the ordering of the complex to be drawn.
#  -allowoverlap: Generate the canonical 2D structure, even if overlaps occur.
#  -docircle: Do circle polymer graph.  This overrides allowoverlap.
#  -svgfile: Allows specification of the name of the SVG file to which the graphic
#            is written.  The required argument following the flag is the file name.
#  -white: Swap black lines for white in the rendering.
#  -prob: Enables shading of the bases based on how likely they are to be in the 
#         depicted state at equilibrium.  The argument for the flag is the name of
#	  a file with the probabilities.  There are three options.  If
#	  the suffix is .prob, the format is as for the standarde input
#	  into ssdraw2d.  This is simply a list (one per line) of the
#	  probability that a given base is in the state in the
#	  structure being predicted.  The number of lines must equal
#	  the number of bases and appear in the order the bases appear
#	  in the specifies ordering of the complex.  If the suffix is
#	  .ocx-ppairs of .ppairs, the format is per standard NUPACK.
#	  These files are parsed to generate a .prob file.
#  -coords: Output a file containing the coordinates of the bases and do not write
#           and SVG file.  The argument is the coordinate file.
#  -material: Either rna or dna.  Uses actual geometry of these
#             materials.  Default is to use standard non-physical drawing.
#  -stacknicks: Nicked helices are stacked
#  -drawbases: draw letters instead of circles for bases
#  -opacity: Probability is indicated with opacity as opposed to colormap.
#  -freeenergy: (no argument) Parse the input file (if ending in .ocx-mfe or .mfe)
#               for free energy.  Overridden by -energy flag.
#  -energy: Argument of flag is printed on plot as free energy.
#  -energypreamble: Argument of flag is printed before the free energy, followed 
#                   immediately by a single whitespace and then the free energy.
#                   Example: -energypreamble "Free energy:"
#  -units: A string for units of free energy.  The default is kcal/mol.
#  -key: Puts a key for color-coding of bases in SVG graphic.
#  -thymine: Used in conjunction with -key to have T instead of default U.
#  -colorstrands: (currently inactive), colors strand backbones
#  -keybuffer: Makes whitespace on the drawing of the same size the key would be on
#              the right side of the drawing.
#  -colorbarbuffer: Makes whitespace on the drawing of the same size the colorbar would
#                   be on the right side of the drawing.
#  -freeenergybuffer: Makes whitespace on the drawing of the same size the free 
#                     energy text would be on the bottom of the drawing
#  -squarebuffer: Makes white space around the entire figure such that the overall
#                 dimension is still square.  The size of the boundary is given by
#                 whichever would be biggest: the space taken by the colorbar, base 
#                 key, or free energy text.  This flag supercedes all other buffering
#                 flags.
#  -ssdrawcommand: Argument is string which is used to call the secondary structure
#                  drawing program.  The default is ssdraw2d, but, for example, if
#                  pwd is not in your path, you might want to use the flag as
#                  -ssdrawcommand "./.ssdraw2d".
#
# 19 October 2007, Justin Bois

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

elif InputFileSuffix != 'mfe' and InputFileSuffix != 'in': # Unrecognized format
    print 'Error: Invalid input file.'
    exit(1)
    


# Read in options and generate command line call to ssdraw2d
drawCommand = 'ssdraw2d '

if '-ssdrawcommand' in argv:
    drawCommand = argv[argv.index('-ssdrawcommand') + 1] + ' '
else:
    drawCommand = 'ssdraw2d '

if '-docircle' in argv:
    drawCommand += '-docircle '

if '-allowoverlap' in argv:
    drawCommand += '-allowoverlap '

if '-opacity' in argv:
    drawCommand += '-opacity '

if '-white' in argv:
    drawCommand += '-white '

if '-key' in argv:
    drawCommand += '-key '

if '-keybuffer' in argv:
    drawCommand += '-keybuffer '

if '-colorbarbuffer' in argv:
    drawCommand += '-colorbarbuffer '

if '-freeenergybuffer' in argv:
    drawCommand += '-freeenergybuffer '

if '-squarebuffer' in argv:
    drawCommand += '-squarebuffer '

if '-thymine' in argv:
    drawCommand += '-thymine '

if '-stacknicks' in argv:
    drawCommand += '-stacknicks '

if '-drawbases' in argv:
    drawCommand += '-drawbases '

if '-colorstrands' in argv:
    drawCommand += '-colorstrands '

if '-energy' in argv:
    argIndex = argv.index('-energy') + 1
    drawCommand += '-energy ' + argv[argIndex] + ' '
    useEnergyString = 1
else:
    useEnergyString = 0

if '-energypreamble' in argv:
    argIndex = argv.index('-energypreamble') + 1
    drawCommand += '-energypreamble "' + argv[argIndex] + '" '

if '-units' in argv:
    argIndex = argv.index('-units') + 1
    drawCommand += '-units "' + argv[argIndex] + '" '

if '-material' in argv:
    argIndex = argv.index('-material') + 1
    drawCommand += '-material ' + argv[argIndex] + ' '

if '-coords' in argv:
    argIndex = argv.index('-coords') + 1
    drawCommand += '-coords ' + argv[argIndex] + ' '

if '-svgfile' in argv:
    argIndex = argv.index('-svgfile') + 1
    drawCommand += '-svgfile ' + argv[argIndex] + ' '

if '-prob' in argv:
    argIndex = argv.index('-prob') + 1
    probFile = argv[argIndex]
    probFileList = probFile.split('.')
    probFileSuffix = probFileList[-1]
    if probFileSuffix == 'prob':
        drawCommand += '-prob ' + argv[argIndex] + ' '
    elif probFileSuffix != 'ocx-ppairs' and probFileSuffix != 'ppairs':
        print 'Error: Invalid file for probabilities.'
        exit(1)

# Find out if we have to parse free energy from file
parseFreeEnergy = 0
if '-freeenergy' in argv: # Must parse free energy from file
    if useEnergyString or InputFileSuffix == 'in':
        print 'Ingnoring -freeenergy flag'
    else:
        parseFreeEnergy = 1


# Parse the file with the structure data
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
    while line != '' and line[0:TargetStringLen] != TargetString:
        line = f.readline()

    if line == '': # Didn't find record
        print '\nComplex/order not found in input file.\nExiting...\n'
        exit(1)

    # Now read in the pair list
    line = f.readline()  # This is the total number of bases
    line = f.readline()  # This is the free energy

    # Grab the free energy
    if parseFreeEnergy:
        freeEnergy = float(line)
    
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

elif InputFileSuffix == 'mfe':
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

    # Grab the free energy
    if parseFreeEnergy:
        freeEnergy = float(line)

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
if InputFileSuffix != 'in':
    inputFileName = prefix + '-drawing-%s.in'%os.getpid()
    f = open(inputFileName,'w')
    f.write('%d\n' % K)
    for x in Sequences:
        f.write('%s\n' % x)
    for x in permKey:
        f.write('%d\t' % x)
    f.write('\n')
    for x in Pairs:
        f.write('%d\t%d\n' % (x[0],x[1]))
    f.close()


# Parse probability file
if '-prob' in argv and (probFileSuffix == 'ppairs' or probFileSuffix == 'ocx-ppairs'):
    if probFileSuffix == 'ppairs':
        f = open(probFile,'r')
        # Blow through comments and newlines until we find the number of strands
        line = f.readline()
        while line != '' and not(line[0].isdigit()):
            line = f.readline()

    elif probFileSuffix == 'ocx-ppairs':
        f = open(probFile,'r')
        # Search for complex and perm id we're after
        TargetString = '% complex' + str(ComplexID) + '-order' + str(PermID)        
        TargetStringLen = len(TargetString)
        line = f.readline()
        while line != '' and line[0:TargetStringLen] != TargetString:
            line = f.readline()

        if line == '': # Didn't find record
            print '\nComplex/order not found in ocx-ppairs file.\nExiting...\n'
            exit(1)

        line = f.readline()

    # We're now at the number of bases
    Nbases = int(line)

    # Construct plist from pairs
    plist = [Nbases+1]*Nbases
    for x in Pairs:
        plist[x[0]-1] = x[1]
        plist[x[1]-1] = x[0]
    
    line = f.readline()
    # Now we're at the pair probability data, read through and get probs
    probs = [0]*Nbases
    while line != '' and line[0].isdigit():
        LineData = WhiteSpaceSearch.split(line)
        i = int(LineData[0])
        j = int(LineData[1])
        p = float(LineData[2])
        if plist[i-1] == j:
            probs[i-1] = p
            if j <= Nbases:
                probs[j-1] = p
        line = f.readline()

    # Close the ppairs file
    f.close()


    # Write data to .prob file
    probInputFile = prefix + '-drawing-%s.prob'%os.getpid()
    f = open(probInputFile,'w')
    for x in probs:
        f.write('%4.3e\n' % x)
    f.close()

    # Update drawing command
    drawCommand += '-prob ' + probInputFile + ' '
    

if parseFreeEnergy:
    drawCommand += '-energy %.2f ' % freeEnergy
drawCommand += prefix
if InputFileSuffix != 'in':
    drawCommand += '-drawing-%s'%os.getpid()


# Make the drawing
print "drawCommand=%s"%drawCommand
system(drawCommand)

# Clean up the input files
cmd = 'rm -f ' + prefix + '-drawing-%s.in'%os.getpid()
system(cmd)
if '-prob' in argv and (probFileSuffix == 'ppairs' or probFileSuffix == 'ocx-ppairs'):
    cmd = 'rm -f ' + prefix + '-drawing-%s.prob'%os.getpid()
    system(cmd)

if '-svgfile' not in argv:
    if InputFileSuffix != 'in':
        cmd = 'mv ' + prefix + '-drawing-mfe_drawing.svg ' + prefix + '.svg'
    else:
        cmd = 'mv ' + prefix + '-mfe_drawing.svg ' + prefix + '.svg'
        
    system(cmd)

exit(0)

