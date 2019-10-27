# Reads input data and returns it as the tuple (PairProb, UnPairProb,
# Sequences, PairList, N, K, permKey, Mirror), which has all the
# information needed as data to generate the dot plot.
#
# Justin Bois,
# 8 September 2006
#
# Updated, 18 October, 2007
import sys
import time

# try opening a file, sleeping and repeating a few times if it's not found
def delay_open(fname, mode='r'):
  i=0
  delay=0.1
  f=None
  while f==None and i<5:
    try:
      f=open(fname, mode)
    except IOError,e:
      if e.errno == 2:
        i+=1
        time.sleep(delay)
        sys.stderr.write("Sleeping %1.2f\n"%delay)
        delay*=2
      else:
        break
  if f == None:
    sys.stderr.write("Could not open file %s\n"%fname)
  return f

def getInputDataComplexes(InputFile,cxFile,ocxFile,permKeyFile,MFEFile,ComplexID,PermID,Params,Errors):

    import re
    WhiteSpaceSearch = re.compile('\s+')
    import sys
    exit = sys.exit

    # Define errors
    FAILED_TO_FIND_PERMUTATION_ERROR = Errors[0]
    FAILED_TO_FIND_INITIAL_CONCENTRATIONS = Errors[1]

    # Define paramters
    Method = Params[0]
    NoPermID = Params[1]
    Asymmetric = Params[2]
    LowerTriangle = Params[3]
    CutUnattached = Params[4]
    UseMFE = Params[5]
    getFreeEnergyFromFile = Params[6]

    # Set free energy so we have something to return
    FreeEnergy = 0.0

    # ################### READ INPUT DATA FROM CX FILE #############
    # Open the input file for reading
    f = delay_open(cxFile,'r')

    # Blow through comments and newlines until we find the number of strands
    line = f.readline()
    if len(line)==0:
      print "EOF detected in %s, exiting..."%cxFile
      exit(0)
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


    # Extract stoichometry for Method == 1 or 2
    if Method == 1 or Method == 2:


        # Blow through the rest of the comments to get to complexes and their data
        while (line[0] == '%' or line[0] == '\0' or line[0] == '\n') and line[0] != '':
            line = f.readline()


        # Extract  file
        FoundLine = 0
        while FoundLine == 0 and line != '':
            if line[0] != '%' and line[0] != '\n' and line[0] != '\0':
                LineData = WhiteSpaceSearch.split(line)
                Comp = int(LineData[0])
                if Comp == ComplexID:
                    FoundLine = 1
            line = f.readline()


        # Check to make sure we found the permutation we were after
        if FoundLine == 0:
            print '\nDesired complex not in .cx file.\nExiting....\n'
            exit(FAILED_TO_FIND_PERMUTATION_ERROR)


        # LineData usually has empty character as last entry from new line; delete it
        # --This might be a Python thing and not necessary in Ruby.
        if LineData[-1] == '':
            del LineData[-1]


        # Convert LineData to numbers
        # First chop off complex ID
        LineData = LineData[1:len(LineData)]

        # Convert the stoichiometry entries
        Aj = []  # Column in stoichiometry matrix representing complex
        for i in range(K):
            Aj.append(int(LineData[i]))

        # Read the free energy if necessary
        if getFreeEnergyFromFile and NoPermID: # This means we're plotting complex information
            FreeEnergy = float(LineData[-1])

    f.close()
    # ##############################################################


    # ####################### READ PERM KEY #######################
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

    else : # Strands are indistinguishable, permKey just sequence list of strands
        permKey = []
        for i in range(K):
            permKey.append(i+1)


    # Get N, total number of bases by adding sequence lengths
    N = 0
    for x in permKey:
        N += Sequences[x-1][2]
    # ##############################################################


    # ############## GET FREE ENERGY FROM .ocx FILE ################
    if getFreeEnergyFromFile and (Method == 1 or Method == 2) and not(NoPermID):
        # Open the input file for reading
        f = open(ocxFile,'r')


        # Blow through comments and newlines
        line = f.readline()
        while line[0] == '%' or line[0] == '\0' or line[0] == '\n':
            line = f.readline()


        # Search for complex and perm id we're after
        FoundLine = 0
        while FoundLine == 0 and line != '':
            while line[0] == '%' or line[0] == '\0' or line[0] == '\n':
                line = f.readline()
            LineList = WhiteSpaceSearch.split(line)
            Comp = int(LineList[0])
            Perm = int(LineList[1])
            if Comp == ComplexID and Perm == PermID:
                FoundLine = 1
            line = f.readline()


        # Close .ocx file
        f.close()


        # Line list usually has empty character as last entry from new line; delete it
        if LineList[-1] == '':
            del LineList[-1]


        # The free energy is now the last entry
        FreeEnergy = float(LineList[-1])
    # ##############################################################


    # ################ READ INPUT DATA FROM INPUT FILE #############
    # Open the input file for reading
    f = open(InputFile,'r')

    if Method == 1 or Method == 2:
        # Find record
        if NoPermID:
            TargetString = '% complex' + str(ComplexID)
        else:
            TargetString = '% complex' + str(ComplexID) + '-order' + str(PermID)
        TargetStringLen = len(TargetString)


        line = f.readline()
        while line[0:TargetStringLen] != TargetString and line != '':
            line = f.readline()


        if line == '': # Didn't find record
            print '\nComplex/order not found in input file.\nExiting...\n'
            exit(FAILED_TO_FIND_PERMUTATION_ERROR)


        # Scan through file to get pairing information
        UnPairProb = []
        PairProb = []
        line = f.readline() # This is N (already computed from permKey)
        line = f.readline() # This is the first entry in pair probabilities
        while line[0] != '%' and line != '':
            # At this point, a row of PairProb is [i,j,P], or [base i, base j, pair prob].
            LineData = WhiteSpaceSearch.split(line)
            i = int(LineData[0]) - 1
            j = int(LineData[1]) - 1
            P = float(LineData[2])
            if j == N: # This is an unpair probability
                UnPairProb.append([i,P])
            else:
                PairProb.append([i,j,P])
                if Asymmetric == 1:  # Must have entry twice for asymmetrical case
                    PairProb.append([j,i,P])
            line = f.readline()


    else:  # Method = 3, data goes i,j,P on each line
        # Blow through comments
        line = f.readline()
        while (line[0] == '%' or line[0] == '\0' or line[0] == '\n') and line[0] != '':
            line = f.readline()
        line = f.readline() # This is the first line of input data (prev. was N)


        UnPairProb = []
        PairProb = []
        # Current line is first line of data
        while line != '':
            LineData = WhiteSpaceSearch.split(line)
            i = int(LineData[0]) - 1
            j = int(LineData[1]) - 1
            P = float(LineData[2])
            if j == N: # This is an unpair probability
                UnPairProb.append([i,P])
            else:
                PairProb.append([i,j,P])
            line = f.readline()


    # Close the input file
    f.close()


    if Method == 3: # Go back through input file and get the inital concentrations
        x0 = []  # x0[i] = 1 if concentration is nonzero and 0 otherwise
        f = open(InputFile,'r')
        line = f.readline()
        while line != '' and line[0:17] != '% Initial monomer':
            line = f.readline()
        # Next line begins the list of concentrations
        for i in range(K):
            line = f.readline()
            LineData = WhiteSpaceSearch.split(line)
            if float(LineData[2]) == 0.0:  # The third element is the concentration
                x0.append(0)
            else:
                x0.append(1)

        if line == '':
            print 'Monomer concentrations not present in .fpairs file.  \n\nExiting....\n\n'
            exit(FAILED_TO_FIND_INITIAL_CONCENTRATIONS)
        f.close()
    # #############################################################


    # ################# GENERATE BASE KEY #########################
    # Each x in BaseKey has x[0] = id of base in numbering system from input file
    # and x[1] = number of base when disconnected strands are cut out
    # and x[2] is the number of strand breaks to the left of that base when disconnected
    # strands are cut out (if desired).  x[3] is the strand id the base belongs to
    BaseKey = []
    InputBaseNumber = 0
    BaseNumber = 0
    nStrandBreaks = 0

    if Method == 1 or Method == 2:
        for StrandNumber in permKey:
            if CutUnattached == 1 and Aj[StrandNumber-1] == 0: # strand is cut out
                # Advance input base counter
                for i in range(Sequences[StrandNumber-1][2]):
                    BaseKey.append([InputBaseNumber,-1,-1,StrandNumber])
                    InputBaseNumber += 1
            else: # Strand stays
                for i in range(Sequences[StrandNumber-1][2]):
                    BaseKey.append([InputBaseNumber,BaseNumber,nStrandBreaks,StrandNumber])
                    InputBaseNumber += 1
                    BaseNumber += 1
                nStrandBreaks += 1
    else: # Method == 3
        for StrandNumber in permKey:
            if CutUnattached == 1 and x0[StrandNumber-1] == 0: # strand is cut out
                # Advance input base counter
                for i in range(Sequences[StrandNumber-1][2]):
                    BaseKey.append([InputBaseNumber,-1,-1,StrandNumber])
                    InputBaseNumber += 1
            else: # Strand stays
                for i in range(Sequences[StrandNumber-1][2]):
                    BaseKey.append([InputBaseNumber,BaseNumber,nStrandBreaks,StrandNumber])
                    InputBaseNumber += 1
                    BaseNumber += 1
                nStrandBreaks += 1
    # #############################################################


    # ################ ADJUST PAIR PROBS FOR METHOD ###############
    if Method == 2:
        for x in UnPairProb:
            x[1] /= Aj[BaseKey[x[0]][3]-1]
        if Asymmetric == 1:
            for x in PairProb:
                x[2] /= Aj[BaseKey[x[0]][3]-1]
        else:
            for x in PairProb:
                x[2] /= min([ Aj[BaseKey[x[0]][3]-1] , Aj[BaseKey[x[1]][3]-1] ])
    # #############################################################


    # ############### PRUNE ARRAYS WITH CUT STRANDS ###############
    # Adjust PairProb to reflect the correct base numbering and strand breaks
    # Do this regardless of where or not we prune
    for x in PairProb:
        x.append(BaseKey[x[0]][2])
        x.append(BaseKey[x[1]][2])
        x[0] = BaseKey[x[0]][1]
        x[1] = BaseKey[x[1]][1]



    # Do the same for UnPairProb
    for x in UnPairProb:
        x.append(BaseKey[x[0]][2])
        x[0] = BaseKey[x[0]][1]


    # Pruning is only done if Method = 2 or 3 because otherwise, all strands are present
    if Method == 2:

        # Prune permKey
        for i in range(K):
            if CutUnattached == 1 and Aj[i] == 0: # Strand is cut out
                for j in range(permKey.count(i+1)):
                    permKey.remove(i+1)

        # Get the new number of bases
        N = 0
        for x in permKey:
            N += Sequences[x-1][2]

    if Method == 3:

        # Prune permKey
        for i in range(K):
            if CutUnattached == 1 and x0[i] == 0: # Strand is cut out
                for j in range(permKey.count(i+1)):
                    permKey.remove(i+1)


        # Get the new number of bases
        N = 0
        for x in permKey:
            N += Sequences[x-1][2]
    # #############################################################


    # ################# GET OTHER STRUCTURE DATA ##################
    # Right now, dot-paren enries are ONLY valid for Method = 1, meaning
    # that the strands are distinguishable.  InputFile2 must have a line
    # that has complex ID, perm ID, stoichiometry, free energy, dot-paren struct.

    # Initialize Mirror
    if Asymmetric == 1:
        Mirror = 0
    else:
        Mirror = 1


    if Asymmetric == 0 and LowerTriangle == 1:  # Only bother if we want to fill the lower triangle

        # In two input files are the same, just mirror upper triangle in lower
        if UseMFE and Method == 1:


            # Since we're calculating this, we're not going to mirror
            Mirror = 0


            # Open file with design target or MFE structure
            f = open(MFEFile,'r')

            # Find record
            TargetString = '% complex' + str(ComplexID) + '-order' + str(PermID)
            TargetStringLen = len(TargetString)


            line = f.readline()
            while line[0:TargetStringLen] != TargetString and line != '':
                line = f.readline()


            if line == '': # Didn't find record
                print '\nComplex/order not found in MFE file.\nExiting...\n'
                exit(FAILED_TO_FIND_PERMUTATION_ERROR)


            DotParen = f.readline()  # This is the dot-paren of the structure


            # Close MFE or design target input file
            f.close()


            # Convert dot-paren to pair list
            # PairList is the same as PairProb, but for the other
            # structure and without P (entry 2)
            PairList = []
            LeftList = []
            RightList = []
            BracketLeftList = []
            i = 0  #  Index in the DotParen string (includes strand breaks, '+')
            BaseIndex = 0  # Index of base in permutation (doesn't include strand breaks)
            # While loop steps through list.  Each left bracket is stored.
            # Whenever we get a right bracket, it is necessarily paired with
            # the last left bracket in leftlist.  This pair is documented
            # and the first entry in LeftList is then deleted.
            while BaseIndex < N:
                if DotParen[i] == '(':
                    LeftList.append(BaseIndex)
                    BaseIndex += 1
                elif DotParen[i] == ')':
                    PairList.append([LeftList[-1],BaseIndex])
                    LeftList.remove(LeftList[-1])
                    BaseIndex += 1
                elif DotParen[i] == '.':
                    BaseIndex += 1
                elif DotParen[i] == "{":
                    BracketLeftList.append(BaseIndex)
                    BaseIndex += 1
                elif DotParen[i] == "}":
                    PairList.append([BracketLeftList[-1],BaseIndex])
                    del BracketLeftList[-1]
                    BaseIndex += 1
                i += 1
                # Do not advance the BaseIndex if there is a strand break ('.')


            # Append PairList to include dummy place holder and strand gaps
            for x in PairList:
                x.append(1)  # Dummy place holder
                # These are the appropriate strand gaps, since no strand cutting was done for Method = 1
                x.append(BaseKey[x[0]][2])
                x.append(BaseKey[x[1]][2])

        else:  # Just mirror the upper triangle
            PairList = []
            for x in PairProb:
                NextPair = []
                NextPair.append(x[1]) # index of first base
                NextPair.append(x[0]) # index of second base
                NextPair.append(x[2]) # pair probability
                NextPair.append(x[4]) # strand gaps before first base
                NextPair.append(x[3]) # strand gaps before second base
                PairList.append(NextPair)

    else: # No lower triangle, so PairList is empty
        PairList = []
    # #############################################################


    # Return results
    return (PairProb, UnPairProb, Sequences, PairList, N, K, permKey, Mirror, FreeEnergy)


# ###########################################################################
#                                                                           #
#                       End of getInputDataComplexes.                       #
#                                                                           #
# ###########################################################################





def getInputDataNuPackLib(InputFile, MFEFile, LowerTriangle, UseMFE,
                          FreeEnergyInFile, Errors):

    import re
    WhiteSpaceSearch = re.compile('\s+')
    plusSearch = re.compile('\+')
    import sys
    exit = sys.exit


    # Set free energy so we have something to return
    FreeEnergy = 0.0

    # ################### READ INPUT DATA #########################
    # Open the input file for reading
    f = delay_open(InputFile,'r')

    # Blow through comments and newlines until we find the number sequence
    line = f.readline()
    while line[0:10].lower() != '% sequence'  and line != '':
        line = f.readline()
    if line == '':
        print 'Error: Failed to find sequence in %s!' % InputFile
        exit(1)

    seqline = line

    # Search for % Structure comment
    line = f.readline()
    while line[0:11].lower() != '% structure' and line != '':
        line = f.readline()
    if line == '':
        print 'Error: Failed to find structure in %s!' % InputFile
        exit(1)

    # The sequence is the third element in the split line
    # Sequences[i] is [sequence ID, sequence, length of sequence]
    LineData = WhiteSpaceSearch.split(seqline)
    StrucLineData = WhiteSpaceSearch.split(line)
    seqs = plusSearch.split(LineData[2])
    strucs = plusSearch.split(StrucLineData[2])
    K = len(seqs)  # Number of sequences
    L = len(strucs)

    if K==1 and L > 1: # We need to split up the sequence
      j=0
      newseqs = []
      seq = seqs[0]
      for i in range(L):
        newseqs.append(seq[j:j+len(strucs[i])])
        j+=len(strucs[i])

      seqs = newseqs
      K = len(seqs)

    Sequences = []
    for i in range(K):
        Sequences.append([i+1,seqs[i],len(seqs[i])])
    permKey = range(1,K+1)


    # ################# GENERATE BASE KEY #########################
    # Each x in BaseKey has x[0] = id of base in numbering system from input file
    # and x[1] = number of base when disconnected strands are cut out
    # and x[2] is the number of strand breaks to the left of that base when disconnected
    # strands are cut out (if desired).  x[3] is the strand id the base belongs to
    BaseKey = []
    InputBaseNumber = 0
    BaseNumber = 0
    nStrandBreaks = 0

    for StrandNumber in permKey:
        for i in range(Sequences[StrandNumber-1][2]):
            BaseKey.append([InputBaseNumber,BaseNumber,nStrandBreaks,StrandNumber])
            InputBaseNumber += 1
            BaseNumber += 1
        nStrandBreaks += 1
    # #############################################################




    # Blow through comments till we get to free energy (if necessary)
    if FreeEnergyInFile:
        while line[0:13].lower() != '% free energy' and line != '' and \
              line[:23].lower() != '% ensemble free energy:':
            line = f.readline()
            if line == '':
                print 'Error: Failed to find free energy in %s!' % InputFile
                exit(1)
        # Fourth entry in the line data is the free energy
        LineData = WhiteSpaceSearch.split(line.split(":")[1].strip())
        FreeEnergy = float(LineData[0])


    # Blow through the rest of the comments to get to pairing data data
    while line[0] == '%' or line[0] == '\0' or line[0] == '\n':
        line = f.readline()

    # First line is number of bases
    LineData = WhiteSpaceSearch.split(line)
    N = int(LineData[0])
    line = f.readline()

    # Scan through file to get pairing information
    UnPairProb = []
    PairProb = []
    while line != '':
        # At this point, a row of PairProb is [i,j,P], or [base i, base j, pair prob].
        LineData = WhiteSpaceSearch.split(line)
        i = int(LineData[0]) - 1
        j = int(LineData[1]) - 1
        P = float(LineData[2])
        if j == N: # This is an unpair probability
            UnPairProb.append([i,P,BaseKey[i][2]])
        else:
            PairProb.append([i,j,P,BaseKey[i][2],BaseKey[j][2]])
        line = f.readline()
    # #############################################################


    # ################# GET OTHER STRUCTURE DATA ##################
    # Right now, MFEFile must contain only a list of base pairs, i,j
    # This is typically an MFE structure

    # Initialzed Mirror
    Mirror = 1

    PairList = []

    if LowerTriangle == 1:  # Only bother if we want to fill the lower triangle

        # In two input files are the same, just mirror upper triangle in lower
        if UseMFE:


            # Since we're calculating this, we're not going to mirror
            Mirror = 0


            # Open file with design target or MFE structure
            f = open(MFEFile,'r')


            # Blow through comments and newlines
            line = f.readline()
            while line != '' and (line[0] == '%' or line[0] == '\0' or line[0] == '\n'):
                line = f.readline()


            # Read in pair data from the file, read till end of first record
            while line != '' and line[0] != '%' and line[0] != '\0' and line[0] != '\n':
                LineData = WhiteSpaceSearch.split(line)
                i = int(LineData[0]) - 1
                j = int(LineData[1]) - 1
                PairList.append([j,i,1.0,0,0]) # 1 is pair prob, and 0s for strand gaps
                line = f.readline()

        else:  # Just mirror the upper triangle
            for x in PairProb:
                NextPair = []
                NextPair.append(x[1]) # index of first base
                NextPair.append(x[0]) # index of second base
                NextPair.append(x[2]) # pair probability
                NextPair.append(x[4]) # strand gaps before first base
                NextPair.append(x[3]) # strand gaps before second base
                PairList.append(NextPair)
    # #############################################################


    # Return results
    return (PairProb, UnPairProb, Sequences, PairList, N, K, permKey, Mirror, FreeEnergy)


# ###########################################################################
#                                                                           #
#                       End of getInputDataNuPackLib.                       #
#                                                                           #
# ###########################################################################




def getInputDataTarget(TargetFile, permKey, Sequences, N, K, ComplexID, PermID,
                       Params, Errors, useMfe=False):

    import re
    WhiteSpaceSearch = re.compile('\s+')
    import sys
    exit = sys.exit

    # Define errors
    FAILED_TO_FIND_PERMUTATION_ERROR = Errors[0]
    FAILED_TO_FIND_INITIAL_CONCENTRATIONS = Errors[1]

    # Define paramters
    Method = Params[0]
    NoPermID = Params[1]
    Asymmetric = Params[2]
    LowerTriangle = Params[3]
    CutUnattached = Params[4]
    UseMFE = Params[5]

    # Get Aj from permKey
    Aj = [0]*K
    for x in permKey:
        Aj[x-1] += 1


    # ################# GENERATE BASE KEY #########################
    # Each x in BaseKey has x[0] = id of base in numbering system from input file
    # and x[1] = number of base when disconnected strands are cut out
    # and x[2] is the number of strand breaks to the left of that base when disconnected
    # strands are cut out (if desired).  x[3] is the strand id the base belongs to
    BaseKey = []
    InputBaseNumber = 0
    BaseNumber = 0
    nStrandBreaks = 0

    for StrandNumber in permKey:
        if CutUnattached == 1 and Aj[StrandNumber-1] == 0: # strand is cut out
            # Advance input base counter
            for i in range(Sequences[StrandNumber-1][2]):
                BaseKey.append([InputBaseNumber,-1,-1,StrandNumber])
                InputBaseNumber += 1
        else: # Strand stays
            for i in range(Sequences[StrandNumber-1][2]):
                BaseKey.append([InputBaseNumber,BaseNumber,nStrandBreaks,StrandNumber])
                InputBaseNumber += 1
                BaseNumber += 1
            nStrandBreaks += 1
    # #############################################################


    # ################ READ INPUT DATA FROM TARGET FILE #############
    # Open the input file for reading
    f = delay_open(TargetFile,'r')

    TargetString = '% complex' + str(ComplexID) + '-order' + str(PermID)
    TargetStringLen = len(TargetString)

    line = f.readline()
    if Method == 1:
        while line[0:TargetStringLen] != TargetString and line != '':
          line = f.readline()

        if line == '': # Didn't find record
            print '\nComplex/order not found in input file.\nExiting...\n'
            exit(FAILED_TO_FIND_PERMUTATION_ERROR)
        else:
          print 'Found target string'
        line = f.readline() # This is N (already given in input)
        line = f.readline() # This is the first entry in pair list or dot-paren string
    elif Method == 4:
        while not(line[0].isdigit()) and line[0] != '(' and line[0] != ')' and line[0] != '.' \
                  and line != '':
            line = f.readline()
        if line == '': # Didn't find record
            print '\nNo record found.\nExiting...\n'
            exit(1)
        elif line[0].isdigit(): # This is number of bases, next entry is either dot paren or pair
            line = f.readline()
    else: # Can only work for method = 1 or 4
        print '\n Target not valid for this type of input file.\nExiting...\n'
        exit(1)

    if useMfe:
      line = f.readline()
    # Right now in the TargetFile, we are at the input data, either a dot-paren or pair list
    # Scan through file to get pairing information
    TargetList = []
    if line[0].isdigit(): # pair list style
        TargetUnpaired = []
        for i in range(N):
            TargetUnpaired.append([i])
        while line[0] != '%' and line[0] != '\n' and line != '':
            # At this point, a row of TargetList is [i,j], or [base i, base j].
            LineData = WhiteSpaceSearch.split(line)
            i = int(LineData[0]) - 1
            j = int(LineData[1]) - 1
            TargetList.append([i,j])
            TargetUnpaired.remove([i])
            TargetUnpaired.remove([j])
            line = f.readline()
    else: # Dot-paren style
        DotParen = line
        # Convert dot-paren to pair list
        # PairList is the same as PairProb, but for the other
        # structure and without P (entry 2)
        TargetUnpaired = []
        LeftList = []
        RightList = []
        BracketLeftList = []
        SBracketLeftList = []
        TBracketLeftList = []
        i = 0  #  Index in the DotParen string (includes strand breaks, '+')
        BaseIndex = 0  # Index of base in permutation (doesn't include strand breaks)
        # While loop steps through list.  Each left bracket is stored.
        # Whenever we get a right bracket, it is necessarily paired with
        # the last left bracket in leftlist.  This pair is documented
        # and the first entry in LeftList is then deleted.
        while BaseIndex < N:
            if DotParen[i] == '(':
                LeftList.append(BaseIndex)
                BaseIndex += 1
            elif DotParen[i] == ')':
                TargetList.append([LeftList[-1],BaseIndex])
                LeftList.remove(LeftList[-1])
                BaseIndex += 1
            elif DotParen[i] == '.':
                TargetUnpaired.append([BaseIndex])
                BaseIndex += 1
            elif DotParen[i] == "{":
                BracketLeftList.append(BaseIndex)
                BaseIndex += 1
            elif DotParen[i] == "}":
                TargetList.append([BracketLeftList[-1],BaseIndex])
                del BracketLeftList[-1]
                BaseIndex += 1
            elif DotParen[i] == "[":
                SBracketLeftList.append(BaseIndex)
                BaseIndex += 1
            elif DotParen[i] == "]":
                TargetList.append([SBracketLeftList[-1],BaseIndex])
                del SBracketLeftList[-1]
                BaseIndex += 1
            elif DotParen[i] == "<":
                TBracketLeftList.append(BaseIndex)
                BaseIndex += 1
            elif DotParen[i] == ">":
                TargetList.append([TBracketLeftList[-1],BaseIndex])
                del TBracketLeftList[-1]
                BaseIndex += 1

            i += 1
            # Do not advance the BaseIndex if there is a strand break ('.')
    # #############################################################


    # #############################################################
    # Adjust TargetList to reflect the correct base numbering and strand breaks
    # Do this regardless of where or not we prune
    for x in TargetList:
        x.append(BaseKey[x[0]][2])
        x.append(BaseKey[x[1]][2])
        x[0] = BaseKey[x[0]][1]
        x[1] = BaseKey[x[1]][1]

    # Do the same for TargetUnpaired
    for x in TargetUnpaired:
        x.append(BaseKey[x[0]][2])
        x[0] = BaseKey[x[0]][1]
    # #############################################################


    # ######### GET LOWER TRIANGLE OF TARGET ######################
    if LowerTriangle == 1 and Asymmetric == 0:
        NpairsTarget = len(TargetList)
        for i in range(NpairsTarget):
            NextPair = []
            NextPair.append(TargetList[i][1]) # index of first base
            NextPair.append(TargetList[i][0]) # index of second base
            NextPair.append(TargetList[i][3]) # strand gaps before first base
            NextPair.append(TargetList[i][2]) # strand gaps before second base
            TargetList.append(NextPair)
    # #############################################################



    # Return results
    return (TargetList, TargetUnpaired)

