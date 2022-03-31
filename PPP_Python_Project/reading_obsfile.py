# Imports of Python Libraries


# Constants defined
SAT_STR_LENGTH = 3
LEN_OBS_IN_RINEX = 16
MAX_RNX2_CHAR_IN_LINE = 80
MAX_NUM_OBS_PER_LINE = 5
MAX_NUM_SAT_PER_LINE = 12

# Maximum number of characters in RINEX v2
MAX_RNX2_CHAR_IN_LINE = 80
LAST_POS_NSAT_RINEXV2 = 32
# First position of RX clock offset in epoch line
FIRST_POS_RX_CLK_OFFSET = 68
MAX_NUM_SAT_PER_LINE = 12


# TODO: ADD HEADER
def getEpochV2LineParam(Epochline):
    h = Epochline[10:12]
    m = Epochline[13:15]
    s = Epochline[16:18]

    # string to int conversion
    Hour = int(float(h))
    Minute = int(float(m))
    Second = int(float(s))

    # Get epoch of day (second of day)
    Epoch = Hour * 3600 + Minute * 60 + Second
    Nsat = int(Epochline[30:32])

    # Estimate Number of satellite line (every 12 satellites, go to next line)
    NSatLines = (Nsat / (MAX_NUM_SAT_PER_LINE + 1)) + 1

    return Epoch, Nsat, NSatLines


##
# Name		: getNobsHeaderRnxV2
# Purpose	: Extract number of observables from header:
#
def getNobsHeaderRnxV2(RnxLine):
    # Number of observables found in header and flag indicating if it is line containing this number
    NObsFlag = False
    Nobs = 0
    try:
        Nobs = int(RnxLine[4:6])
        NObsFlag = True
    except:
        NObsFlag = False
    # End of try

    # Return, flag and number of observables
    return NObsFlag, Nobs


##
# Name		: getStrOfSat_V2
# Purpose	: get the string containing satellite name as in epoch line of RINEX v2
#
def getStrOfSat_V2(EpochLine, NSatLines, Nsat):
    # Definition of local variables
    FirstPos = 0
    LastPos = 0
    SAT_STR_LENGTH = 3
    SatLine = ''

    # If epoch line is only one
    if NSatLines == 1:

        # Eg:  20 09 04 00 00 00.0000000  0 10G10G13G17G32G15G24G12G25G20G19        .000000000
        # Firs position of the satellites
        FirstPos = LAST_POS_NSAT_RINEXV2

        # Last position of the satellites
        LastPos = LAST_POS_NSAT_RINEXV2 + Nsat * SAT_STR_LENGTH

        # get the string with all Satellites
        SatLine = EpochLine[FirstPos:LastPos]

    else:

        # loop in all lines containing Satellite ID
        for jLine in range(int(NSatLines)):

            # If index is not the last value of the loop
            if jLine != (int(NSatLines) - 1):

                # Firs position of the satellites
                FirstPos = LAST_POS_NSAT_RINEXV2 + jLine * MAX_RNX2_CHAR_IN_LINE

                # Last position of the satellites
                LastPos = FIRST_POS_RX_CLK_OFFSET + jLine * MAX_RNX2_CHAR_IN_LINE

            # If last value of the loop
            else:

                # Firs position of the satellites
                FirstPos = LAST_POS_NSAT_RINEXV2 + jLine * MAX_RNX2_CHAR_IN_LINE

                # Last position of the satellites
                LastPos = FirstPos + (Nsat - MAX_NUM_SAT_PER_LINE * jLine) * SAT_STR_LENGTH

            # End of  if (jLine != (Parameter.NSatLines-INT_ONE)):

            # IF first line
            if jLine == 0:

                # Copy the satellites string
                SatLine = EpochLine[FirstPos:LastPos]

            else:

                # Concatenate the satellites strings
                SatLine += EpochLine[FirstPos:LastPos]
            # end of if

        # End of for jLine in

    # End of if Parameter.NSatLines == 1

    return SatLine


# End of def  getStrOfSat_V2


# TODO: ADD HEADER
def getObservablesFromheader(Fname):
    # Define dictionary
    Header = {}

    Header['System'] = {}
    Header['NSystem'] = {}

    # Input RINEX OBS File handle
    FidRnxObsIn = open(Fname, 'r')

    # Parse each line of the input RINEX file
    RnxObsreadlines = FidRnxObsIn.readlines()

    # Initialization of the RINEX file lines counter
    Nline = 0

    # Number of GNSS constellatinons in file
    NSysCounter = 0

    # Flag indicating end of header has been reached or not
    OutOfHeaderFlag = False

    # Loop over all lines in the input RINEX file
    for RnxLine in RnxObsreadlines:

        # If line is empty, break the loop
        if not RnxLine: break
        # end if not line

        # If out of Header flag is false for input
        if not OutOfHeaderFlag:

            # If line contains RINEX vs 2 header SYS / # / OBS TYPES
            if RnxLine.find('# / TYPES OF OBSERV') != -1:

                # get RINEX v2 observation type in char form
                # Check if number of observables has been defined
                NObsFlag, NObs = getNobsHeaderRnxV2(RnxLine)

                # if the Number of observation is defined (2nd line)
                if NObsFlag:
                    # Fill the Header structure Number of observables, no GNSS ID in version 2
                    # Header['System'].update({NSysCounter: {'SysCode': SysCode, 'GnssAcronym': GnssAcronym, 'NObs': NObs,'Observables': {}}})
                    # Header['System'][NSysCounter]['NObs'] = NObs
                    Header['System'].update({NSysCounter: {'NObs': NObs, 'Observables': {}}})

                    # Loop over all observables
                    # a range from 0 to Nobs at interval of 1
                    for jObs in range(0, int(Header['System'][NSysCounter]['NObs'])):
                        # Get the observable ID in RINEX V2  e.g. L1
                        ObsId = RnxLine[10 + 6 * jObs: 12 + 6 * jObs]

                        # Fulfill the Header structure System Observables
                        # Header['System'][NSysCounter]['Observables'][jObs] = ObsId
                        Header['System'][NSysCounter]['Observables'].update({jObs: ObsId})

                    # End of for jObs in range

                    # Add one to the system counter
                    NSysCounter += 1

                    # Save the total number of GNSS systems in file
                    Header['NSystem'] = NSysCounter

            elif RnxLine.find('END OF HEADER') != -1:
                # Change Out of RINEX Header indicator flag
                OutOfHeaderFlag = True
                break

        # End of if strstr
        else:
            break
            # End of if OutOfHeaderFlag

        # End of for RnxLine in RnxObsreadlines:

    # Close file
    FidRnxObsIn.close()

    return Header


# TODO: ADD HEADER
def readRinexFileV2(FullFileName):
    print("Reading: Rinex OBS Version 2 --> ", FullFileName)
    print("\n")

    OBSstruct = {}
    Header = getObservablesFromheader(FullFileName)

    # RINEX OBS File handle
    FidRnxObs = open(FullFileName, 'r')

    # Flags for RINEX IN and RINEX OUT
    OutOfHeaderFlag = False

    # Read line in RINEX input
    RnxLine = FidRnxObs.readline().rstrip()

    # Loop over all lines in the input RINEX file
    # for line in RnxObsreadlinesIN:
    while True:

        # If line is empty, break the loop
        if not RnxLine: break
        # end if not line

        # If out of Header flag is false for input
        if not OutOfHeaderFlag:

            # If line contains end of header
            if RnxLine.find('END OF HEADER') != -1:
                # Change Out of RINEX Header indicator flag
                OutOfHeaderFlag = True

            # End of if

            # Read line in RINEX input
            RnxLine = FidRnxObs.readline().rstrip()

            continue

        else:

            # EPOCH LINE + PRN LISTS

            # Get Epoch line information: Second of the day, Epoch flag and Num of sat
            Epoch, Nsat, NSatLines = getEpochV2LineParam(RnxLine)

            val1 = int(MAX_RNX2_CHAR_IN_LINE * NSatLines)
            mySpace = " "

            # Create the epoch lines block
            EpochLines = mySpace * val1
            # print(EpochLines)
            # exit()

            # Loop in all lines containing satellite list (epoch and continuation)
            # NSatLines - Line per satellite
            for jLineSat in range(int(NSatLines)):

                # initial positions of the line in string (0,81,162,...)
                InitialSatLinePos = 0

                # final positions of the line in string
                FinalSatLinePos = InitialSatLinePos + len(RnxLine)

                # If first epoch line
                if jLineSat == 0:

                    # Copy the first line of epoch (without '\n' before due to rstrip no need of -1)
                    EpochLines = RnxLine[InitialSatLinePos:FinalSatLinePos] + \
                                 ' ' * (MAX_RNX2_CHAR_IN_LINE - FinalSatLinePos)

                else:

                    # Concatenate following lines (without '\n'before due to rstrip no need of -1)
                    EpochLines += RnxLine[InitialSatLinePos: FinalSatLinePos]

                # End of if jLineSat == 0

                # Read the following line
                RnxLine = FidRnxObs.readline().rstrip()

            # End of for jLineSat in range(NSatLines):

            # END OF EPOCH LINES + PRN LISTS

            # String of satellites for the current epoch (RINEX version 2)
            SatellitesStr = getStrOfSat_V2(EpochLines, int(NSatLines), Nsat)

            ## LOOP IN ALL SATELLITES IN BLOCK
            # For each satellite in the block
            for jSat in range(Nsat):

                # Get satellite ID
                SatId = SatellitesStr[SAT_STR_LENGTH * jSat: SAT_STR_LENGTH * (jSat + 1)]

                # Reset the line in each satellite (there may be many lines with observations per satellite)
                jLinePerSat = 0

                # Loop in all observables : Nobs
                for jObs in range(Header['System'][0]['NObs']):

                    # Get Obs Id
                    ObsId = Header['System'][0]['Observables'][jObs]

                    # First position of the observable string
                    PosFirstCharOfObs = LEN_OBS_IN_RINEX * (jObs - jLinePerSat * MAX_NUM_OBS_PER_LINE)

                    # Last position of the observable string
                    PosLastCharOfObs = LEN_OBS_IN_RINEX * (jObs - jLinePerSat * MAX_NUM_OBS_PER_LINE + 1) - 2

                    # If position exceeds line length
                    if len(RnxLine) < PosLastCharOfObs:

                        # Define observation as empty
                        obsIn = ' '

                    # else, if last observable of line is missing
                    else:

                        # Get observations string (if lol whitespaces characters)
                        obsIn = RnxLine[PosFirstCharOfObs: PosLastCharOfObs]

                    # End of if len(RnxLine) < PosLastCharOfObs:

                    # If dictionary already contains epoch info
                    if OBSstruct.get(Epoch):

                        # If dictionary already contains PRN info
                        if OBSstruct[Epoch].get(SatId):

                            OBSstruct[Epoch][SatId].update({ObsId: obsIn})

                        # else, if dictionary did not contain PRN info
                        else:

                            OBSstruct[Epoch].update({SatId: {ObsId: obsIn}})

                        # End of if OBSstruct[Epoch].get(SatId):

                    # If no info for epoch
                    else:

                        OBSstruct.update({Epoch: {SatId: {ObsId: obsIn}}})

                    # End of if OBSstruct.get(Epoch)

                    # If last observable of line, or last observable of list
                    if (((jObs + 1) % MAX_NUM_OBS_PER_LINE) == 0) or \
                            ((jObs + 1) == Header['System'][0]['NObs']):
                        # Read the following line
                        RnxLine = FidRnxObs.readline().rstrip()

                        # Add one to line counter per satellite
                        jLinePerSat += 1

                    # End of if (jObs + 1) ....

                # End of for jObs in range(Header['System'][0]['NObs'])

            # End of for jSat in range(Nsat):

        # End of if OutOfHeaderFlag == FLAG_FALSE:

    # End of while True:

    FidRnxObs.close()

    return OBSstruct

    # End of function readRinexFileV2()


# Define path and file names
# To read Rinex OBS file V2
ScenarioPath = "D:/PPP/"
Path2Data = ScenarioPath + "DATA/"
FileName = "TLSA2480.20o"



"""
GeneratedOBS.close()

# Well Sort the Epochs : In Ascending order
# Note: The Header is written at the last line instead of the first line
# Make sure to Read only upto (LastLine minus 1)
with open(ScenarioPath + 'OUTPUT/TempOBS.20o') as fin, open(ScenarioPath + 'OUTPUT/generatedOBS.20o', 'w') as fout:
    fout.writelines(sorted(fin))
"""

def run():
    print("\n\n<><><><><><><><><><><><><><><><><><><><><><><>")
    print("Starting Processing: Rinex OBS Version 2")
    print("<><><><><><><><><><><><><><><><><><><><><><><><>\n\n")

    """
    # Write Header of the output Rinex file (NAV)
    GeneratedOBS  =   open(ScenarioPath + 'OUTPUT/TempOBS.20o',"w")
    #GeneratedOBS.write(" Epoch \tSatId  \tL1  \tL2  \tC1  \tP2  \tS1  \tS2 \n")
    """

    print("\nWriting Rinex OBS Version 2 into a file ...")

    # Invoke a function that will read all the Rinex OBS file
    # This function returns a dictionary with:
    # Each Epoch, Each PRN in that Epoch, and Each OBS in that PRN
    OBSstruct = readRinexFileV2(Path2Data + FileName)

    # Loop through the Dictionary and print values for each epoch and GPS PRN
    response = {}
    for key, value in OBSstruct.items():
        CurrentEpoch = key
        for i, j in value.items():
            CurrentPRN = i
            MyListOfOBS = []
            for m, n in j.items():
                if key == CurrentEpoch and i == CurrentPRN:

                    if len(n.strip()) != 0:  # Check if a field is not empty
                        # print  key , ':' , i, ':', m , ':' , n, len(n.strip())
                        MyListOfOBS.append(float(n))
                    else:
                        # If value is not present: assign -99999.0
                        MyListOfOBS.append(-99999.0)

            #print(CurrentEpoch, " : ", CurrentPRN, " : ", MyListOfOBS)
            response[CurrentEpoch] = {"sat": CurrentPRN, "obs": MyListOfOBS}

    return response

    """
    #Epoch SatId P2 L2 L1 S2 C1 S1 
    GeneratedOBS.write('{0:8} {1:8} {2: .3f} {3: .3f}   {4: .3f}   {5: .3f}   {6: .3f}   {7: .3f} \n'.format(CurrentEpoch,CurrentPRN,\
                MyList[2],MyList[1],MyList[4],MyList[0],MyList[5],MyList[3]) )
    """