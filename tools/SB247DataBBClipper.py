#-------------------------------------------------------------------------------
# Name:        SB247DataClipper
# Purpose:     Extract data for a BB from SB247 files, correcting header accordingly
#
# Author:      ajph
#
# Created:     10/01/2019
# Copyright:   (c) ajph 2019
# Licence:     <your licence>
#-------------------------------------------------------------------------------

# Expecting src and dest filenames as arguments, plus LLX, LLY, URX, URY coordinates
# of clipping BB.

# Sample command:
## python T:\UC1450_pop247nrt\Code\SurfaceBuilder247\SB247DataClipper\SB247DataClipper.py
##    T:\pop247\Pop247Open\Pop247OpenEng2011\SurfaceBuilder247_v2\Data\Dests\Dest_Eng_Accom_OTT_2011.csv
##    T:\UC1450_pop247nrt\ESRC_IAA_Secondments\PHE\Training\ForDelegates\Data\SB247\Dests\Dest_Eng_Accom_OTT_2011.csv
##    365000 152000 389000 176000
##


import sys, os, string, csv

def main():

    print("\nSB247 Data file clipper, AJPH Jan 2019\n")

    # Parse and check commandline arguments

    if len(sys.argv) != 7:
        print("""
    Incorrect number of command line arguments, expecting:
        <Infile> <Outfile> <BB LL X> <BB LL Y> <BB UR X> <BB UR Y>
        """)
        sys.exit()

    # Get file names

    srcPath = sys.argv[1]
    dstPath = sys.argv[2]

    # Check file exists
    if not os.path.exists(srcPath):
        print("  Source file not found, exiting")
        sys.exit()

    if os.path.exists(dstPath):
        print("  Destination file already exists, exiting")
        sys.exit()

    # Read BB
    BB = []

    for b in range(3,7):
        if not isNumber(sys.argv[b]):
            print("  Non-numeric BB coordinates, exiting")
            sys.exit()
        else:
            BB.append(float(sys.argv[b]))

    # Check for coordinate ordering

    if BB[0] >= BB[2] or BB[1] >= BB[3]:
        print("  BB coords are incorrectly ordered, exiting")
        sys.exit()


    # Now read header from source file

    with open(srcPath) as srcFile:

        srcCSV = csv.reader(srcFile)

        header = {}
        headerOrder = []
        inTxt = srcCSV.next()
        linePos = 1

        while inTxt[0] != '':
            header[inTxt[0]] = inTxt[1:]
            headerOrder.append(inTxt[0])
            inTxt = srcCSV.next()
            linePos += 1

        if header.has_key('DataBlock'):
            dataStRow,dataRows = header['DataBlock'][1:3]
        else:
            print("  No DataBlock header line found, exiting")
            sys.exit()

        print('  %s data rows in source file' % (dataRows))


        # Read the lines between the last header text and the first data line.
        fillerRows = []

        while linePos < int(dataStRow) -1:
            fillerRows.append(inTxt)
            inTxt = srcCSV.next()
            linePos += 1

        # Insert last line read, as otherwise this will be dropped
        fillerRows.append(inTxt)

        # Now the file is at the start of the data rows.
        # Read and filter each in turn, retaining them in a temp file.

        with open(dstPath + '.tmp','w') as outFile_tmp:
            tmpCSV = csv.writer(outFile_tmp,lineterminator='\n')

            xCol = int(header['X'][0])
            yCol = int(header['Y'][0])
            filtRowCt = 0

            for inTxt in srcCSV:
                ptX = float(inTxt[xCol - 1])
                ptY = float(inTxt[yCol - 1])

                if ptX >= BB[0] and ptX <= BB[2] and ptY >= BB[1] and ptY <= BB[3]:
                    # Keep this row - falls within the BB
                    tmpCSV.writerow(inTxt)
                    filtRowCt += 1

    if filtRowCt > 0:
        print('  %s filtered data rows in output file' % (filtRowCt))
        header['DataBlock'][2] = str(filtRowCt)

        # Now need to write the updated header and append the filtered data rows
        import shutil

        with open(dstPath,'w') as dstFile:
            dstCSV = csv.writer(dstFile,lineterminator='\n')
            # Header first
            for hdLine in headerOrder:
                header[hdLine].insert(0,hdLine)
                dstCSV.writerow(header[hdLine])

            # Now filler rows, no trimming of EOL char done for these
            for l in fillerRows:
                dstCSV.writerow(l)

            # Finally the data rows
            with open(dstPath + '.tmp') as outFile_tmp:
                shutil.copyfileobj(outFile_tmp,dstFile)

        print('  Complete, see %s' % (dstPath))

    else:
        print('  No rows fell within the BB coordinates')


    # Delete the temp file.
    os.remove(dstPath + '.tmp')


def isNumber(s):
    try:
        float(s)
    except ValueError:
        return False
    else:
        return True


if __name__ == '__main__':
    main()
