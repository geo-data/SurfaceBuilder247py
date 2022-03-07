#-------------------------------------------------------------------------------
# Name:        SB247DataClipper
# Purpose:     Extract data for list of IDs from SB247 files, correcting header accordingly
#
# Author:      ajph
#
# Created:     05/02/2019
# Copyright:   (c) ajph 2019
# Licence:     <your licence>
#
# Changelog:
#
#   ajph 28/1/22 Updated to run in Python 3
#
#-------------------------------------------------------------------------------

# Expecting src and dest filenames plus a file containing a list of IDs to retain as arguments

# Sample command:
## python T:\UC1450_pop247nrt\Code\SurfaceBuilder247\SB247DataClipper\SB247DataClipper.py
##    T:\pop247\Pop247Open\Pop247OpenEng2011\SurfaceBuilder247_v2\Data\AccomFilterList.txt
##    T:\pop247\Pop247Open\Pop247OpenEng2011\SurfaceBuilder247_v2\Data\Dests\Dest_Eng_Accom_OTT_2011.csv
##    T:\UC1450_pop247nrt\ESRC_IAA_Secondments\PHE\Training\ForDelegates\Data\SB247\Dests\Dest_Eng_Accom_OTT_2011_Subset.csv
##


import sys, os, string, csv

def main():

    print("\nSB247 Data file ID subsetter, AJPH Feb 2019\n")

    # Parse and check commandline arguments
    invertSelect = False
    fileOffset = 0

    if len(sys.argv) == 4:
        pass
    elif len(sys.argv) == 5 and sys.argv[1] == '-i':
        invertSelect = True
        fileOffset = 1
    else:
        print("""
    Incorrect command line arguments, expecting:
        [-i] <IDListFile> <Infile> <Outfile>
        """)
        sys.exit()

    # Get file names

    idListPath = sys.argv[1 + fileOffset]
    srcPath = sys.argv[2 + fileOffset]
    dstPath = sys.argv[3 + fileOffset]

# Testing
##    idListPath = r"T:\UC1450_pop247nrt\Code\SurfaceBuilder247\SB247DataClipper\Testing\England_WZ.txt"
##    srcPath = r"T:\UC1450_pop247nrt\Code\SurfaceBuilder247\SB247DataClipper\Testing\Dest_EW_Mine+Transp_2011.csv"
##    dstPath = r"T:\UC1450_pop247nrt\Code\SurfaceBuilder247\SB247DataClipper\Testing\Dest_EW_Mine+Transp_2011_Out.csv"
##    invertSelect = False


    # Check file exists
    if not os.path.exists(srcPath):
        print("  Source file not found, exiting")
        sys.exit()

    if not os.path.exists(idListPath):
        print("  ID list file not found, exiting")
        sys.exit()

    if os.path.exists(dstPath):
        print("  Destination file already exists, exiting")
        sys.exit()

    idList = set()
    matchedIDList = []

    # Read the ID list
    with open(idListPath) as idListFile:
        for line in idListFile:
            idList.add(line[:-1])


    # Now read header from source file

    with open(srcPath) as srcFile:

        srcCSV = csv.reader(srcFile)

        header = {}
        headerOrder = []
        inTxt = next(srcCSV)
        linePos = 1

        while inTxt[0] != '':
            header[inTxt[0]] = inTxt[1:]
            headerOrder.append(inTxt[0])
            inTxt = next(srcCSV)
            linePos += 1

        if 'DataBlock' in header:
        # if header.has_key('DataBlock'):
            dataStRow,dataRows = header['DataBlock'][1:3]
        else:
            print("  No DataBlock header line found, exiting")
            sys.exit()

        print('  %s data rows in source file' % (dataRows))


        # Read the lines between the last header text and the first data line.
        fillerRows = []

        while linePos < int(dataStRow) -1:
            fillerRows.append(inTxt)
            inTxt = next(srcCSV)
            linePos += 1

        # Insert last line read, as otherwise this will be dropped
        fillerRows.append(inTxt)

        # Now the file is at the start of the data rows.
        # Read and filter each in turn, retaining them in a temp file.

        with open(dstPath + '.tmp','w') as outFile_tmp:
            tmpCSV = csv.writer(outFile_tmp,lineterminator='\n')

            idCol = int(header['IdentifyingCode'][0])
            filtRowCt = 0


            for (row,inTxt) in enumerate(srcCSV):
                #print ('  Processing row: {}'.format(int(dataStRow) + row))
                
                curID = inTxt[idCol - 1]

                if curID in idList:
                    # Keep this row - matches an ID
                    if not invertSelect:
                        tmpCSV.writerow(inTxt)
                        filtRowCt += 1

                    matchedIDList.append(curID)

                else:
                    if invertSelect:
                        tmpCSV.writerow(inTxt)
                        filtRowCt += 1

    print('\n  %s IDs read' % (len(idList)))
    print('    %s matched a data row' % (len(matchedIDList)))
    print('    %s unmatched' % (len(idList) - len(matchedIDList)))

    print('\n  %s filtered data rows to write out%s' % (filtRowCt,' (inverted - unmatched)'[:invertSelect * 100]))

    if filtRowCt > 0:
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
