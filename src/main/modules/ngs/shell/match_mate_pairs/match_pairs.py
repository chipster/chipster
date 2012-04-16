#!/usr/bin/env python

"""
This program will clean out missing mate pairs from a FASTQ file.
Matching mates will be written into a new FASTQ file.

@author: Mikael Karlsson <i8myshoes@gmail.com>
@copyright: 2012
@version: 0.1
@since: 2012-04-13
"""
import sys

print sys.argv

# Open input file
inFile = open(sys.argv[1], "r")

# Open output file
outFile = open(sys.argv[2], "w")

prevRead = None
currRead = None
inCount = 0
outCount = 0

# Match pairs in file
while True:
    # Read in one read 
    currRead = []
    for i in range(4):
        line = inFile.readline()
        if not line:
            break
        currRead.append(line)
    if not line:
        break
    inCount += 1

    # Compare prevRead to currRead
    if prevRead:
        prevHead = prevRead[0].rstrip()
        currHead = currRead[0].rstrip()
        if prevHead[-1] == '1' and currHead[-1] == '2':
            if prevHead[:-2] == currHead[:-2]:
                # Write out to file
                for line in prevRead + currRead:
                    outFile.write(line)
                outCount += 2
            # Unset currRead
            currRead = None

    # Update prevRead
    prevRead = currRead

outFile.close()
inFile.close()

print "In count:", inCount
print "Out count:", outCount
