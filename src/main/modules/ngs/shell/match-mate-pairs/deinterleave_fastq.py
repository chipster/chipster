#!/usr/bin/env python

"""
This program will deinterleave one FASTQ file into two FASTQ files

@author: Mikael Karlsson <i8myshoes@gmail.com>
@copyright: 2012
@version: 0.1
@since: 2012-03-16
"""
import sys

print sys.argv

# Open input file
inFile = open(sys.argv[1], "r")

# Open output files
outFile1 = open(sys.argv[2], "w")
outFile2 = open(sys.argv[3], "w")

# Deinterleave files
while True:
    for file in outFile1, outFile2:
        for i in range(4):
            line = inFile.readline()
            if not line:
                break
            file.write(line)
    if not line:
        break

outFile2.close()
outFile1.close()
inFile.close()
