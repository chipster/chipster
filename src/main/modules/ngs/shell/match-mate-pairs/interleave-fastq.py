#!/usr/bin/env python

"""
This program will interleave two FASTQ files into one FASTQ file

@author: Mikael Karlsson <i8myshoes@gmail.com>
@copyright: 2012
@version: 0.1
@since: 2012-03-09
"""
import sys

print sys.argv

# Open input files
inFile1 = open(sys.argv[1], "r")
inFile2 = open(sys.argv[2], "r")

# Open output file
outFile = open(sys.argv[3], "w")

# Interleave files
while True:
    for file in inFile1, inFile2:
        for i in range(4):
            line = file.readline()
            if not line:
                break
            outFile.write(line)
    if not line:
        break

outFile.close()
inFile2.close()
inFile1.close()
