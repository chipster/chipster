#!/bin/bash

# Script for downloading sheep genome and converting it into correct format for Chipster. 
# Please verify the filesize (SEE THE LAST COMMENT), wget/server seems to cut off files sometimes.

# download assembly report
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCA_000298735.1.assembly.txt

# keep only chromosome identifiers
cat GCA_000298735.1.assembly.txt |grep chromosome |grep -v "#" |cut -d "	" -f 4 > chrs.txt

# download sequence for each chromosome
while read chr; do
	wget "http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?tool=portal&sendto=on&log$=seqview&db=nuccore&dopt=fasta&sort=&val=$chr&extrafeat=0&maxplex=1" -O "$chr.fa"
done < chrs.txt

# append all sequences together
cat CM*.fa > toplevel.fa

# clean headers
sed 's/>gi.*chromosome />/g' toplevel.fa | sed 's/, whole genome shotgun sequence//g' > toplevel-rename.fa

# remove empty rows
grep . toplevel-rename.fa > toplevel-clean.fa

# Result length should equal to "Total sequence length" minus "Total length" of "Unpaced scaffolds", in this case 2587507083
BYTES=$(grep -v ">" toplevel-clean.fa | wc -c); LINES=$(grep -v ">" toplevel-clean.fa | wc -l);expr $BYTES - $LINES
echo $BYTES

mv toplevel-clean.fa Ovis_aries.Oar_v3.1.dna.toplevel.fa
