#!/bin/bash

# Script for downloading chicken genome and converting it into correct format for Chipster. 

# download assembly report
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000002315.3.assembly.txt

# keep only chromosome identifiers
cat GCF_000002315.3.assembly.txt |grep assembled-molecule |grep -v "MT" |cut -d "	" -f 3 > chrs.txt

# download sequence for each chromosome
while read chr; do
	wget "ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_other/Gallus_gallus/Gallus_gallus-4.0/Primary_Assembly/assembled_chromosomes/FASTA/chr$chr.fa.gz" -O "$chr.fa.gz"
	gunzip "$chr.fa.gz"
done < chrs.txt

# append all sequences together
cat *.fa > toplevel.fa

# clean headers
sed 's/>gi.*chromosome />/g' toplevel.fa | sed 's/, whole genome shotgun sequence//g' > toplevel-rename.fa

mv toplevel-rename.fa Gallus_gallus.Gallus_gallus-4.0.pre.fa

