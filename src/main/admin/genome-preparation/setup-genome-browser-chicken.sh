#!/bin/bash

# Script for downloading sheep genome and converting it into correct format for Chipster. 
# Please verify the filesize (SEE THE LAST COMMENT), wget/server seems to cut off files sometimes.

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

# remove empty rows
grep . toplevel-rename.fa > toplevel-clean.fa

# Result length should equal to "Total sequence length" minus "Total length" of "Unpaced scaffolds", in this case 2587507083
BYTES=$(grep -v ">" toplevel-clean.fa | wc -c); LINES=$(grep -v ">" toplevel-clean.fa | wc -l);expr $BYTES - $LINES
echo $BYTES

mv toplevel-clean.fa Ovis_aries.Oar_v3.1.dna.toplevel.fa
