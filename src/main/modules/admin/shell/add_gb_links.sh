#!/bin/bash 

TMPDIR=$1

#cd /opt/chipster/tools/genomes
cd $TMPDIR

FILE=./genomebrowser/Arabidopsis_thaliana/TAIR10.22/genome.yaml
if [ -e $FILE ] 
then
	rm $FILE
	echo $FILE >> genomes.tmp
	echo "species: Arabidopsis_thaliana" >> $FILE
	echo "version: (TAIR10.22)" >> $FILE
	echo "ensemblBrowserUrl: http://plants.ensembl.org/Arabidopsis_thaliana/Location/View?r=[CHR]%3A[START]-[END]" >> $FILE
	echo "ucscBrowserUrl: " >> $FILE
	echo "sortId: plant" >> $FILE
fi

FILE=./genomebrowser/Bos_taurus/UMD3.1.75/genome.yaml
if [ -e $FILE ] 
then
	rm $FILE
	echo $FILE >> genomes.tmp
	echo "species: Cow" >> $FILE
	echo "version: (UMD3.1.75)" >> $FILE
	echo "ensemblBrowserUrl: http://www.ensembl.org/Bos_taurus/Location/View?r=[CHR]%3A[START]-[END]" >> $FILE
	echo "ucscBrowserUrl: http://genome.ucsc.edu/cgi-bin/hgTracks?clade=mammal&org=Cow&db=bosTau6&position=chr[CHR]%3A[START]-[END]" >> $FILE
	echo "sortId: other" >> $FILE
fi

FILE=./genomebrowser/Canis_familiaris/CanFam3.1.75/genome.yaml
if [ -e $FILE ] 
then
	rm $FILE
	echo $FILE >> genomes.tmp
	echo "species: Dog" >> $FILE
	echo "version: (CanFam3.1.75)" >> $FILE
	echo "ensemblBrowserUrl: http://www.ensembl.org/Canis_familiaris/Location/View?r=[CHR]%3A[START]-[END]" >> $FILE
	echo "ucscBrowserUrl: http://genome.ucsc.edu/cgi-bin/hgTracks?clade=mammal&org=Dog&db=canFam3&position=chr[CHR]%3A[START]-[END]" >> $FILE
	echo "sortId: other" >> $FILE
fi

FILE=./genomebrowser/Drosophila_melanogaster/BDGP5.75/genome.yaml
if [ -e $FILE ] 
then
	rm $FILE
	echo $FILE >> genomes.tmp
	echo "species: Drosophila melanogaster" >> $FILE
	echo "version: (BDGP5.75)" >> $FILE
	echo "ensemblBrowserUrl: http://www.ensembl.org/Drosophila_melanogaster/Location/View?r=[CHR]%3A[START]-[END]" >> $FILE
	echo "ucscBrowserUrl: http://genome.ucsc.edu/cgi-bin/hgTracks?clade=insect&org=D.+melanogaster&db=dm3&position=chr[CHR]%3A[START]-[END]" >> $FILE
	echo "sortId: other" >> $FILE
fi

FILE=./genomebrowser/Gallus_gallus/Galgal4.75/genome.yaml
if [ -e $FILE ] 
then
	rm $FILE
	echo $FILE >> genomes.tmp
	echo "species: Chicken" >> $FILE
	echo "version: Galgal4.75" >> $FILE
	echo "ensemblBrowserUrl: http://www.ensembl.org/Gallus_gallus/Location/View?r=[CHR]%3A[START]-[END]" >> $FILE
	echo "ucscBrowserUrl: http://genome.ucsc.edu/cgi-bin/hgTracks?db=galGal4&position=chr[CHR]%3A[START]-[END]" >> $FILE
	echo "sortId: other" >> $FILE
fi

FILE=./genomebrowser/Gasterosteus_aculeatus/BROADS1.75/genome.yaml
if [ -e $FILE ] 
then
	rm $FILE
	echo $FILE >> genomes.tmp
	echo "species: Stickleback" >> $FILE
	echo "version: (BROADS1.75)" >> $FILE
	echo "ensemblBrowserUrl: http://www.ensembl.org/Gasterosteus_aculeatus/Location/View?r=[CHR]%3A[START]-[END]" >> $FILE
	echo "ucscBrowserUrl: https://genome.ucsc.edu/cgi-bin/hgTracks?db=gasAcu1&position=chr[CHR]%3A[START]-[END]" >> $FILE
	echo "sortId: other" >> $FILE
fi

FILE=./genomebrowser/Halorubrum_lacusprofundi_ATCC_49239/GCA_000022205.1.22/genome.yaml
if [ -e $FILE ] 
then
	rm $FILE
	echo $FILE >> genomes.tmp
	echo "species: Halorubrum_lacusprofundi_ATCC_49239" >> $FILE
	echo "version: (GCA_000022205.1.22)" >> $FILE
	echo "ensemblBrowserUrl: http://bacteria.ensembl.org/halorubrum_lacusprofundi_atcc_49239/Location/View?r=[CHR]%3A[START]-[END]" >> $FILE
	echo "ucscBrowserUrl: " >> $FILE
	echo "sortId: other" >> $FILE
fi

FILE=./genomebrowser/Homo_sapiens/GRCh37.75/genome.yaml
if [ -e $FILE ] 
then
	rm $FILE
	echo $FILE >> genomes.tmp
	echo "species: Human" >> $FILE
	echo "version: hg19 (GRCh37.75)" >> $FILE
	echo "ensemblBrowserUrl: http://www.ensembl.org/Homo_sapiens/Location/View?r=[CHR]%3A[START]-[END]" >> $FILE
	echo "ucscBrowserUrl: http://genome.ucsc.edu/cgi-bin/hgTracks?clade=mammal&org=Human&db=hg19&position=chr[CHR]%3A[START]-[END]" >> $FILE
	echo "sortId: main" >> $FILE
fi

FILE=./genomebrowser/Mus_musculus/GRCm38.75/genome.yaml
if [ -e $FILE ] 
then
	rm $FILE
	echo $FILE >> genomes.tmp
	echo "species: Mouse" >> $FILE
	echo "version: mm10 (GRCm38.75)" >> $FILE
	echo "ensemblBrowserUrl: http://www.ensembl.org/Mus_musculus/Location/View?r=[CHR]%3A[START]-[END]" >> $FILE
	echo "ucscBrowserUrl: http://genome.ucsc.edu/cgi-bin/hgTracks?clade=mammal&org=Mouse&db=mm10&position=chr[CHR]%3A[START]-[END]" >> $FILE
	echo "sortId: main" >> $FILE
fi

FILE=./genomebrowser/Ovis_aries/Oar_v3.1.75/genome.yaml
if [ -e $FILE ] 
then
	rm $FILE
	echo $FILE >> genomes.tmp
	echo "species: Sheep" >> $FILE
	echo "version: (Oar_v3.1.75)" >> $FILE
	echo "ensemblBrowserUrl: http://www.ensembl.org/Ovis_aries/Location/View?r=[CHR]%3A[START]-[END]" >> $FILE
	echo "ucscBrowserUrl: https://genome.ucsc.edu/cgi-bin/hgTracks?db=oviAri3&position=chr[CHR]%3A[START]-[END]" >> $FILE
	echo "sortId: other" >> $FILE
fi

FILE=./genomebrowser/Rattus_norvegicus/Rnor_5.0.75/genome.yaml
if [ -e $FILE ] 
then
	rm $FILE
	echo $FILE >> genomes.tmp
	echo "species: Rat" >> $FILE
	echo "version: rn5 (Rnor_5.0.75)" >> $FILE
	echo "ensemblBrowserUrl: http://www.ensembl.org/Rattus_norvegicus/Location/View?r=[CHR]%3A[START]-[END]" >> $FILE
	echo "ucscBrowserUrl: http://genome.ucsc.edu/cgi-bin/hgTracks?clade=mammal&org=Rat&db=rn5&position=chr[CHR]%3A[START]-[END]" >> $FILE
	echo "sortId: main" >> $FILE
fi

FILE=./genomebrowser/Sus_scrofa/Sscrofa10.2.75/genome.yaml
if [ -e $FILE ] 
then
	rm $FILE
	echo $FILE >> genomes.tmp
	echo "species: Pig" >> $FILE
	echo "version: (Sscrofa10.2.75)" >> $FILE
	echo "ensemblBrowserUrl: http://www.ensembl.org/Sus_scrofa/Location/View?r=[CHR]%3A[START]-[END]" >> $FILE
	echo "ucscBrowserUrl: http://genome.ucsc.edu/cgi-bin/hgTracks?clade=mammal&org=Pig&db=susScr2&position=chr[CHR]%3A[START]-[END]" >> $FILE
	echo "sortId: other" >> $FILE
fi

FILE=./genomebrowser/Schizosaccharomyces_pombe/ASM294v2.22/genome.yaml
if [ -e $FILE ] 
then
	rm $FILE
	echo $FILE >> genomes.tmp
	echo "species: S. pombe" >> $FILE
	echo "version: (ASM294v2.22)" >> $FILE
	echo "ensemblBrowserUrl: http://fungi.ensembl.org/Schizosaccharomyces_pombe/Location/View?r=[CHR]%3A[START]-[END]" >> $FILE
	echo "ucscBrowserUrl: " >> $FILE
	echo "sortId: other" >> $FILE
fi

FILE=./genomebrowser/Vitis_vinifera/IGGP_12x.22/genome.yaml
if [ -e $FILE ] 
then
	rm $FILE
	echo $FILE >> genomes.tmp
	echo "species: Vitis vinifera" >> $FILE
	echo "version: (IGGP_12x.22)" >> $FILE
	echo "ensemblBrowserUrl: http://plants.ensembl.org/Vitis_vinifera/Location/View?r=[CHR]%3A[START]-[END]" >> $FILE
	echo "ucscBrowserUrl: " >> $FILE
	echo "sortId: plant" >> $FILE
fi

cat genomes.tmp | sort > processed-genomes.tmp
find -name genome.yaml | sort > all-genomes.tmp

DIFF=$(diff processed-genomes.tmp all-genomes.tmp) 
if [ "$DIFF" != "" ] 
then
  echo ""
  echo "WARNING: following genomes were processed, but script add_gb_links.sh didn't write genome.yaml file for them"
  echo "diff processed-genomes.tmp all-genomes.tmp"
  echo $DIFF
fi


rm genomes.tmp processed-genomes.tmp all-genomes.tmp

exit 0


./Canis_familiaris/BROADD2.67/genome_v1.yaml
species: Dog
version: (BROADD2.67)
ensemblBrowserUrl: http://may2012.archive.ensembl.org/Canis_familiaris/Location/View?r=[CHR]%3A[START]-[END]
ucscBrowserUrl: http://genome.ucsc.edu/cgi-bin/hgTracks?clade=mammal&org=Dog&db=canFam2&position=chr[CHR]%3A[START]-[END]
sortId: other

./Rattus_norvegicus/RGSC3.4.69/genome_v1.yaml
species: Rat
version: rn4 (RGSC3.4.69)
ensemblBrowserUrl: http://oct2012.archive.ensembl.org/Rattus_norvegicus/Location/View?r=[CHR]%3A[START]-[END]
ucscBrowserUrl: http://genome.ucsc.edu/cgi-bin/hgTracks?clade=mammal&org=Rat&db=rn4&position=chr[CHR]%3A[START]-[END]
sortId: main


./Homo_sapiens/NCBI36.54/genome_v1.yaml
species: Human
version: hg18 (NCBI36.54)
ensemblBrowserUrl: http://may2009.archive.ensembl.org/Homo_sapiens/Location/View?db=core;r=[CHR]%3A[START]-[END]
ucscBrowserUrl: http://genome.ucsc.edu/cgi-bin/hgTracks?clade=mammal&org=Human&db=hg18&position=chr[CHR]%3A[START]-[END]
sortId: main

./Human-MT/NC_012920.1/genome_v1.yaml
species: Human mitoch.
version: (NC_012920)
ensemblBrowserUrl: 
ucscBrowserUrl: 
sortId: rest

./Mus_musculus/NCBIM37.67/genome_v1.yaml
species: Mouse
version: mm9 (NCBIM37.67)
ensemblBrowserUrl: http://may2012.archive.ensembl.org/Mus_musculus/Location/View?r=[CHR]%3A[START]-[END]
ucscBrowserUrl: http://genome.ucsc.edu/cgi-bin/hgTracks?clade=mammal&org=Mouse&db=mm9&position=chr[CHR]%3A[START]-[END]
sortId: main
./Ovis_aries/Oar_v3.1/genome_v1.yaml




