#!/usr/bin/env bash

# Script for downloading and extracting annotation files for Chipster genome browser.
# File contents2.txt is created for class AnnotationManager to describe available data.
# Existing files won't be downloaded again, so it is easy to add some new annotations without
# downloading everything. However, any broken or partial files have to be removed manually
# to get them downloaded again. 

# Species and releases are listed in the end of this script.

# Debug tips:
# - Add "exit 0" after problematic species' urls (at the end of the file)
# - Comment out temp file removal in the problematic part of the script
# - Copy commnands from script to shell and run them one by one or investigate temp files to see where it goes wrong
# - Note: Many commands include tab-character, which is easily lost when doing copy-paste or other editing.
#   To type the tab character in shell press ctrl+v first and then tab.



##########################################################################################
# download and extract .gz file from url unless the extracted file exists already
# second parameter has to match file name in url (excluding .gz or .tar.gz ending)

download_and_extract () # parameters 1:url 2:file
{
	if [ ! -e "$2" ] # if doesn't exist
	then
		#wget --no-verbose "$1"
		wget "$1"

		if [[ $1 == *.tar.gz ]]
		then
			tar -zxf "$2.tar.gz"

			if [[ "$(ls -A genomes/fasta/nochr)" ]] #dir exists
			then							
				mv genomes/fasta/nochr/* .
				rmdir -p genomes/fasta/nochr

			elif [[ "$(ls -A genomes/fasta)" ]] #dir exists
			then							
				mv genomes/fasta/* .
				rmdir -p genomes/fasta
			fi			
		else
			gzip -dq "$2.gz"	
		fi
	else 
		echo "   Existing file $2 skipped"
	fi
}


# adds a new row to the contents file

contents_append () # parameters 1:type 2:chr 3:file
{
	FILE_SIZE=$(stat -c%s "$3")
	echo -e "$SPECIES\t$VERSION\t$1\t$2\t$3\t$FILE_SIZE" >> contents2.txt
}

# adds a new row to the contents file

contents_append_url () # parameters 1:type 2:url
{
	echo -e "$SPECIES\t$VERSION\t$1\t*\t$2\t0" >> contents2.txt
}

index_fasta () # parameters: 1:file
{
	FILE=$1
	if [ ! -e "$FILE.fai" ] # if doesn't exist
	then
		# create index
		samtools faidx $FILE
	fi
	contents_append "Reference sequence" "*" "$FILE"
	contents_append "Reference sequence index" "*" "$FILE.fai"
}

# Chipster fasta files

download_fasta () # parameters: 1:url 2:file
{
	if [ ! -e "$2" ] # if doesn't exist
	then
		ARCHIVE=$(basename $1 .tar.gz)
		download_and_extract "$1" "$ARCHIVE"

		# arhive file path
		PATHNAME=$(tar -tf "$ARCHIVE.tar.gz")
		# file name part of the path
		FILE="${PATHNAME##*/}"

		rm "$ARCHIVE.tar.gz"
	else
		FILE=$2
	fi

	index_fasta "$FILE"
}


# process file from specified url and optionally rename the file

download_and_rename () # parameters 1:url 2:new-name
{
	FILE=$(basename $1 .gz)

	if [ "$2" ]
	then
		if [ ! -e "$2" ] # if doesn't exist
		then
			download_and_extract "$1" "$FILE"
			mv "$FILE" "$2"
		else 
			echo "   Existing file $2 skipped"
		fi
		FILE="$2"
	else
		download_and_extract "$1" "$FILE"
	fi
}


process_gtf () # parameters 1:url
{
	FILE_BODY=$(basename $1 .gtf.gz)

	if [ ! -e "$FILE_BODY.tabix.gtf.gz.tbi" ] # if doesn't exist
	then	
		download_and_rename "$1" # no rename needed

		#generate list of chromosomes of genes
		#Read file  		Take only chr and name columns     	Filter out other names    	Remove duplicates   Replace useless chars with tab Or remove        And write to file
		cat "$FILE_BODY.gtf" | 	cut -f 1,9 --output-delimiter=';' | 	cut -d ';' -f 1,5      | 	uniq              | sed -e 's/; gene_name "/	/' | sed -e 's/\"//' > "$FILE_BODY.gene.tsv"

		#tabix installation folder hast to be in $PATH to find bgzip and tabix programs
		#don't exit even if grep exits with error (when there isn't any comments)
		set +e
		grep "^#" "$FILE_BODY.gtf" > "$FILE_BODY-1.gtf"
		set -e

		grep -v "^#" "$FILE_BODY.gtf" >> "$FILE_BODY-1.gtf"
		cat $FILE_BODY-1.gtf | sort -k1,1 -k4,4n -t "	" > "$FILE_BODY-sorted.gtf"		
		cat "$FILE_BODY-sorted.gtf" | bgzip > "$FILE_BODY.tabix.gtf.gz"

		rm "$FILE_BODY.gtf"
		rm "$FILE_BODY-1.gtf"
		rm "$FILE_BODY-sorted.gtf"

		#generate index
		tabix -p gff "$FILE_BODY.tabix.gtf.gz"; 

	
	else 
		echo "   Existing files $FILE_BODY.tabix.gtf* skipped"
	fi
	
	contents_append "Gene name" "*" "$FILE_BODY.gene.tsv"
	contents_append "Transcript" "*" "$FILE_BODY.tabix.gtf.gz"
	contents_append "Transcript index" "*" "$FILE_BODY.tabix.gtf.gz.tbi"	
}

# process ensembl mysql files

ensembl_mysql () # parameters 1:url 2:new-name
{

	if [ ! -e "$2repeat-tabix.bed.gz.tbi" ] || [ ! -e "$2cytoband-chr.txt" ] # if doesn't exist
	then 

		# Download database dump files
	

		download_and_rename "$1seq_region.txt.gz" "$2seq_region.txt"
		download_and_rename "$1coord_system.txt.gz" "$2coord_system.txt"

		download_and_rename "$1karyotype.txt.gz" "$2cytoband-tmp.txt" 

		download_and_rename "$1repeat_feature.txt.gz" "$2repeat_feature.txt"
		download_and_rename "$1analysis.txt.gz" "$2analysis.txt"


		# Prepare chromosome name and identifier mapping


		# search for chomosome (or group) coordinate systems (group for stickleback)
		cat "$2coord_system.txt" | grep "chromosome\|group" > coord_system-chr.txt

		# join requires sorted input
		LANG=en_EN sort -k 1 coord_system-chr.txt > coord_system-sorted.txt
		LANG=en_EN sort -k 3 "$2seq_region.txt" > seq_region-sorted.txt

		# join chromosome names and seq_region identifiers to create map of chromosome identifiers
		LANG=en_EN join -t '	' -1 1 -2 3 coord_system-sorted.txt seq_region-sorted.txt > chr_map-join.txt
	
		# remove extra columns
		cat chr_map-join.txt | cut -f 7,8 > chr_map.txt

		# join requires sorted input
		LANG=en_EN sort -k 1 chr_map.txt > chr_map-sorted.txt
	
		# clean
		rm "$2coord_system.txt" "$2seq_region.txt"
		rm coord_system-chr.txt coord_system-sorted.txt seq_region-sorted.txt chr_map-join.txt chr_map.txt

	
		# Low complexity region data


		# search for RepeatMasker analysis id, actually we should grep only from column 3
		# ignore case is required, because at least Vitis vinifera and Human have different forms
		#cat "$2analysis.txt" | grep -i "	RepeatMask	" > repeat-masker-row.txt

		# keep only the RepeatMasker id
		#cut -f 1 repeat-masker-row.txt > repeat-masker-id.txt

		# only one row in repeat-masker-id.txt, no need to sort it
		# sort the actual data
		#LANG=en_EN sort -k 9 "$2repeat_feature.txt" > repeat-sorted.txt

		# join data with RepeatMasker id to filter out data of any other analysis tools
		#LANG=en_EN join -t '	' -1 9 -2 1 repeat-sorted.txt repeat-masker-id.txt > repeat-masker-join.txt

		# remove extra columns
		#cat repeat-masker-join.txt | cut -f 3,4,5 > repeat-masker.txt		
		cat "$2repeat_feature.txt" | cut -d '	' -f 2,3,4 > repeat-masker.txt
		

		# join requires sorted input
		LANG=en_EN sort -k 1 repeat-masker.txt > repeat-masker-sorted.txt

		# join chromosome identifiers and the data
		LANG=en_EN join -t '	' -1 1 -2 1 chr_map-sorted.txt repeat-masker-sorted.txt > repeat-join.txt

		# remove extra columns
		# FIXME Ensembl uses 1-based coordinates, whereas standard bed file must have 0-based coordinates
		cat repeat-join.txt | cut -d '	' -f 2,3,4 > repeat.bed
	
		# bed to tabix
		cat repeat.bed | sort -k1,1 -k2,2n > repeat-sorted.bed		
		cat repeat-sorted.bed | bgzip > "$2repeat-tabix.bed.gz"

		#generate index
		tabix -p bed "$2repeat-tabix.bed.gz"

		# clean
		rm "$2analysis.txt" "$2repeat_feature.txt"
		#rm repeat-masker-row.txt repeat-masker-id.txt repeat-sorted.txt repeat-masker-join.txt
		rm repeat-masker.txt  repeat-masker-sorted.txt repeat-join.txt
		rm repeat.bed repeat-sorted.bed


		# Cytoband data


		# sort the actual data
		LANG=en_EN sort -k 2 "$2cytoband-tmp.txt" > cytoband-sorted.txt

		# join chromosome identifiers and the data
		LANG=en_EN join -t '	' -1 1 -2 2 chr_map-sorted.txt cytoband-sorted.txt > cytoband-join.txt

		# remove extra columns
		cat cytoband-join.txt | cut -d '	' -f 2,3,4,5,6,7 > "$2cytoband-chr.txt"

		# clean
		rm cytoband-sorted.txt cytoband-join.txt "$2cytoband-tmp.txt"
		rm chr_map-sorted.txt

	else 
		echo "   Existing cytoband and repeat files ($2) skipped"
	fi

	contents_append "Cytoband" "*" "$2cytoband-chr.txt"
	contents_append "Repeat" "*" "$2repeat-tabix.bed.gz"
	contents_append "Repeat index" "*" "$2repeat-tabix.bed.gz.tbi"
}


# report success

exit-trap ()
{
	CODE=$? # keep the relevant exit code 
	if [ ! $CODE -eq "0" ] # there is probably easier way to do this
	then
		echo "Genome Browser annotations FAILED with exit code: $CODE" # how to report problematic command or script line?
	else
		echo "Genome Browser annotations done"
	fi
}


trap exit-trap INT TERM EXIT # execute on exit

echo "Downloading Genome Browser annotations..."

mv contents2.txt contents2.txt.old # make space for new contents

set -e # exit on errors (let contents2.txt backup above still fail, if this is the first run)

echo "CHIPSTER ANNOTATION CONTENTS FILE VERSION 2" > contents2.txt

##########################################################################################

# Species and releases
# Humans
SPECIES="Human"
VERSION="hg19 (GRCh37.70)"

process_gtf "ftp://ftp.ensembl.org/pub/release-70/gtf/homo_sapiens/Homo_sapiens.GRCh37.70.gtf.gz"
ensembl_mysql "ftp://ftp.ensembl.org/pub/release-70/mysql/homo_sapiens_core_70_37/" "Homo_sapiens.GRCh37.70."
download_fasta "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_hg19.tar.gz" "hg19.fa"
contents_append_url "Ensembl" "http://www.ensembl.org/Homo_sapiens/Location/View?r=[CHR]%3A[START]-[END]"
contents_append_url "UCSC" "http://genome.ucsc.edu/cgi-bin/hgTracks?clade=mammal&org=Human&db=hg19&position=chr[CHR]%3A[START]-[END]"

SPECIES="Human"
VERSION="hg18 (NCBI36.54)"

process_gtf "ftp://ftp.ensembl.org/pub/release-54/gtf/homo_sapiens/Homo_sapiens.NCBI36.54.gtf.gz"
ensembl_mysql "ftp://ftp.ensembl.org/pub/release-54/mysql/homo_sapiens_core_54_36p/" "Homo_sapiens.NCBI36.54."
download_fasta "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Homo_sapiens.NCBI36.54-v2.tar.gz" "Homo_sapiens.NCBI36.54.dna.toplevel.fa"
contents_append_url "Ensembl" "http://may2009.archive.ensembl.org/Homo_sapiens/Location/View?db=core;r=[CHR]%3A[START]-[END]"
contents_append_url "UCSC" "http://genome.ucsc.edu/cgi-bin/hgTracks?clade=mammal&org=Human&db=hg18&position=chr[CHR]%3A[START]-[END]"


# Regularly used animals
SPECIES="Mouse"
VERSION="mm10 (GRCm38.70)"

process_gtf "ftp://ftp.ensembl.org/pub/release-70/gtf/mus_musculus/Mus_musculus.GRCm38.70.gtf.gz"
ensembl_mysql "ftp://ftp.ensembl.org/pub/release-70/mysql/mus_musculus_core_70_38/" "Mus_musculus.GRCm38.70."
download_fasta "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_mm10.tar.gz" "mm10.fa"
contents_append_url "Ensembl" "http://www.ensembl.org/Mus_musculus/Location/View?r=[CHR]%3A[START]-[END]"
contents_append_url "UCSC" "http://genome.ucsc.edu/cgi-bin/hgTracks?clade=mammal&org=Mouse&db=mm10&position=chr[CHR]%3A[START]-[END]"

SPECIES="Mouse"
VERSION="mm9 (NCBIM37.67)"

ensembl_mysql "ftp://ftp.ensembl.org/pub/release-67/mysql/mus_musculus_core_67_37/" "Mus_musculus.NCBIM37.67."
process_gtf "ftp://ftp.ensembl.org/pub/release-67/gtf/mus_musculus/Mus_musculus.NCBIM37.67.gtf.gz"
download_fasta "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_mm9.tar.gz" "mm9.fa"
contents_append_url "Ensembl" "http://may2012.archive.ensembl.org/Mus_musculus/Location/View?r=[CHR]%3A[START]-[END]"
contents_append_url "UCSC" "http://genome.ucsc.edu/cgi-bin/hgTracks?clade=mammal&org=Mouse&db=mm9&position=chr[CHR]%3A[START]-[END]"

SPECIES="Rat"
VERSION="rn5 (Rnor_5.0.70)"

process_gtf "ftp://ftp.ensembl.org/pub/release-70/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_5.0.70.gtf.gz"
ensembl_mysql "ftp://ftp.ensembl.org/pub/release-70/mysql/rattus_norvegicus_core_70_5/" "Rattus_norvegicus.Rnor_5.0.70."
download_fasta "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Rattus_norvegicus.Rnor_5.0.70.tar.gz" "Rattus_norvegicus.Rnor_5.0.70.dna.toplevel.fa"
contents_append_url "Ensembl" "http://www.ensembl.org/Rattus_norvegicus/Location/View?r=[CHR]%3A[START]-[END]"
contents_append_url "UCSC" "http://genome.ucsc.edu/cgi-bin/hgTracks?clade=mammal&org=Rat&db=rn5&position=chr[CHR]%3A[START]-[END]"

SPECIES="Rat"
VERSION="rn4 (RGSC3.4.69)"

process_gtf "ftp://ftp.ensembl.org/pub/release-69/gtf/rattus_norvegicus/Rattus_norvegicus.RGSC3.4.69.gtf.gz"
ensembl_mysql "ftp://ftp.ensembl.org/pub/release-69/mysql/rattus_norvegicus_core_69_34/" "Rattus_norvegicus.RGSC3.4.69."
download_fasta "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_rn4.tar.gz" "rn4.fa"
contents_append_url "Ensembl" "http://oct2012.archive.ensembl.org/Rattus_norvegicus/Location/View?r=[CHR]%3A[START]-[END]"
contents_append_url "UCSC" "http://genome.ucsc.edu/cgi-bin/hgTracks?clade=mammal&org=Rat&db=rn4&position=chr[CHR]%3A[START]-[END]"


# Other animals in alphabetical order


SPECIES="Dog"
VERSION="(CanFam3.1.70)"

process_gtf "ftp://ftp.ensembl.org/pub/release-70/gtf/canis_familiaris/Canis_familiaris.CanFam3.1.70.gtf.gz"
ensembl_mysql "ftp://ftp.ensembl.org/pub/release-70/mysql/canis_familiaris_core_70_31/" "Canis_familiaris.CanFam3.1.70."
download_fasta "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Canis_familiaris.CanFam3.1.70.tar.gz" "Canis_familiaris.CanFam3.1.70.dna.toplevel.fa"
contents_append_url "Ensembl" "http://www.ensembl.org/Canis_familiaris/Location/View?r=[CHR]%3A[START]-[END]"
contents_append_url "UCSC" "http://genome.ucsc.edu/cgi-bin/hgTracks?clade=mammal&org=Dog&db=canFam3&position=chr[CHR]%3A[START]-[END]"

SPECIES="Dog"
VERSION="(BROADD2.67)"

ensembl_mysql "ftp://ftp.ensembl.org/pub/release-67/mysql/canis_familiaris_core_67_2/" "Canis_familiaris.BROADD2.67."
process_gtf "ftp://ftp.ensembl.org/pub/release-67/gtf/canis_familiaris/Canis_familiaris.BROADD2.67.gtf.gz"
download_fasta "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Canis_familiaris.BROADD2.67-v2.tar.gz" "Canis_familiaris.BROADD2.67.dna.toplevel.fa"
contents_append_url "Ensembl" "http://may2012.archive.ensembl.org/Canis_familiaris/Location/View?r=[CHR]%3A[START]-[END]"
contents_append_url "UCSC" "http://genome.ucsc.edu/cgi-bin/hgTracks?clade=mammal&org=Dog&db=canFam2&position=chr[CHR]%3A[START]-[END]"

SPECIES="Chicken"
VERSION="(Gallus_gallus-4.0)"

process_gtf "ftp://ftp.ensembl.org/pub/pre/gtf/gallus_gallus/Gallus_gallus.Gallus_gallus-4.0.pre.gtf.gz"
download_fasta "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Gallus_gallus.Gallus_gallus-4.0.pre-v2.tar.gz" "Gallus_gallus.Gallus_gallus-4.0.pre.fa"
contents_append_url "Ensembl" ""
contents_append_url "UCSC" "" 

SPECIES="Cow"
VERSION="(UMD3.1.70)"

process_gtf "ftp://ftp.ensembl.org/pub/release-70/gtf/bos_taurus/Bos_taurus.UMD3.1.70.gtf.gz"
ensembl_mysql "ftp://ftp.ensembl.org/pub/release-70/mysql/bos_taurus_core_70_31/" "Bos_taurus.UMD3.1.70."
download_fasta "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Bos_taurus.UMD3.1.70.tar.gz" "Bos_taurus.UMD3.1.70.dna.toplevel.fa"
contents_append_url "Ensembl" "http://www.ensembl.org/Bos_taurus/Location/View?r=[CHR]%3A[START]-[END]"
contents_append_url "UCSC" "http://genome.ucsc.edu/cgi-bin/hgTracks?clade=mammal&org=Cow&db=bosTau6&position=chr[CHR]%3A[START]-[END]"

SPECIES="Fruitfly"
VERSION="(BDGP5.70)"

process_gtf "ftp://ftp.ensembl.org/pub/release-70/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP5.70.gtf.gz"
ensembl_mysql "ftp://ftp.ensembl.org/pub/release-70/mysql/drosophila_melanogaster_core_70_546/" "Drosophila_melanogaster.BDGP5.70."
download_fasta "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Drosophila_melanogaster.BDGP5.70.tar.gz" "Drosophila_melanogaster.BDGP5.70.dna.toplevel.fa"
contents_append_url "Ensembl" "http://www.ensembl.org/Drosophila_melanogaster/Location/View?r=[CHR]%3A[START]-[END]"
contents_append_url "UCSC" "http://genome.ucsc.edu/cgi-bin/hgTracks?clade=insect&org=D.+melanogaster&db=dm3&position=chr[CHR]%3A[START]-[END]"

SPECIES="Sheep"
VERSION="(Oar_v3.1)"

download_fasta "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_ovis_aries_texel.tar.gz" "ovis_aries_texel.fa"
contents_append_url "Ensembl" ""
contents_append_url "UCSC" ""

SPECIES="Stickleback"
VERSION="(BROADS1.70)"

process_gtf "ftp://ftp.ensembl.org/pub/release-70/gtf/gasterosteus_aculeatus/Gasterosteus_aculeatus.BROADS1.70.gtf.gz"
ensembl_mysql "ftp://ftp.ensembl.org/pub/release-70/mysql/gasterosteus_aculeatus_core_70_1/" "Gasterosteus_aculeatus.BROADS1.70."
download_fasta "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Gasterosteus_aculeatus.BROADS1.70.tar.gz" "Gasterosteus_aculeatus.BROADS1.70.dna.toplevel.fa"
contents_append_url "Ensembl" "http://www.ensembl.org/Gasterosteus_aculeatus/Location/View?r=[CHR]%3A[START]-[END]"
contents_append_url "UCSC" ""

SPECIES="Pig"
VERSION="(Sscrofa10.2.70)"

process_gtf "ftp://ftp.ensembl.org/pub/release-70/gtf/sus_scrofa/Sus_scrofa.Sscrofa10.2.70.gtf.gz"
ensembl_mysql "ftp://ftp.ensembl.org/pub/release-70/mysql/sus_scrofa_core_70_102/" "Sus_scrofa.Sscrofa10.2.70."
download_fasta "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Sus_scrofa.Sscrofa10.2.69.tar.gz" "Sus_scrofa.Sscrofa10.2.69.dna.toplevel.fa"
contents_append_url "Ensembl" "http://www.ensembl.org/Sus_scrofa/Location/View?r=[CHR]%3A[START]-[END]"
contents_append_url "UCSC" "http://genome.ucsc.edu/cgi-bin/hgTracks?clade=mammal&org=Pig&db=susScr2&position=chr[CHR]%3A[START]-[END]"

# Plants in alphabetical order
SPECIES="Arabidopsis lyrata"
VERSION="(v.1.0.17)"

process_gtf "ftp://ftp.ensemblgenomes.org/pub/plants/release-17/gtf/arabidopsis_lyrata/Arabidopsis_lyrata.v.1.0.17.gtf.gz"
ensembl_mysql "ftp://ftp.ensemblgenomes.org/pub/plants/release-17/mysql/arabidopsis_lyrata_core_17_70_10/" "Arabidopsis_lyrata.v.1.0.17."
download_fasta "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_Arabidopsis_lyrata.v.1.0.16.tar.gz" "Arabidopsis_lyrata.v.1.0.16.fa"
contents_append_url "Ensembl" "http://plants.ensembl.org/Arabidopsis_lyrata/Location/View?db=core;r=[CHR]%3A[START]-[END]"
contents_append_url "UCSC" ""

SPECIES="Arabidopsis thaliana"
VERSION="(TAIR10.17)"

process_gtf "ftp://ftp.ensemblgenomes.org/pub/plants/release-17/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.17.gtf.gz"
ensembl_mysql "ftp://ftp.ensemblgenomes.org/pub/plants/release-17/mysql/arabidopsis_thaliana_core_17_70_10/" "Arabidopsis_thaliana.TAIR10.17."
download_fasta "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Arabidopsis_thaliana.TAIR10.17.tar.gz" "Arabidopsis_thaliana.TAIR10.17.dna.toplevel.fa"
contents_append_url "Ensembl" "http://plants.ensembl.org/Arabidopsis_thaliana/Location/View?r=[CHR]%3A[START]-[END]"
contents_append_url "UCSC" ""

SPECIES="Vitis vinifera"
VERSION="(IGGP_12x.17)"

process_gtf "ftp://ftp.ensemblgenomes.org/pub/plants/release-17/gtf/vitis_vinifera/Vitis_vinifera.IGGP_12x.17.gtf.gz"
ensembl_mysql "ftp://ftp.ensemblgenomes.org/pub/plants/release-17/mysql/vitis_vinifera_core_17_70_3/" "Vitis_vinifera.IGGP_12x.17."
download_fasta "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Vitis_vinifera.IGGP_12x.17.tar.gz" "Vitis_vinifera.IGGP_12x.17.dna.toplevel.fa"
contents_append_url "Ensembl" "http://plants.ensembl.org/Vitis_vinifera/Location/View?r=[CHR]%3A[START]-[END]"
contents_append_url "UCSC" ""

# Other
SPECIES="Yersinia similis"
VERSION="(N916Ysi)"

process_gtf "http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/setup/N916Ysi.gtf.gz"
download_fasta "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_N916Ysi.tar.gz" "N916Ysi.fa"
contents_append_url "Ensembl" ""
contents_append_url "UCSC" ""

SPECIES="Y. phage phiR1-RT"
VERSION="(HE956709.1)"

process_gtf "http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/setup/R1-RT.gtf.gz"
download_fasta "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_R1-RT.tar.gz" "R1-RT.fa"
contents_append_url "Ensembl" ""
contents_append_url "UCSC" ""

SPECIES="Human mitoch."
VERSION="(NC_012920)"

process_gtf "http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/setup/Human-MT.NC_012920.1.gtf.gz"
download_fasta "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Human-MT.NC_012920.1-v2.tar.gz" "Human-MT.NC_012920.1.fa"
contents_append_url "Ensembl" ""
contents_append_url "UCSC" ""

exit 0 
