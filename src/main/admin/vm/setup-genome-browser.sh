#!/usr/bin/env bash

# Script for downloading and extracting annotation files for Chipster genome browser.
# File contents2.txt is created for class AnnotationManager to describe available data.
# Existing files won't be downloaded again, so it is easy to add some new annotations without
# downloading everything. However, any broken or partial files have to be removed manually
# to get them downloaded again. 

# Species and releases are listed in the end of this script.


# download and extract .gz file from url unless the extracted file exists already
# second parameter has to match file name in url (excluding .gz ending)

download_and_extract () # parameters 1:url 2:file
{
	if [ ! -e "$2" ] # if doesn't exist
	then
		wget --no-verbose "$1"
		gzip -dq "$2.gz"	
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


# create url from its parts and chromosome name and process that file

download_chr () # parameters: 1:chromosome
{
	URL="$URL_PREFIX$1$URL_POSTFIX"
	FILE=$(basename $URL .gz)
	download_and_extract "$URL" "$FILE"
	contents_append "$TYPE" "$1" "$FILE"
}


# process arbitrary number of chromosomes
# numerical value is interpreted as all chromosomes between 1 and the value, 
# non-numbers are interpreted directly as a chromosome name

download_chrs () # parameters 1:type 2:url-prefix 3:url-postfix 4...:chrs
{
	TYPE=$1
	URL_PREFIX=$2
	URL_POSTFIX=$3

	while [ $4 ] # as long as there are more chromosome definitions
	do
		if [ "$4" -eq "$4" ] 2>/dev/null; # if the chromosome definition is a number
		then # iterate every chromosome from 1 to number
			LIMIT=$4
			for ((NCHR=1; NCHR <= LIMIT ; NCHR++))
			do
				download_chr "$NCHR"
			done
		else # just process requested chromosome
			download_chr "$4"
		fi
		shift # shift next chromosome definition to $4
	done
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

contents_append_all_chr () # parameters 1:type 2:file-name
{
	contents_append "$1" "*" "$2"
}

# process file from specified url and optionally rename the file

ensembl_mysql () # parameters 1:url 2:new-name 3:gene-name-database
{

	# These files are just downloaded
	download_and_rename "$1karyotype.txt.gz" "$2cytoband.txt" 
	contents_append_all_chr "Cytoband" "$2cytoband.txt"

	download_and_rename "$1seq_region.txt.gz" "$2seq_region.txt"
	contents_append_all_chr "Cytoband seq_region" "$2seq_region.txt"

	download_and_rename "$1coord_system.txt.gz" "$2coord_system.txt"
	contents_append_all_chr "Cytoband coord_system" "$2coord_system.txt"

	download_and_rename "$1gene.txt.gz" "$2gene.txt"
	contents_append_all_chr "Gene" "$2gene.txt"


	# xref.txt is needed to get the gene names, but it is too big for frequent read through
	# use information from external_db.txt to filter out everything useless:

	download_and_rename $1"external_db.txt.gz" $2"external_db.txt"
	download_and_rename $1"xref.txt.gz" $2"xref.txt"

	HGNC_ID=$(cat "$2external_db.txt" |grep "	$3	" |cut -f 1)
	
	# Filtering according to second tab-delimited column:
	# starts with number (one or more), followed by <tab> and value $HGNC_ID and <tab>
	# '^          [0-9]    \{1,\}                   <tab>           $HGNC_ID     <tab>'

	cat $2xref.txt | grep "^[0-9]\{1,\}	$HGNC_ID	" > "$2xref2.txt"

	rm "$2xref.txt"
	
	# external_db.txt is not needed anymore, but keep it avoid further downloads

	mv "$2xref2.txt" "$2xref.txt"

	contents_append_all_chr "Gene xref" "$2xref.txt"
}


# report success

exit-trap ()
{
	CODE=$? # keep the relevant exit code 
	if [ ! $CODE -eq "0" ] # there is probably easier way to do this
	then
		echo "Genome Browser annotations FAILED with exit code: $CODE" # how about problematic command or line?
	else
		echo "Genome Browser annotations done"
	fi
}


trap exit-trap INT TERM EXIT # execute on exit

echo "Downloading Genome Browser annotations..."

mv contents2.txt contents2.txt.old # make space for new contents

set -e # exit on errors (let contents2.txt backup above still fail, if this is the first run)

echo "CHIPSTER ANNOTATION CONTENTS FILE VERSION 2" > contents2.txt


# Species and releases

SPECIES="Human"
VERSION="hg18 (NCBI36.54)"

ensembl_mysql "ftp://ftp.ensembl.org/pub/release-54/mysql/homo_sapiens_core_54_36p/" "Homo_sapiens.NCBI36.54." "HGNC"
download_and_rename "ftp://ftp.ensembl.org/pub/release-54/gtf/homo_sapiens/Homo_sapiens.NCBI36.54.gtf.gz" # no rename needed
contents_append_all_chr "Transcript" "Homo_sapiens.NCBI36.54.gtf"
download_chrs "Reference sequence" "ftp://ftp.ensembl.org/pub/release-54/fasta/homo_sapiens/dna/Homo_sapiens.NCBI36.54.dna.chromosome." ".fa.gz" "22" "X" "Y" 

SPECIES="Human"
VERSION="hg19 (GRCh37.66)"

ensembl_mysql "ftp://ftp.ensembl.org/pub/release-66/mysql/homo_sapiens_core_66_37/" "Homo_sapiens.GRCh37.66." "HGNC"
download_and_rename "ftp://ftp.ensembl.org/pub/release-66/gtf/homo_sapiens/Homo_sapiens.GRCh37.66.gtf.gz" # no rename needed
contents_append_all_chr "Transcript" "Homo_sapiens.GRCh37.66.gtf"
download_chrs "Reference sequence" "ftp://ftp.ensembl.org/pub/release-66/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.66.dna.chromosome." ".fa.gz" "22" "X" "Y" 

SPECIES="Mouse"
VERSION="mm9 (NCBIM37.66)"

ensembl_mysql "ftp://ftp.ensembl.org/pub/release-66/mysql/mus_musculus_core_66_37/" "Mus_musculus.NCBIM37.66." "MGI"
download_and_rename "ftp://ftp.ensembl.org/pub/release-66/gtf/mus_musculus/Mus_musculus.NCBIM37.66.gtf.gz" # no rename needed
contents_append_all_chr "Transcript" "Mus_musculus.NCBIM37.66.gtf"
download_chrs "Reference sequence" "ftp://ftp.ensembl.org/pub/release-66/fasta/mus_musculus/dna/Mus_musculus.NCBIM37.66.dna.chromosome." ".fa.gz" "19" "X" "Y" 

SPECIES="Rat"
VERSION="rn4 (RGSC3.4.66)"

ensembl_mysql "ftp://ftp.ensembl.org/pub/release-66/mysql/rattus_norvegicus_core_66_34/" "Rattus_norvegicus.RGSC3.4.66." "RGD"
download_and_rename "ftp://ftp.ensembl.org/pub/release-66/gtf/rattus_norvegicus/Rattus_norvegicus.RGSC3.4.66.gtf" # no rename needed
contents_append_all_chr "Transcript" "Rattus_norvegicus.RGSC3.4.66.gtf"
download_chrs "Reference sequence" "ftp://ftp.ensembl.org/pub/release-66/fasta/rattus_norvegicus/dna/Rattus_norvegicus.RGSC3.4.66.dna.chromosome." ".fa.gz" "20" "X" 

SPECIES="Arabidopsis lyrata"
VERSION="(Araly1.13)"

ensembl_mysql "ftp://ftp.ensemblgenomes.org/pub/plants/release-13/mysql/arabidopsis_lyrata_core_13_66_10/" "Arabidopsis_lyrata.Araly1.13." "EntrezGene"
download_and_rename "ftp://ftp.ensemblgenomes.org/pub/plants/release-13/gtf/arabidopsis_lyrata/Arabidopsis_lyrata.Araly1.13.gtf" # no rename needed
contents_append_all_chr "Transcript" "Arabidopsis_lyrata.Araly1.13.gtf"
download_chrs "Reference sequence" "ftp://ftp.ensemblgenomes.org/pub/plants/release-13/fasta/arabidopsis_lyrata/dna/Arabidopsis_lyrata.Araly1.13.dna.chromosome." ".fa.gz" "8"

exit 0 
