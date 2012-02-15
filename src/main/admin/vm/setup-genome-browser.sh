#!/usr/bin/env bash

# Script for downloading and extracting annotation files for Chipster genome browser.
# File contents2.txt is created for class AnnotationManager to describe available data.
# Existing files won't be downloaded again, so it is easy to add some new annotations without
# downloading everything. However, any broken or partial files have to be removed manually
# to make them downloaded again. 


# download and extract .gz file from url unless the extracted file exists already
# second parameter has to match file name in url (excluding .gz ending)

download_and_extract () # 1:url 2:file
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

contents_append () # 1:type 2:chr 3:file
{
	FILE_SIZE=$(stat -c%s "$3")
	echo -e "$SPECIES\t$VERSION\t$1\t$2\t$3\t$FILE_SIZE" >> contents2.txt
}


# create url from its parts and chromosome name and process that file

download_chr () # 1:chromosome
{
	URL="$URL_PREFIX$1$URL_POSTFIX"
	FILE=$(basename $URL .gz)
	download_and_extract "$URL" "$FILE"
	contents_append "$TYPE" "$1" "$FILE"
}


# process arbitrary number of chromosomes
# numerical value is interpreted as all chromosomes between 1 and the value, 
# non-numbers are interreted as chromosomes directly

download_chrs () # 1:type 2:url-prefix 3:url-postfix 4...:chrs
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

download_and_rename () # 1:type 2:url 3:new-name
{
	FILE=$(basename $2 .gz)
	download_and_extract "$2" "$FILE"
	if [ "$3" ]
	then
		mv "$FILE" "$3"
		FILE="$3"
	fi	
	contents_append "$1" "*" "$FILE"
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

set -e # exit on errors

trap exit-trap INT TERM EXIT # execute on exit

echo "Downloading Genome Browser annotations..."

mv contents2.txt contents2.txt.old # make space for new contents
echo "CHIPSTER ANNOTATION CONTENTS FILE VERSION 2" > contents2.txt

SPECIES="Human"
VERSION="hg19 (GRCh37.65)"

download_and_rename "Cytobands" "ftp://ftp.ensembl.org/pub/release-65/mysql/homo_sapiens_core_65_37/karyotype.txt.gz" "Homo_sapiens.GRCh37.65.cytobands.txt" 
download_and_rename "Cytobands seq_region" "ftp://ftp.ensembl.org/pub/release-65/mysql/homo_sapiens_core_65_37/seq_region.txt.gz" "Homo_sapiens.GRCh37.65.seq_region.txt"
download_and_rename "ENSEMBL Gtf" "ftp://ftp.ensembl.org/pub/release-65/gtf/homo_sapiens/Homo_sapiens.GRCh37.65.gtf.gz" # no rename needed
download_chrs "Reference sequence" "ftp://ftp.ensembl.org/pub/release-65/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.65.dna.chromosome." ".fa.gz" "3" "X" "Y" 

exit 0 
