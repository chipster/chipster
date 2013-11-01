#!/bin/bash

# I wasn't able to get cut to work with tab delimiter, so CONVERT THEM FIRST
# (note the tab character, copy-paste will lose it, use ctrl+v and tab to type it
# on command line.
#
# cat contents2.txt | sed "s/	/,/g" > contents.csv

ls *.gtf.gz > genomes.txt
cat genomes.txt | cut -d "." -f 1 > species.txt
cat genomes.txt | cut -d "." -f 2- | sed "s/.tabix.gtf.gz//" > versions.txt
paste -d " " species.txt versions.txt > species-and-versions.txt

while read LINE;
do
        SPECIES=$(echo $LINE | cut -d " " -f 1)
        VERSION=$(echo $LINE | cut -d " " -f 2)
        echo "$SPECIES/$VERSION"
        mkdir -p "$SPECIES/$VERSION"
        mv $SPECIES.$VERSION* $SPECIES/$VERSION/

        cat contents.csv | grep ".gtf.gz.tbi" | grep $SPECIES | grep $VERSION > display.txt
        DISPLAY_SPECIES=$(cat display.txt | cut -d "," -f 1)
        DISPLAY_VERSION=$(cat display.txt | cut -d "," -f 2)

        echo "species: $DISPLAY_SPECIES" > genome_v1.yaml
        echo "version: $DISPLAY_VERSION" >> genome_v1.yaml

        ENSEMBL=$(cat contents.csv | grep "$DISPLAY_SPECIES" | grep "$DISPLAY_VERSION" | grep Ensembl | cut -d "," -f 5)
        UCSC=$(cat contents.csv | grep "$DISPLAY_SPECIES" | grep "$DISPLAY_VERSION" | grep UCSC | cut -d "," -f 5)

        echo "ensemblBrowserUrl: $ENSEMBL" >> genome_v1.yaml
        echo "ucscBrowserUrl: $UCSC" >> genome_v1.yaml
        echo "sortId: other" >> genome_v1.yaml
         
        mv genome_v1.yaml $SPECIES/$VERSION/
done <  species-and-versions.txt
