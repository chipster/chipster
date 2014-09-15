#!/bin/bash

# exit on error
set -e

if [[ -z $3 ]]
then 
    echo "usage: run-tool-for-all-genomes.bash input-file tool-id output-session"
    exit 1
fi

input="$1"
tool="$2"
session="$3"

echo "Clear session..."
bash chipster-cli.bash @credentials clear-session

echo "List genomes..."
bash get-genomes.bash "$tool" genomes


bash chipster-cli.bash @credentials --quiet import $input

while read genome
do
    echo "Run $tool for genome $genome..."    
    bash chipster-cli.bash @credentials --quiet run "$tool" --dataset "$input" --parameter organism="$genome"
done < genomes
mv genomes $tool.genomes

echo "Save session $session..."
bash chipster-cli.bash @credentials save-session "$session"