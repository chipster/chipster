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
bash chipster-cli.bash @credentials --yaml tool "$tool" > tool.yaml

rm -f genomes

# iterate through a sequence of parameters
cat tool.yaml | shyaml get-values-0 parameters | \
  while read -r -d $'\0' parameter
  do
    # each parameter is a yaml struct stored in variable $parameter
    # get the parameter id (stored with key "parameter")
    id=$(echo "$parameter" | shyaml get-value parameter)
    # if this is a parameter with id "organism"
    if [ $id = "organism" ]
    then
      # iterate through a sequence of options
      echo "$parameter" | shyaml get-values-0 options | \
        while read -r -d $'\0' option_struct
        do
	  # each option is a yaml struct stored in variable $option_struct
          # get the optio id (stored with key "option")
          genome=$(echo "$option_struct" | shyaml get-value option)
          if [ $genome != "other" ]
          then
            echo $genome >> genomes
          fi
        done
    fi
  done

rm tool.yaml

bash chipster-cli.bash @credentials --quiet import $input

while read genome
do
    echo "Run $tool for genome $genome..."    
    bash chipster-cli.bash @credentials --quiet run "$tool" --dataset "$input" --parameter organism="$genome"
done < genomes
mv genomes $tool.genomes

echo "Save session $session..."
bash chipster-cli.bash @credentials save-session "$session"