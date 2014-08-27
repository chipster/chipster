#!/bin/bash

# exit on error
set -e

echo "Get a list of all tools..."

#bash chipster-cli.bash @credentials --yaml list-tools > tools.yaml

rm -f tool-list

cat tools.yaml | shyaml keys | while read key
do
  cat tools.yaml | shyaml keys "$key" | while read category
  do
    cat tools.yaml | shyaml get-values-0 "$key.$category" | while read -r -d $'\0' tool
    do
      echo "$tool" | shyaml get-value tool >> tool-list
      echo "" >> tool-list
    done
  done
done

cat tool-list | sort | uniq > uniq-tool-list
rm tool-list

echo "Search tools with parameter \"organism\"..."

tool_count=$(cat uniq-tool-list | wc -l)
searched=0
found=0

rm -f genome-tools

cat uniq-tool-list | while read tool
do
  bash chipster-cli.bash @credentials --yaml tool "$tool" > tool.yaml
  cat tool.yaml | shyaml get-values-0 parameters | while read -r -d $'\0' parameter
  do
    parameter=$(echo "$parameter" | shyaml get-value parameter)
    if [ "$parameter" = "organism" ]
    then
      echo $parameter >> genome-tools
      found=$(($found+1))
    fi
  done
  echo "($searched / $tool_count), found $found, $tool"
  searched=$(($searched+1))
  rm tool.yaml
done





