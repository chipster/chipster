#!/bin/bash
# Command 'onevm list' cuts template names too short. This script shows the full names alongside the original output
onevm list > list.txt
cat list.txt | cut -d " " -f 2 > id.txt
while read ID; do
	onevm show $ID | grep "NAME                :" >> name.txt
done < id.txt

echo "" > name2.txt
cat name.txt | cut -d ":" -f 2 >> name2.txt

paste list.txt name2.txt > paste.txt
cat paste.txt

rm list.txt id.txt name.txt name2.txt paste.txt