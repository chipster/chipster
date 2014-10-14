#!/usr/bin/env bash

# Script compares files on the disk and files defined in contents2.txt. Any differences in file names or file sizes are reported.

echo "Checking genome browser annotations..."

# Filter file names and file sizes from file contents2.txt
cat contents2.txt |grep -v CHIPSTER |grep -v Ensembl |grep -v UCSC |cut -d "	" -f 5,6 |sort > contents-files.txt

# Save file names and file sizes of all existing files
for FILE in *
do
 	FILE_SIZE=$(stat -c%s "$FILE")
	echo -e "$FILE\t$FILE_SIZE" >> dir-files.txt
done

# Don't complain about files that are there on purpose
cat dir-files.txt |grep -v "contents-files.txt" |grep -v "contents2.txt" |grep -v "gbrowser-annotation-check.sh" |sort > dir-files-sort.txt
rm dir-files.txt

# If the lists match, the check was succesfull
DIFF_ROWS=$(diff contents-files.txt dir-files-sort.txt |wc -l)
if [ $DIFF_ROWS -eq "0" ]
then
	echo "Genome browser annotations are OK"
else
	# Check failed, report differences
	echo "Genome browser annotation check FAILED"
	echo "Files missing:"
	diff contents-files.txt dir-files-sort.txt |grep "<"
	echo "Extra files:"
	diff contents-files.txt dir-files-sort.txt |grep ">"
fi

rm contents-files.txt dir-files-sort.txt

exit 0 
