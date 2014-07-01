#!/bin/bash

set -e

TMPDIR=$1

#cd /opt/chipster/tools/genomes
cd $TMPDIR

find ../tmp/ -name *.files > packages.tmp

while read FILE_LIST           
do           
	PACKAGE=$(basename $FILE_LIST .files)
	echo "Packaging $PACKAGE.tar.gz"
	# multi-core
	tar -I pigz -T $FILE_LIST -cf $PACKAGE.tar.gz
	# single-core
	# tar -T $FILE_LIST -czf $PACKAGE.tar.gz
done < packages.tmp 

rm packages.tmp

echo "Packages done."
