#!/bin/bash

set -e

cd /opt/chipster/tools/genomes

rm -f install-chipster-commands.sh

find ../tmp/ -name *.files > packages.tmp

while read FILE_LIST           
do           
	PACKAGE=$(basename $FILE_LIST .files)
	echo "Packaging $PACKAGE.tar.gz"
  	echo '  curl -s http://${NIC_MIRROR}/pub/sci/molbio/chipster/dist/tools_extras/genomes/3.0.0/'${PACKAGE}'.tar.gz | tar -xz -C ${TOOLS_PATH}/genomes/' >> install-chipster-commands.sh
	# multi-core
	tar -I pigz -T $FILE_LIST -cf $PACKAGE.tar.gz
	# single-core
	# tar -T $FILE_LIST -czf $PACKAGE.tar.gz
done < packages.tmp 

rm packages.tmp

echo "Packaging done. Installation commands in file install-chipster-commands.sh"
