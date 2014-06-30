#!/bin/bash

add_genome () # parameters 1:genome
{
	echo "started $1: 	Writing output to /opt/chipster/tools/tmp/$1.log"
	{ time /opt/chipster/comp/modules/admin/shell/add_genome.sh -chipster_path /opt/chipster $1 -karyotype ; echo "done $1" ; } &> /opt/chipster/tools/tmp/$1.log &
}

# sudo su - chipster

# run in parallel, takes about 10 hours
add_genome Arabidopsis_thaliana
#add_genome Bos_taurus -karyotype
add_genome Canis_familiaris
#add_genome Drosophila_melanogaster
#add_genome Gallus_gallus
#add_genome Gasterosteus_aculeatus
#add_genome Halorubrum_lacusprofundi_ATCC_49239
#add_genome Homo_sapiens
#add_genome Mus_musculus
add_genome Ovis_aries
add_genome Rattus_norvegicus
#add_genome Schizosaccharomyces_pombe
add_genome Sus_scrofa
#add_genome Vitis_vinifera

# wait all background tasks to finish
wait

# this script will tell you if you need to update genome.yaml files, so check its output
./add_gb_links.sh

# package everything
./package_genomes.sh

# copy to nic
