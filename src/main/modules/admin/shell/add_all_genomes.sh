#!/bin/bash

# cd into script's directory
cd ${0%/*}

add_genome () # parameters 1:genome
{
	echo "job started, writing output to /opt/chipster/tools/tmp/$1.log"
	{ time /opt/chipster/comp/modules/admin/shell/add_genome.sh -chipster_path /opt/chipster $1 -karyotype ; } &> /opt/chipster/tools/tmp/$1.log
	echo "done $1"
}

# sudo su - chipster

# Pouta medium instance doesn't have enough memory (30GB) to run everything in parallel

add_genome Arabidopsis_thaliana &
add_genome Bos_taurus &
add_genome Canis_familiaris &
add_genome Drosophila_melanogaster &
add_genome Gallus_gallus &
wait
add_genome Gasterosteus_aculeatus &
add_genome Halorubrum_lacusprofundi_ATCC_49239 &
add_genome Homo_sapiens &
add_genome Mus_musculus &
add_genome Ovis_aries &
wait
add_genome Rattus_norvegicus &
add_genome Schizosaccharomyces_pombe &
add_genome Sus_scrofa &
add_genome Vitis_vinifera &
wait

# this script will tell you if you need to update genome.yaml files, so check its output
bash add_gb_links.sh

# package everything
bash package_genomes.sh

# copy to nic
