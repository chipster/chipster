#!/bin/bash

# sudo su - chipster

# run in parallel, takes about 10 hours
time /opt/chipster/comp/modules/admin/shell/add_genome.sh -chipster_path /opt/chipster Arabidopsis_thaliana -karyotype &
time /opt/chipster/comp/modules/admin/shell/add_genome.sh -chipster_path /opt/chipster Bos_taurus -karyotype &
time /opt/chipster/comp/modules/admin/shell/add_genome.sh -chipster_path /opt/chipster Canis_familiaris -karyotype &
time /opt/chipster/comp/modules/admin/shell/add_genome.sh -chipster_path /opt/chipster Drosophila_melanogaster -karyotype &
time /opt/chipster/comp/modules/admin/shell/add_genome.sh -chipster_path /opt/chipster Gallus_gallus -karyotype &
time /opt/chipster/comp/modules/admin/shell/add_genome.sh -chipster_path /opt/chipster Gasterosteus_aculeatus -karyotype &
time /opt/chipster/comp/modules/admin/shell/add_genome.sh -chipster_path /opt/chipster Halorubrum_lacusprofundi_ATCC_49239 -karyotype &
time /opt/chipster/comp/modules/admin/shell/add_genome.sh -chipster_path /opt/chipster Homo_sapiens -karyotype &
time /opt/chipster/comp/modules/admin/shell/add_genome.sh -chipster_path /opt/chipster Mus_musculus -karyotype &
time /opt/chipster/comp/modules/admin/shell/add_genome.sh -chipster_path /opt/chipster Ovis_aries -karyotype &
time /opt/chipster/comp/modules/admin/shell/add_genome.sh -chipster_path /opt/chipster Rattus_norvegicus -karyotype &
time /opt/chipster/comp/modules/admin/shell/add_genome.sh -chipster_path /opt/chipster Schizosaccharomyces_pombe -karyotype &
time /opt/chipster/comp/modules/admin/shell/add_genome.sh -chipster_path /opt/chipster Sus_scrofa -karyotype &
time /opt/chipster/comp/modules/admin/shell/add_genome.sh -chipster_path /opt/chipster Vitis_vinifera -karyotype &

# wait all background tasks to finish
wait

# this script will tell you if you need to update genome.yaml files, so check its output
add_gb_links.sh

# package everything
package_genomes.sh

# copy to nic
