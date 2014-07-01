##depends:external_genomes/data_for_tools.bash

# Fasta files for tools

  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_mm9.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_rn4.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_athaliana.TAIR10.tar.gz | tar xz -C ${TOOLS_PATH}/
