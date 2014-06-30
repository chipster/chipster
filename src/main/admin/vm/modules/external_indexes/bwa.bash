##depends:none

# bwa indexes, built for Chipster
  cd ${TMPDIR_PATH}/
  mkdir -p ${TOOLS_PATH}/bwa_indexes/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_mm9.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_rn4.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_Rattus_norvegicus.Rnor_5.0.70.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_athaliana.TAIR10.tar.gz | tar xz -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_Canis_familiaris.CanFam3.1.71.tar.gz | tar xz -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_Gasterosteus_aculeatus.BROADS1.71.tar.gz | tar xz -C ${TOOLS_PATH}/

