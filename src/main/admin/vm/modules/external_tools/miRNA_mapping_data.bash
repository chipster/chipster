##depends:none

# miRNA mapping data
  cd ${TMPDIR_PATH}/
  mkdir -p ${TOOLS_PATH}/miRNA_mappings/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/miRNA_mappings/All_miRNA_mappings_v1.tar.gz | tar -xz -C ${TOOLS_PATH}/miRNA_mappings/

