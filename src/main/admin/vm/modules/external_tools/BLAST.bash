##depends:none

# BLAST, public domain
  #curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.29/ncbi-blast-2.2.29+-x64-linux.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/ncbi-blast-2.2.29+-x64-linux.tar.gz | tar -xz -C ${TOOLS_PATH}/
  ln -s ncbi-blast-2.2.29+ ${TOOLS_PATH}/blast    
