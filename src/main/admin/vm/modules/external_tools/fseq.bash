 ##depends:none
 
 # fseq, GPLv3
  curl -s http://fureylab.med.unc.edu/fseq/fseq_1.84.tgz | tar -xz -C ${TMPDIR_PATH}/ 
  mv ${TMPDIR_PATH}/fseq ${TOOLS_PATH}/fseq-1.84
  ln -s fseq-1.84 ${TOOLS_PATH}/fseq
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/misc/read_extend_bed.pm -o ${TOOLS_PATH}/fseq/bin/read_extend_bed.pm
  chmod 775 ${TOOLS_PATH}/fseq/bin/read_extend_bed.pm
