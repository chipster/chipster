 ##depends:none
 
 # fseq, GPLv3
  # the original link is broken 2015-12-08
  # http://trackhubs.its.unc.edu/fureylab/fseq/fseq_1.84.tgz
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/fseq/fseq_1.84.tgz | tar -xz -C ${TMPDIR_PATH}/ 
  mv ${TMPDIR_PATH}/fseq ${TOOLS_PATH}/fseq-1.84
  ln -s fseq-1.84 ${TOOLS_PATH}/fseq
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/fseq/read_extend_bed.pm -o ${TOOLS_PATH}/fseq/bin/read_extend_bed.pm
  chmod 775 ${TOOLS_PATH}/fseq/bin/read_extend_bed.pm
