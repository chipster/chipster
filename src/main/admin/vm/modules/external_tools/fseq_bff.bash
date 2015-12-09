##depends:none

source ../installation_files/functions.bash

# Fseq bff files
  cd ${TMPDIR_PATH}/
  
  #wget_retry -nv -O unique20bp_hg19.tgz http://fureylab.med.unc.edu/fseq/unique20bp_hg19.tgz
  #wget_retry -nv -O unique35bp_hg19.tgz http://fureylab.med.unc.edu/fseq/unique35bp_hg19.tgz
  wget_retry -nv -O unique20bp_hg19.tgz http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/fseq_bff/unique20bp_hg19.tgz
  wget_retry -nv -O unique35bp_hg19.tgz http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/fseq_bff/unique35bp_hg19.tgz
  
  mkdir -p unique20bp_hg19
  mkdir -p unique35bp_hg19
  
  tar zxf unique20bp_hg19.tgz -C unique20bp_hg19
  tar zxf unique35bp_hg19.tgz -C unique35bp_hg19
  
  mkdir -p ${TOOLS_PATH}/fseq_bff
  
  mv unique20bp_hg19/bff_20 ${TOOLS_PATH}/fseq_bff/unique20bp_hg19
  mv unique35bp_hg19/bff_35 ${TOOLS_PATH}/fseq_bff/unique35bp_hg19
   
  rm -f unique20bp_hg19.tgz
  rm -f unique35bp_hg19.tgz
  
  rmdir unique20bp_hg19
  rmdir unique35bp_hg19
