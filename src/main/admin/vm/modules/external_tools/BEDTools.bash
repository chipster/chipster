##depends:none

source ../installation_files/functions.bash

# BEDTools, GNU GPL v2
  cd ${TMPDIR_PATH}/
  #wget_retry -nv https://github.com/arq5x/bedtools2/archive/v2.25.0.tar.gz
  wget_retry -nv http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bedtools_v2.25.0.tar.gz

  tar xf bedtools_v2.25.0.tar.gz
  cd bedtools2-2.25.0
    make clean
  make all
  cd ../
  mv bedtools2-2.25.0 ${TOOLS_PATH}/
  rm bedtools_v2.25.0.tar.gz
  cd ${TOOLS_PATH}
  ln -s bedtools2-2.25.0 ${TOOLS_PATH}/bedtools
  

  #curl -L http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bedtools-2.17.0-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/

