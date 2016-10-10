##depends:none

source ../installation_files/functions.bash

# Picard tools, Apache License V2.0, MIT
  cd ${TMPDIR_PATH}/
  #wget_retry -nv -O picard-tools-1.138.zip https://github.com/broadinstitute/picard/releases/download/1.138/picard-tools-1.138.zip
  #wget_retry -nv -O picard-tools-1.138.zip http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/picard-tools-1.138.zip
  wget_retry -nv https://github.com/broadinstitute/picard/releases/download/2.6.0/picard.jar
  mkdir ${TOOLS_PATH}/picard-tools-2.6.0
  mv picard.jar ${TOOLS_PATH}/picard-tools-2.6.0
  cd ${TOOLS_PATH}
  ln -s picard-tools-2.6.0 picard-tools
