##depends:none

source ../installation_files/functions.bash

# Picard tools, Apache License V2.0, MIT
  cd ${TMPDIR_PATH}/
  #wget_retry -nv -O picard-tools-1.138.zip https://github.com/broadinstitute/picard/releases/download/1.138/picard-tools-1.138.zip
  wget_retry -nv -O picard-tools-1.138.zip http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/picard-tools-1.138.zip
  unzip -q picard-tools-1.138.zip
  rm picard-tools-1.138.zip
  # remove this optional jar because it's in the root of the zip
  #rm snappy-java-1.0.3-rc3.jar
  mv picard-tools-1.138/ ${TOOLS_PATH}
  cd ${TOOLS_PATH}
  ln -s picard-tools-1.138 picard-tools
