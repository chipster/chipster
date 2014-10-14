##depends:none

source ../installation_files/functions.bash

# Picard tools, Apache License V2.0, MIT
  cd ${TMPDIR_PATH}/
  wget_retry -nv -O picard-tools-1.105.zip http://sourceforge.net/projects/picard/files/picard-tools/1.105/picard-tools-1.105.zip/download
  unzip -q picard-tools-1.105.zip
  rm picard-tools-1.105.zip
  # remove this optional jar because it's in the root of the zip
  rm snappy-java-1.0.3-rc3.jar
  mv picard-tools-1.105/ ${TOOLS_PATH}
  cd ${TOOLS_PATH}
  ln -s picard-tools-1.105 picard-tools
