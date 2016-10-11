##depends:none

source ../installation_files/functions.bash

# Drop seq tools
  cd ${TMPDIR_PATH}/
  wget_retry -nv -O Drop-seq_tools-1.12.zip http://mccarrolllab.com/download/922/
  unzip Drop-seq_tools-1.12.zip
  mv Drop-seq_tools-1.12 ${TOOLS_PATH}
  cd ${TOOLS_PATH}
  ln -s Drop-seq_tools-1.12 drop-seq_tools 
