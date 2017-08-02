##depends:none

source ../installation_files/functions.bash

#  TODO: is this TMPDIR_PATH the same as tools?
  cd ${TMPDIR_PATH}/
  wget_retry -nv -0 stringtie-1.3.3b.Linux_x86_64.tar.gz  http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.3b.Linux_x86_64.tar.gz
  tar -xzf stringtie-1.3.3b.Linux_x86_64.tar.gz
#  mv stringtie-1.3.3b.Linux_x86_64/ ${TOOLS_PATH}
  ln -s stringtie-1.3.3b.Linux_x86_64 ${TOOLS_PATH}/stringtie
  # Clean up
#  rm stringtie-1.3.3b.Linux_x86_64.tar.gz
  rm -rf stringtie-1.3.3b.Linux_x86_64
