##depends:none

source ../installation_files/functions.bash

# HISAT2, GPLv3
  cd ${TMPDIR_PATH}/
  wget_retry -O hisat2-2.1.0-Linux_x86_64.zip ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip
 
  unzip hisat2-2.1.0-Linux_x86_64.zip -d ${TOOLS_PATH}/
  rm -f hisat2-2.1.0-Linux_x86_64.zip
  cd ${TOOLS_PATH}
  ln -s hisat2-2.1.0 hisat2
