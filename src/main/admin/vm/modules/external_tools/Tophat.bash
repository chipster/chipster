##depends:none

source ../installation_files/functions.bash

# Tophat, The Artistic License
  cd ${TMPDIR_PATH}/
  #wget_retry -O tophat-1.3.2.Linux_x86_64.tar.gz http://ccb.jhu.edu/software/tophat/downloads/tophat-1.3.2.Linux_x86_64.tar.gz
  wget_retry -O tophat-1.3.2.Linux_x86_64.tar.gz http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/tophat-1.3.2.Linux_x86_64.tar.gz
  tar -xf tophat-1.3.2.Linux_x86_64.tar.gz
  rm -f tophat-1.3.2.Linux_x86_64.tar.gz
  mv tophat-1.3.2.Linux_x86_64 ${TOOLS_PATH}/
  ln -s tophat-1.3.2.Linux_x86_64 ${TOOLS_PATH}/tophat
