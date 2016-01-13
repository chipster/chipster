##depends:none

source ../installation_files/functions.bash

# Tophat 2, The Artistic License
  cd ${TMPDIR_PATH}/
  #wget_retry -O tophat-2.0.10.Linux_x86_64.tar.gz http://ccb.jhu.edu/software/tophat/downloads/tophat-2.0.10.Linux_x86_64.tar.gz
  wget_retry -O tophat-2.0.10.Linux_x86_64.tar.gz http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/tophat-2.0.10.Linux_x86_64.tar.gz
  tar -xf tophat-2.0.10.Linux_x86_64.tar.gz -C ${TOOLS_PATH}/
  rm -f tophat-2.0.10.Linux_x86_64.tar.gz
  ln -s tophat-2.0.10.Linux_x86_64 ${TOOLS_PATH}/tophat2
