##depends:none

source ../installation_files/functions.bash

# ClusterBuster, no license
  cd ${TMPDIR_PATH}/
  wget_retry -O cbust-src.tar.gz http://zlab.bu.edu/~mfrith/downloads/cbust-src.tar.gz
  tar xf cbust-src.tar.gz
  rm -rf cbust-src.tar.gz
  cd cbust-src/
  make
  mkdir ${TOOLS_PATH}/ClusterBuster/
  mv cbust ${TOOLS_PATH}/ClusterBuster/
  cd ../
  rm -rf cbust-src/
