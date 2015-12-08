##depends:none

source ../installation_files/functions.bash

# tagcleaner, GPLv3
  cd ${TMPDIR_PATH}/
  #wget_retry -nv -O tagcleaner-standalone-0.12.tar.gz http://downloads.sourceforge.net/project/tagcleaner/standalone/tagcleaner-standalone-0.12.tar.gz
  wget_retry -nv -O tagcleaner-standalone-0.12.tar.gz http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/tagcleaner-standalone-0.12.tar.gz
  tar xf tagcleaner-standalone-0.12.tar.gz -C ${TOOLS_PATH}/
  rm -f tagcleaner-standalone-0.12.tar.gz
  ln -s tagcleaner-standalone-0.12 ${TOOLS_PATH}/tagcleaner
