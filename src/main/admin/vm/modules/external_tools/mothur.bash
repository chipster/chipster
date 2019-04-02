##depends:none

source ../installation_files/functions.bash

 # mothur GPLv3
  cd ${TMPDIR_PATH}/
  
  # Mothur
  # Retain version 1.39.5 as backup
  wget_retry -nv  https://github.com/mothur/mothur/releases/download/v1.39.5/Mothur.linux_64.noReadLine.zip
  unzip -q Mothur.linux_64.noReadLine.zip
  mv mothur ${TOOLS_PATH}/mothur-1.39.5
  rm -rf  __MACOSX
  
  # Make version 1.41.3 the default
  cd ${TMPDIR_PATH}/
  wget_retry -nv https://github.com/mothur/mothur/releases/download/v1.41.3/Mothur.linux_64.zip
  unzip -q Mothur.linux_64.zip
  mv mothur ${TOOLS_PATH}/mothur-1.41.3
  rm -rf  __MACOSX
  cd ${TOOLS_PATH}
  ln -s mothur-1.41.3 mothur
  
  mkdir -p ${TOOLS_PATH}/mothur-silva-reference/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/mothur/silva/v102.tar.lz4  | lz4 -d | tar -x -C ${TOOLS_PATH}/mothur-silva-reference/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/mothur/silva/v132.tar.lz4  | lz4 -d | tar -x -C ${TOOLS_PATH}/mothur-silva-reference/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/mothur/silva/silva-gold.tar.lz4  | lz4 -d | tar -x -C ${TOOLS_PATH}/mothur-silva-reference/
  
  ln -s v132 ${TOOLS_PATH}/mothur-silva-reference/silva
 
