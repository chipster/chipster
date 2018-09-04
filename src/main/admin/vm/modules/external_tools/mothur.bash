##depends:none

source ../installation_files/functions.bash

 # mothur GPLv3
  cd ${TMPDIR_PATH}/
  
  # Mothur
  wget_retry -nv  https://github.com/mothur/mothur/releases/download/v1.39.5/Mothur.linux_64.noReadLine.zip
  unzip -q Mothur.linux_64.noReadLine.zip
  mv mothur ${TOOLS_PATH}/mothur-1.39.5
  ln -s ${TOOLS_PATH}/mothur-1.39.5 ${TOOLS_PATH}/mothur
  rm -rf  __MACOSX
  
  wget_retry -nv  http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/mothur/Mothur-1.40.5.zip
  unzip -q Mothur-1.40.5.zip
  mv mothur ${TOOLS_PATH}/mothur-1.40.5

  mkdir -p ${TOOLS_PATH}/mothur-silva-reference/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/mothur/silva/v102.tar.lz4  | lz4 -d | tar -x -C ${TOOLS_PATH}/mothur-silva-reference/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/mothur/silva/v132.tar.lz4  | lz4 -d | tar -x -C ${TOOLS_PATH}/mothur-silva-reference/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/mothur/silva/silva-gold.tar.lz4  | lz4 -d | tar -x -C ${TOOLS_PATH}/mothur-silva-reference/
  
  ln -s v132 ${TOOLS_PATH}/mothur-silva-reference/silva
 
