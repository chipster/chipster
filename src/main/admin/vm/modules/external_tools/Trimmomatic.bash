##depends:none

source ../installation_files/functions.bash

# Trimmomatic, GPL
  cd ${TMPDIR_PATH}/
  #wget_retry -nv http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip
  wget_retry -nv http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/Trimmomatic-0.33.zip
  unzip Trimmomatic-0.33.zip
  mv Trimmomatic-0.33 ${TOOLS_PATH}/
  rm Trimmomatic-0.33.zip
  cd ${TOOLS_PATH}
  ln -s Trimmomatic-0.33 ${TOOLS_PATH}/trimmomatic
  
