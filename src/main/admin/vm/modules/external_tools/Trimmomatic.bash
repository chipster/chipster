##depends:none

source ../installation_files/functions.bash

# Trimmomatic, GPL
  cd ${TMPDIR_PATH}/
  wget_retry -nv http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip
  unzip Trimmomatic-0.32.zip
  mv Trimmomatic-0.32 ${TOOLS_PATH}/
  ln -s Trimmomatic-0.32 ${TOOLS_PATH}/trimmomatic
  rm Trimmomatic-0.32.zip
