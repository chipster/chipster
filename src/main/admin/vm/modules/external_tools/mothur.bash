##depends:none

source ../installation_files/functions.bash

 # mothur GPLv3
  cd ${TMPDIR_PATH}/
  wget_retry -nv http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/mothur/Mothur-1.28.cen_64.noReadLine.zip
  unzip -q Mothur-1.28.cen_64.noReadLine.zip
  mv mothur ${TOOLS_PATH}/mothur-1.28
  ln -s mothur-1.28 ${TOOLS_PATH}/mothur
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/mothur/mothur-data.tar.gz | tar -xz -C ${TOOLS_PATH}/

