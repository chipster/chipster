##depends:none

source ../installation_files/functions.bash

# RSeQC, GPLv3
  cd ${TMPDIR_PATH}/
  #wget_retry -nv -O RSeQC-2.3.7.tar.gz http://sourceforge.net/projects/rseqc/files/RSeQC-2.3.7.tar.gz/download
  wget_retry -nv -O RSeQC-2.3.7.tar.gz http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/RSeQC-2.3.7.tar.gz
  tar xf RSeQC-2.3.7.tar.gz
  mv RSeQC-2.3.7 ${TOOLS_PATH}/RSeQC-2.3.7
  rm -f RSeQC-2.3.7.tar.gz
  cd ${TOOLS_PATH}
  rm -f RSeQC
  ln -s RSeQC-2.3.7 RSeQC
  cd RSeQC
  python setup.py install #sudo
