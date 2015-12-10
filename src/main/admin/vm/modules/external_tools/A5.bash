##depends:none

source ../installation_files/functions.bash

# Fseq bff files
  cd ${TMPDIR_PATH}/
  
  #wget_retry -nv -O a5_miseq_linux_20140604.tar.gz http://downloads.sourceforge.net/project/ngopt/a5_miseq_linux_20140604.tar.gz
  wget_retry -nv -O a5_miseq_linux_20140604.tar.gz http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/a5_miseq_linux_20140604.tar.gz
  
  tar zxvf a5_miseq_linux_20140604.tar.gz
  
  mv a5_miseq_linux_20140604 ${TOOLS_PATH}/a5_miseq
   
  rm -f a5_miseq_linux_20140604.tar.gz
  
