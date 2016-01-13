##depends:none

source ../installation_files/functions.bash

# BWA, GPL v3 or later, MIT License
  cd ${TMPDIR_PATH}/
  #wget_retry -O bwa-0.7.12.tar.bz2 http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2                                   
  wget_retry -O bwa-0.7.12.tar.bz2 http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bwa-0.7.12.tar.bz2                                   
  tar -xf bwa-0.7.12.tar.bz2 --use-compress-program=pbzip2
  rm -f bwa-0.7.12.tar.bz2
  cd bwa-0.7.12/
  make
  cd ../
  mv bwa-0.7.12/ ${TOOLS_PATH}/
  cd ${TOOLS_PATH}
  ln -s bwa-0.7.12 ${TOOLS_PATH}/bwa
