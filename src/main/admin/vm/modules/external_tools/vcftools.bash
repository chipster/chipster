##depends:none

source ../installation_files/functions.bash


# vcftools, GPLv3
  cd ${TMPDIR_PATH}/
  wget_retry -O vcftools_0.1.11.tar.gz http://sourceforge.net/projects/vcftools/files/vcftools_0.1.11.tar.gz/download
  tar xf vcftools_0.1.11.tar.gz
  cd vcftools_0.1.11/
  make
  cd ../
  rm -f vcftools_0.1.11.tar.gz
  mv vcftools_0.1.11/ ${TOOLS_PATH}/
  #curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/vcftools_0.1.11-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/
  ln -s vcftools_0.1.11 ${TOOLS_PATH}/vcftools
