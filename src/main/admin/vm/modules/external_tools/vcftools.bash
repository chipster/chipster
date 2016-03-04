##depends:none

source ../installation_files/functions.bash


# vcftools, GPLv3
  cd ${TMPDIR_PATH}/
  #wget_retry  https://github.com/vcftools/vcftools/releases/download/v0.1.14/vcftools-0.1.14.tar.gz
  wget_retry  http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/vcftools-0.1.14.tar.gz
  tar xf vcftools-0.1.14.tar.gz
  cd vcftools-0.1.14/
  ./configure --prefix=$PWD
  make
  make install
  cd ../
  rm -f vcftools-0.1.14.tar.gz
  mv vcftools-0.1.14/ ${TOOLS_PATH}/
  cd ${TOOLS_PATH}
  ln -s vcftools-0.1.14 vcftools
