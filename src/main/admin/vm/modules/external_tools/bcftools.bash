##depends:none
# BCF tools, BSD License, MIT License
  cd ${TMPDIR_PATH}/
  #wget_retry https://github.com/samtools/bcftools/releases/download/1.2/bcftools-1.2.tar.bz2
  wget_retry http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bcftools-1.2.tar.bz2

  tar xf bcftools-1.2.tar.bz2
  cd bcftools-1.2
  make
  cd ..
  mv bcftools-1.2/ ${TOOLS_PATH}/
  rm bcftools-1.2.tar.bz2
  cd ${TOOLS_PATH}
  ln -s bcftools-1.2 bcftools