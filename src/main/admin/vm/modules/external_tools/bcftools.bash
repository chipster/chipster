##depends:none
# SAM tools, BSD License, MIT License
  cd ${TMPDIR_PATH}/
  wget_retry https://github.com/samtools/bcftools/releases/download/1.2/bcftools-1.2.tar.bz2
  tar xf bcftools-1.2.tar.bz2
  cd bcftools-1.2
  make
  cd ..
  mv bcftools-1.2/ ${TOOLS_PATH}/
  rm bcftools-1.2.tar.bz2
  #ln -s ${TOOLS_PATH}/bcftools-1.2 ${TOOLS_PATH}/bcftools