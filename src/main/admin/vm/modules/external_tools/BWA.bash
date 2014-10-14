##depends:none

source ../installation_files/functions.bash

# BWA, GPL v3 or later, MIT License
  cd ${TMPDIR_PATH}/
  wget_retry -O bwa-0.7.10.tar.bz2 http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.10.tar.bz2/download
  tar -xf bwa-0.7.10.tar.bz2 --use-compress-program=pbzip2
  rm -f bwa-0.7.10.tar.bz2
  cd bwa-0.7.10/
  make
  cd ../
  mv bwa-0.7.10/ ${TOOLS_PATH}/
  ln -s bwa-0.7.10 ${TOOLS_PATH}/bwa
