##depends:none
# STAR, GPLv3 License
  cd ${TMPDIR_PATH}/
  curl -sL https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz | tar -xz
  cd STAR-2.5.3a/source
  make STAR
  mkdir ${TOOLS_PATH}/STAR-2.5.3a
  mv STAR ${TOOLS_PATH}/STAR-2.5.3a
  cd ${TOOLS_PATH}
  ln -s STAR-2.5.3a STAR
  