##depends:none

# prinseq
  cd ${TMPDIR_PATH}/
  curl -sL http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz/download | tar -xz
  chmod a+x prinseq-lite-0.20.4/*.pl
  mv prinseq-lite-0.20.4 ${TOOLS_PATH}/
  ln -s prinseq-lite-0.20.4 ${TOOLS_PATH}/prinseq

