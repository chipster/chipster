##depends:none

# prinseq
  cd ${TMPDIR_PATH}/
  #curl -sL http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz/download | tar -xz
  curl -sL http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/prinseq-lite-0.20.4.tar.gz | tar -xz
  chmod a+x prinseq-lite-0.20.4/*.pl
  mv prinseq-lite-0.20.4 ${TOOLS_PATH}/
  ln -s prinseq-lite-0.20.4 ${TOOLS_PATH}/prinseq

