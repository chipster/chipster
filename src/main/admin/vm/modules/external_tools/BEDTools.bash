##depends:none

# BEDTools, GNU GPL v2
  cd ${TMPDIR_PATH}/
  curl -s http://bedtools.googlecode.com/files/BEDTools.v2.17.0.tar.gz | tar -xz
  cd bedtools-2.17.0
  make clean
  make all
  cd ../
  mv bedtools-2.17.0 ${TOOLS_PATH}/
  ln -s bedtools-2.17.0 ${TOOLS_PATH}/bedtools

  #curl -L http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bedtools-2.17.0-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/

