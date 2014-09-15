##depends:none
# SAM tools, BSD License, MIT License
  cd ${TMPDIR_PATH}/
  #curl -sL http://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2/download | tar -xj
  #cd samtools-0.1.19/
  #make
  #cd ../
  #mv samtools-0.1.19/ ${TOOLS_PATH}
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/samtools-0.1.19-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/
  ln -s samtools-0.1.19 ${TOOLS_PATH}/samtools
