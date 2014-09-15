##depends:none

# Mafft, BSD
  #cd ${TMPDIR_PATH}/
  #curl http://mafft.cbrc.jp/alignment/software/mafft-7.130-without-extensions-src.tgz | tar -xz
  #cd mafft-7.130-without-extensions/core
  #sed -i 's/PREFIX = \/usr\/local/PREFIX = \/opt\/chipster\/tools\/mafft-7.130-without-extensions/g' Makefile
  #make clean
  #make
  #make install
  curl http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/misc/mafft-7.130-without-extensions-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/
  ln -s mafft-7.130-without-extensions ${TOOLS_PATH}/mafft
