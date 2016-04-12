##depends:none
# SAM tools, BSD License, MIT License
  cd ${TMPDIR_PATH}/
  #wget_retry https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2
  wget_retry http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/samtools-1.2.tar.bz2
  tar xf samtools-1.2.tar.bz2
  cd samtools-1.2
  perl -pi -e 's/-D_CURSES_LIB=1/-D_CURSES_LIB=0/g' Makefile
  perl -pi -e 's/LIBCURSES=/#LIBCURSES=/g' Makefile
  make
  cd ..
  mv samtools-1.2/ ${TOOLS_PATH}/
  rm samtools-1.2.tar.bz2
  #ln -s ${TOOLS_PATH}/samtools-1.2 ${TOOLS_PATH}/samtools