##depends:none

#  Stringtie, OSI Artistic License

  cd ${TMPDIR_PATH}/
  wget_retry -nv -O stringtie-1.3.3b.Linux_x86_64.tar.gz  http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.3b.Linux_x86_64.tar.gz
  tar -xzf stringtie-1.3.3b.Linux_x86_64.tar.gz  -C ${TOOLS_PATH}/
  cd ${TOOLS_PATH}
  ln -s stringtie-1.3.3b.Linux_x86_64 stringtie
  # Clean up
  cd ${TMPDIR_PATH}/
  rm stringtie-1.3.3b.Linux_x86_64.tar.gz
  
