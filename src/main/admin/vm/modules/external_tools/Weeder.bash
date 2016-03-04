##depends:none

# Weeder, custom license, according to developers VM bundling is ok
  cd ${TMPDIR_PATH}/
  #curl -s http://159.149.160.51/modtools/downloads/weeder1.4.2.tar.gz | tar -xz
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/weeder1.4.2.tar.gz | tar -xz
  cd Weeder1.4.2/
  ./compileall
  cd ../
  mkdir ${TOOLS_PATH}/weeder/
  mv Weeder1.4.2/ ${TOOLS_PATH}/weeder/

# Promoter sequence files, license?
  cd ${TMPDIR_PATH}/
  mkdir ${TOOLS_PATH}/weeder/seqs/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/weeder_seqs/All_Weeder_sequence_files_v1.tar.gz | tar -xz -C ${TOOLS_PATH}/weeder/seqs/

