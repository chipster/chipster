##depends:none

# rtg-tools, BSD License
  cd ${TMPDIR_PATH}/
  wget https://github.com/RealTimeGenomics/rtg-tools/releases/download/3.8.3/rtg-tools-3.8.3-linux-x64.zip
  unzip rtg-tools-3.8.3-linux-x64.zip -d ${TOOLS_PATH}
  cd ${TOOLS_PATH}
  ln -s  rtg-tools-3.8.3 rtg
  