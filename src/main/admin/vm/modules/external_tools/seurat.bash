##depends:none

# cell cycle markers-file
  cd ${TMPDIR_PATH}/
  mkdir seurat
  cd seurat
  wget http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/seurat/cell_cycle_vignette_files.zip
  unzip cell_cycle_vignette_files.zip
  cd ..
  mv seurat ${TOOLS_PATH}/
