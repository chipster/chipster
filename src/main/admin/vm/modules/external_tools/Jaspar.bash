##depends:external_tools/ClusterBuster.bash

# Jaspar, no license
  cd ${TMPDIR_PATH}/
  wget_retry -nv http://zlab.bu.edu/clover/jaspar2
  mv jaspar2 ${TOOLS_PATH}/ClusterBuster/jaspar2005core.txt
