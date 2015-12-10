##depends:none

source ../installation_files/functions.bash

# FastQC, GPL v3 or later
  cd ${TMPDIR_PATH}/
  #wget_retry -nv http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.3.zip 
  wget_retry -nv http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/fastqc_v0.11.3.zip 
  unzip -q fastqc_v0.11.3.zip
  chmod a+x FastQC/fastqc
  mv FastQC/ ${TOOLS_PATH}/
  rm fastqc_v0.11.3.zip
