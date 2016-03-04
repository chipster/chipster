##depends:none

source ../installation_files/functions.bash

# Bowtie, Artistic License
  cd ${TMPDIR_PATH}/
  #wget_retry -nv -O bowtie-0.12.7-linux-x86_64.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-linux-x86_64.zip/download
  wget_retry -nv -O bowtie-0.12.7-linux-x86_64.zip http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie-0.12.7-linux-x86_64.zip
  unzip -q bowtie-0.12.7-linux-x86_64.zip
  mv bowtie-0.12.7/ ${TOOLS_PATH}
  ln -s bowtie-0.12.7 ${TOOLS_PATH}/bowtie
  rm bowtie-0.12.7-linux-x86_64.zip
  # remove example index
  rm ${TOOLS_PATH}/bowtie/indexes/e_coli.*
