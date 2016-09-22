##depends:none

source ../installation_files/functions.bash

# Bowtie 2, Artistic License
  cd ${TMPDIR_PATH}/
  # wget_retry -nv -O bowtie2-2.1.0-linux-x86_64.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.1.0/bowtie2-2.1.0-linux-x86_64.zip/download
  # wget_retry -nv -O bowtie2-2.1.0-linux-x86_64.zip http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie2-2.1.0-linux-x86_64.zip
  wget_retry -nv -O bowtie2-2.2.9-linux-x86_64.zip https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip/download
  unzip -q bowtie2-2.2.9-linux-x86_64.zip
  mv bowtie2-2.2.9 ${TOOLS_PATH}
  ln -s bowtie2-2.2.9 ${TOOLS_PATH}/bowtie2
  rm bowtie2-2.2.9-linux-x86_64.zip
  