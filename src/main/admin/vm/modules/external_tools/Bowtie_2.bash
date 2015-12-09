##depends:none

source ../installation_files/functions.bash

# Bowtie 2, Artistic License
  cd ${TMPDIR_PATH}/
  # wget_retry -nv -O bowtie2-2.1.0-linux-x86_64.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.1.0/bowtie2-2.1.0-linux-x86_64.zip/download
  wget_retry -nv -O bowtie2-2.1.0-linux-x86_64.zip http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie2-2.1.0-linux-x86_64.zip
  unzip -q bowtie2-2.1.0-linux-x86_64.zip
  mv bowtie2-2.1.0 ${TOOLS_PATH}
  ln -s bowtie2-2.1.0 ${TOOLS_PATH}/bowtie2
  rm bowtie2-2.1.0-linux-x86_64.zip
  
  mkdir -p ${TOOLS_PATH}/bowtie2/indexes/

# GRCh37_74 ensembl transcripts   
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/misc/GRCh37_74_ensembl_transcripts.tar.gz | tar -xzv -C ${TOOLS_PATH}/bowtie2/indexes/
