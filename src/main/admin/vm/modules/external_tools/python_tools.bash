##depends:none

# Python
  wget_retry -nv https://www.python.org/ftp/python/2.7.12/Python-2.7.12.tgz
  tar xf Python-2.7.12.tgz
  cd Python-2.7.12
  mkdir ${TOOLS_PATH}/Python-2.7.12
  ./configure --prefix ${TOOLS_PATH}/Python-2.7.12
  make
  make install

# install pip
  wget_retry -nv https://bootstrap.pypa.io/get-pip.py
  ${TOOLS_PATH}/Python-2.7.12/bin/python get-pip.py

# Python tools
#
# MACS
  ${TOOLS_PATH}/Python-2.7.12/bin/pip install MACS==1.4.2
  mkdir ${TOOLS_PATH}/macs/
  ln -s ${TOOLS_PATH}/Python-2.7.12/bin/macs14 ${TOOLS_PATH}/macs/macs14
  
# MACS2
  ${TOOLS_PATH}/Python-2.7.12/bin/pip install numpy
  ${TOOLS_PATH}/Python-2.7.12/bin/pip install MACS2==2.1.1.20160309
  ln -s ${TOOLS_PATH}/Python-2.7.12/bin/macs2 ${TOOLS_PATH}/macs/macs2

# RSeQC
  {TOOLS_PATH}/Python-2.7.12/bin/pip install RSeQC==2.6.4

# HTSeq
  ${TOOLS_PATH}/Python-2.7.12/bin/pip install matplotlib
  ${TOOLS_PATH}/Python-2.7.12/bin/pip install HTSeq==0.6.1
  # htseq-count_chr is not part of distribution
  wget_retry -O ${TOOLS_PATH}/Python-2.7.12/bin/htseq-count_chr http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/htseq/htseq-count_chr 
  chmod 755 ${TOOLS_PATH}/Python-2.7.12/bin/htseq-count_chr
  wget_retry -O ${TOOLS_PATH}/Python-2.7.12/lib/python2.7/site-packages/HTSeq/scripts/count_chr.py http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/htseq/count_chr_v2.py
  # make links
  mkdir -p ${TOOLS_PATH}/htseq/
  ln -s ${TOOLS_PATH}/Python-2.7.12/bin/htseq-qa ${TOOLS_PATH}/htseq/htseq-qa
  ln -s ${TOOLS_PATH}/Python-2.7.12/bin/htseq-count ${TOOLS_PATH}/htseq/htseq-count
  ln -s ${TOOLS_PATH}/Python-2.7.12/bin/htseq-count_chr ${TOOLS_PATH}/htseq/htseq-count_chr