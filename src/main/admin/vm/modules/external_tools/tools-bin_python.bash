##depends:none

# Python
  cd ${TMPDIR_PATH}/
  wget_retry -nv https://www.python.org/ftp/python/2.7.12/Python-2.7.12.tgz
  tar xf Python-2.7.12.tgz
  cd Python-2.7.12
  mkdir ${TOOLS_PATH}/Python-2.7.12
  ./configure --prefix ${TOOLS_PATH}/Python-2.7.12
  make
  make install
  PYTHON=Python-2.7.12

# install pip
  wget_retry -nv https://bootstrap.pypa.io/get-pip.py
  ${TOOLS_PATH}/${PYTHON}/bin/python get-pip.py

# Python tools
#
# S3
  ${TOOLS_PATH}/${PYTHON}/bin/pip install s3cmd

# HTSeq
  ${TOOLS_PATH}/${PYTHON}/bin/pip install matplotlib
  ${TOOLS_PATH}/${PYTHON}/bin/pip install HTSeq==0.6.1
  # htseq-count_chr is not part of distribution
  #wget_retry -O ${TOOLS_PATH}/${PYTHON}/bin/htseq-count_chr http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/htseq/htseq-count_chr 
  cp ${TOOLS_PATH}/${PYTHON}/bin/htseq-count ${TOOLS_PATH}/${PYTHON}/bin/htseq-count_chr
  sed -i 's/HTSeq.scripts.count/HTSeq.scripts.count_chr/' ${TOOLS_PATH}/${PYTHON}/bin/htseq-count_chr
  chmod 755 ${TOOLS_PATH}/${PYTHON}/bin/htseq-count_chr
  wget_retry -O ${TOOLS_PATH}/${PYTHON}/lib/python2.7/site-packages/HTSeq/scripts/count_chr.py http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/htseq/count_chr_v2.py
  # make links
  mkdir -p ${TOOLS_PATH}/htseq/
  cd ${TOOLS_PATH}/htseq/
  ln -s ../${PYTHON}/bin/htseq-qa htseq-qa
  ln -s ../${PYTHON}/bin/htseq-count htseq-count
  ln -s ../${PYTHON}/bin/htseq-count_chr htseq-count_chr
  cd ${TMPDIR_PATH}

# MACS
  ${TOOLS_PATH}/${PYTHON}/bin/pip install MACS==1.4.2
  mkdir ${TOOLS_PATH}/macs/
  cd ${TOOLS_PATH}/macs/
  ln -s ../${PYTHON}/bin/macs14 macs14
  cd ${TMPDIR_PATH}
  
# MACS2
  ${TOOLS_PATH}/${PYTHON}/bin/pip install numpy
  ${TOOLS_PATH}/${PYTHON}/bin/pip install MACS2==2.1.1.20160309
  cd ${TOOLS_PATH}/macs/
  ln -s ../${PYTHON}/bin/macs2 macs2
  cd ${TMPDIR_PATH}

# MultiQC
  ${TOOLS_PATH}/${PYTHON}/bin/pip install multiqc
  mkdir ${TOOLS_PATH}/multiqc
  cd ${TOOLS_PATH}/multiqc
  ln -s ../${PYTHON}/bin multiqc
  cd ${TMPDIR_PATH}
  
# RSeQC
  ${TOOLS_PATH}/${PYTHON}/bin/pip install RSeQC==2.6.4
  mkdir ${TOOLS_PATH}/rseqc
  cd ${TOOLS_PATH}/rseqc
  ln -s ../${PYTHON}/bin rseqc
  cd ${TMPDIR_PATH}

# ZIFA
# Requires: scipy, pandas, numpy,  scikits.learn, matplotlib, pandas
# Some of these modules are already installed, but lets still list all the requirements for ZIFA
${TOOLS_PATH}/${PYTHON}/bin/pip install numpy scipy scikits.learn matplotlib
# pandas >0.21 needs currently some extra care
${TOOLS_PATH}/${PYTHON}/bin/pip install --no-build-isolation pandas