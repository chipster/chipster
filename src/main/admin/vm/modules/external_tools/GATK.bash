##depends:none

# GATK, MIT
#  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/GenomeAnalysisTKLite-latest.tar.bz2 | tar -xj -C ${TOOLS_PATH}/
#  ln -s GenomeAnalysisTKLite-2.1-11-gfb37f33 ${TOOLS_PATH}/GenomeAnalysisTK2

# GATK 3.8
cd ${TMPDIR_PATH}/
wget https://software.broadinstitute.org/gatk/download/auth?package=GATK -O GenomeAnalysisTK-.tar
tar xf GenomeAnalysisTK-.tar -C ${TOOLS_PATH}/
cd ${TOOLS_PATH}
ln -s GenomeAnalysisTK-3.8-0-ge9d806836 GATK

# GATK4
cd ${TMPDIR_PATH}/
wget https://software.broadinstitute.org/gatk/download/auth?package=BETA -O gatk-4.beta.3.zip
unzip gatk-4.beta.3.zip -d ${TOOLS_PATH}
cd ${TOOLS_PATH}
ln -s gatk-4.beta.3-SNAPSHOT GATK4
  