##depends:none

#curl -L http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/R/R-3.3.2-vmbin/R-3.3.2_ubuntu-16.04_2016-11-28.tar.gz | tar -xz -C ${TOOLS_PATH}/
#curl -L http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/R/R-3.3.2-vmbin/R-3.3.2_ubuntu-16.04_2017-06-16.tar.lz4 | lz4 -d | tar -x -C ${TOOLS_PATH}/
lz4 -d /mnt/artefacts/downloads/R-3.3.2_ubuntu-16.04_2017-06-16.tar.lz4 |tar x -C ${TOOLS_PATH}/

# For GATK
R_VER=3.3.2
wget https://raw.githubusercontent.com/broadinstitute/gatk/master/scripts/docker/gatkbase/install_R_packages.R
# Temporary fix because cran.mtu.edu is down
sed -i s_' "http://cran.mtu.edu",'__ install_R_packages.R
${TOOLS_PATH}/R-${R_VER}/bin/Rscript --vanilla install_R_packages.R